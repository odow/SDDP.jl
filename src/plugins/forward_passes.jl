#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    DefaultForwardPass(; include_last_node::Bool = true)

The default forward pass.

If `include_last_node = false` and the sample terminated due to a cycle, then
the last node (which forms the cycle) is omitted. This can be useful option to
set when training, but it comes at the cost of not knowing which node formed the
cycle (if there are multiple possibilities).
"""
struct DefaultForwardPass <: AbstractForwardPass
    include_last_node::Bool
    function DefaultForwardPass(; include_last_node::Bool = true)
        return new(include_last_node)
    end
end

function forward_pass(
    model::PolicyGraph{T},
    options::Options,
    pass::DefaultForwardPass,
) where {T}
    # First up, sample a scenario. Note that if a cycle is detected, this will
    # return the cycle node as well.
    @_timeit_threadsafe model.timer_output "sample_scenario" begin
        scenario_path, terminated_due_to_cycle =
            sample_scenario(model, options.sampling_scheme)
    end
    final_node = scenario_path[end]
    if terminated_due_to_cycle && !pass.include_last_node
        pop!(scenario_path)
    end
    # Storage for the list of outgoing states that we visit on the forward pass.
    sampled_states = Dict{Symbol,Float64}[]
    # Storage for the belief states: partition index and the belief dictionary.
    belief_states = Tuple{Int,Dict{T,Float64}}[]
    current_belief = initialize_belief(model)
    # Our initial incoming state.
    incoming_state_value = copy(options.initial_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0
    # Objective state interpolation.
    objective_state_vector, N =
        initialize_objective_state(model[scenario_path[1][1]])
    objective_states = NTuple{N,Float64}[]
    # Iterate down the scenario.
    for (depth, (node_index, noise)) in enumerate(scenario_path)
        node = model[node_index]
        lock(node.lock)
        try
            # Objective state interpolation.
            objective_state_vector = update_objective_state(
                node.objective_state,
                objective_state_vector,
                noise,
            )
            if objective_state_vector !== nothing
                push!(objective_states, objective_state_vector)
            end
            # Update belief state, etc.
            if node.belief_state !== nothing
                belief = node.belief_state::BeliefState{T}
                partition_index = belief.partition_index
                current_belief = belief.updater(
                    belief.belief,
                    current_belief,
                    partition_index,
                    noise,
                )
                push!(belief_states, (partition_index, copy(current_belief)))
            end
            # ===== Begin: starting state for infinite horizon =====
            starting_states = options.starting_states[node_index]
            if length(starting_states) > 0
                # There is at least one other possible starting state. If our
                # incoming state is more than δ away from the other states, add it
                # as a possible starting state.
                if distance(starting_states, incoming_state_value) >
                   options.cycle_discretization_delta
                    push!(starting_states, incoming_state_value)
                end
                # TODO(odow):
                # - A better way of randomly sampling a starting state.
                # - Is is bad that we splice! here instead of just sampling? For
                #   convergence it is probably bad, since our list of possible
                #   starting states keeps changing, but from a computational
                #   perspective, we don't want to keep a list of discretized points
                #   in the state-space δ distance apart...
                incoming_state_value =
                    splice!(starting_states, rand(1:length(starting_states)))
            end
            # ===== End: starting state for infinite horizon =====
            # Solve the subproblem, note that `duality_handler = nothing`.
            @_timeit_threadsafe model.timer_output "solve_subproblem" begin
                subproblem_results = solve_subproblem(
                    model,
                    node,
                    incoming_state_value,
                    noise,
                    scenario_path[1:depth];
                    duality_handler = nothing,
                )
            end
            # Cumulate the stage_objective.
            cumulative_value += subproblem_results.stage_objective
            # Set the outgoing state value as the incoming state value for the next
            # node.
            incoming_state_value = copy(subproblem_results.state)
            # Add the outgoing state variable to the list of states we have sampled
            # on this forward pass.
            push!(sampled_states, incoming_state_value)
        finally
            unlock(node.lock)
        end
    end
    if terminated_due_to_cycle
        # We terminated due to a cycle. Here is the list of possible
        # starting states for that node:
        starting_states = options.starting_states[final_node[1]]
        # We also need the incoming state variable to the final node, which
        # is the outgoing state value of the second to last node:
        incoming_state_value = if pass.include_last_node
            sampled_states[end-1]
        else
            sampled_states[end]
        end
        # If this incoming state value is more than δ away from another
        # state, add it to the list.
        if distance(starting_states, incoming_state_value) >
           options.cycle_discretization_delta
            push!(starting_states, incoming_state_value)
        end
    end
    # ===== End: drop off starting state if terminated due to cycle =====
    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        objective_states = objective_states,
        belief_states = belief_states,
        cumulative_value = cumulative_value,
    )
end

mutable struct RevisitingForwardPass <: AbstractForwardPass
    period::Int
    sub_pass::AbstractForwardPass
    archive::Vector{Any}
    last_index::Int
    counter::Int
end

"""
    RevisitingForwardPass(
        period::Int = 500;
        sub_pass::AbstractForwardPass = DefaultForwardPass(),
    )

A forward pass scheme that generate `period` new forward passes (using
`sub_pass`), then revisits all previously explored forward passes. This can
be useful to encourage convergence at a diversity of points in the
state-space.

Set `period = typemax(Int)` to disable.

For example, if `period = 2`, then the forward passes will be revisited as
follows: `1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 1, 2, ...`.
"""
function RevisitingForwardPass(
    period::Int = 500;
    sub_pass::AbstractForwardPass = DefaultForwardPass(),
)
    @assert period > 0
    return RevisitingForwardPass(period, sub_pass, Any[], 0, 0)
end

function forward_pass(
    model::PolicyGraph,
    options::Options,
    fp::RevisitingForwardPass,
)
    fp.counter += 1
    if fp.counter - fp.period > fp.last_index
        fp.counter = 1
        fp.last_index = length(fp.archive)
    end
    if fp.counter <= length(fp.archive)
        return fp.archive[fp.counter]
    else
        pass = forward_pass(model, options, fp.sub_pass)
        push!(fp.archive, pass)
        return pass
    end
end

mutable struct RiskAdjustedForwardPass{F,T} <: AbstractForwardPass
    forward_pass::F
    risk_measure::T
    resampling_probability::Float64
    rejection_count::Int
    objectives::Vector{Float64}
    nominal_probability::Vector{Float64}
    adjusted_probability::Vector{Float64}
    archive::Vector{Any}
    resample_count::Vector{Int}
end

"""
    RiskAdjustedForwardPass(;
        forward_pass::AbstractForwardPass,
        risk_measure::AbstractRiskMeasure,
        resampling_probability::Float64,
        rejection_count::Int = 5,
    )

A forward pass that resamples a previous forward pass with
`resampling_probability` probability, and otherwise samples a new forward pass
using `forward_pass`.

The forward pass to revisit is chosen based on the risk-adjusted (using
`risk_measure`) probability of the cumulative stage objectives.

Note that this objective corresponds to the _first_ time we visited the
trajectory. Subsequent visits may have improved things, but we don't have the
mechanisms in-place to update it. Therefore, remove the forward pass from
resampling consideration after `rejection_count` revisits.
"""
function RiskAdjustedForwardPass(;
    forward_pass::AbstractForwardPass,
    risk_measure::AbstractRiskMeasure,
    resampling_probability::Float64,
    rejection_count::Int = 5,
)
    if !(0 < resampling_probability < 1)
        throw(ArgumentError("Resampling probability must be in `(0, 1)`"))
    end
    return RiskAdjustedForwardPass{typeof(forward_pass),typeof(risk_measure)}(
        forward_pass,
        risk_measure,
        resampling_probability,
        rejection_count,
        Float64[],
        Float64[],
        Float64[],
        Any[],
        Int[],
    )
end

function forward_pass(
    model::PolicyGraph,
    options::Options,
    fp::RiskAdjustedForwardPass,
)
    if length(fp.archive) > 0 && rand() < fp.resampling_probability
        r = rand()
        for i in 1:length(fp.adjusted_probability)
            r -= fp.adjusted_probability[i]
            if r > 1e-8
                continue
            end
            pass = fp.archive[i]
            if fp.resample_count[i] >= fp.rejection_count
                # We've explored this pass too many times. Kick it out of the
                # archive.
                splice!(fp.objectives, i)
                splice!(fp.nominal_probability, i)
                splice!(fp.adjusted_probability, i)
                splice!(fp.archive, i)
                splice!(fp.resample_count, i)
            else
                fp.resample_count[i] += 1
            end
            return pass
        end
    end
    pass = forward_pass(model, options, fp.forward_pass)
    push!(fp.objectives, pass.cumulative_value)
    push!(fp.nominal_probability, 0.0)
    fill!(fp.nominal_probability, 1 / length(fp.nominal_probability))
    push!(fp.adjusted_probability, 0.0)
    push!(fp.archive, pass)
    push!(fp.resample_count, 1)
    adjust_probability(
        fp.risk_measure,
        fp.adjusted_probability,
        fp.nominal_probability,
        fp.objectives,
        fp.objectives,
        model.objective_sense == MOI.MIN_SENSE,
    )
    return pass
end

"""
    RegularizedForwardPass(;
        rho::Float64 = 0.05,
        forward_pass::AbstractForwardPass = DefaultForwardPass(),
    )

A forward pass that regularizes the outgoing first-stage state variables with an
L-infty trust-region constraint about the previous iteration's solution.
Specifically, the bounds of the outgoing state variable `x` are updated from
`(l, u)` to `max(l, x^k - rho * (u - l)) <= x <= min(u, x^k + rho * (u - l))`,
where `x^k` is the optimal solution of `x` in the previous iteration. On the
first iteration, the value of the state at the root node is used.

By default, `rho` is set to 5%, which seems to work well empirically.

Pass a different `forward_pass` to control the forward pass within the
regularized forward pass.

This forward pass is largely intended to be used for investment problems in
which the first stage makes a series of capacity decisions that then influence
the rest of the graph. An error is thrown if the first stage problem is not
deterministic, and states are silently skipped if they do not have finite
bounds.
"""
mutable struct RegularizedForwardPass{T<:AbstractForwardPass} <:
               AbstractForwardPass
    forward_pass::T
    trial_centre::Dict{Symbol,Float64}
    old_bounds::Dict{Symbol,Tuple{Float64,Float64}}
    ρ::Float64

    function RegularizedForwardPass(;
        rho::Float64 = 0.05,
        forward_pass::AbstractForwardPass = DefaultForwardPass(),
    )
        centre = Dict{Symbol,Float64}()
        old_bounds = Dict{Symbol,Tuple{Float64,Float64}}()
        return new{typeof(forward_pass)}(forward_pass, centre, old_bounds, rho)
    end
end

function forward_pass(
    model::PolicyGraph,
    options::Options,
    fp::RegularizedForwardPass,
)
    if length(model.root_children) != 1
        error(
            "RegularizedForwardPass cannot be applied because first-stage is " *
            "not deterministic",
        )
    end
    node = model[model.root_children[1].term]
    if length(node.noise_terms) > 1
        error(
            "RegularizedForwardPass cannot be applied because first-stage is " *
            "not deterministic",
        )
    end
    # It's safe to lock this node while we modify the forward pass and the node
    # because there is only one node in the first stage and it is deterministic.
    lock(node.lock) do
        for (k, v) in node.states
            if !(has_lower_bound(v.out) && has_upper_bound(v.out))
                continue  # Not a finitely bounded state. Ignore for now
            end
            if !haskey(fp.old_bounds, k)
                fp.old_bounds[k] = (lower_bound(v.out), upper_bound(v.out))
            end
            l, u = fp.old_bounds[k]
            x = get(fp.trial_centre, k, model.initial_root_state[k])
            set_lower_bound(v.out, max(l, x - fp.ρ * (u - l)))
            set_upper_bound(v.out, min(u, x + fp.ρ * (u - l)))
        end
        return
    end
    pass = forward_pass(model, options, fp.forward_pass)
    # We're locking the node to reset the variable bounds back to their default.
    # There are some potential scheduling issues to be aware of:
    #
    #  * Thread A might have obtained the lock, modified the bounds to the trial
    #    centre, released the lock, and then entered `forward_pass`
    #  * But before it can re-obtain a lock on the first node to solve the
    #    problem, Thread B has come along below and modified the bounds back to
    #    their original value.
    #  * This means that Thread A is solving the unregularized problem, but it
    #    doesn't really matter because it doesn't change the validity; this is a
    #    performance optimization.
    #  * It might also be that we "skip" some of the starting trial points,
    #    because Thread A sets a trial centre, then thread B comes along and
    #    sets a new one before A can start the forward pass. Again, this doesn't
    #    really matter; this regularization is just a performance optimization.
    lock(node.lock) do
        for (k, (l, u)) in fp.old_bounds
            fp.trial_centre[k] = pass.sampled_states[1][k]
            set_lower_bound(node.states[k].out, l)
            set_upper_bound(node.states[k].out, u)
        end
        return
    end
    return pass
end

struct ImportanceSamplingForwardPass <: AbstractForwardPass end

function forward_pass(
    model::PolicyGraph{T},
    options::Options,
    pass::ImportanceSamplingForwardPass,
) where {T}
    @assert isempty(model.belief_partition)
    scenario_path = Tuple{T,Any}[]
    sampled_states = Dict{Symbol,Float64}[]
    cumulative_value = 0.0
    incoming_state_value = copy(options.initial_state)
    node_index = sample_noise(model.root_children)
    node = model[node_index]
    noise = sample_noise(node.noise_terms)
    while node_index !== nothing
        node = model[node_index]
        lock(node.lock)
        try
            push!(scenario_path, (node_index, noise))
            subproblem_results = solve_subproblem(
                model,
                node,
                incoming_state_value,
                noise,
                scenario_path;
                duality_handler = nothing,
            )
            cumulative_value += subproblem_results.stage_objective
            incoming_state_value = copy(subproblem_results.state)
            push!(sampled_states, incoming_state_value)
            if isempty(node.bellman_function.local_thetas)
                # First iteration with no multi-cuts, or a node with no children
                node_index = SDDP.sample_noise(node.children)
                if node_index !== nothing
                    new_node = model[node_index]
                    noise = SDDP.sample_noise(new_node.noise_terms)
                end
            else
                objectives = map(node.bellman_function.local_thetas) do t
                    return JuMP.value(t.theta)
                end
                adjusted_probability = fill(NaN, length(objectives))
                nominal_probability = Float64[]
                support = Any[]
                for child in node.children
                    for noise in model[child.term].noise_terms
                        push!(nominal_probability, child.probability * noise.probability)
                        push!(support, (child.term, noise.term))
                    end
                end
                @assert length(nominal_probability) == length(objectives)
                _ = SDDP.adjust_probability(
                    options.risk_measures[node_index],
                    adjusted_probability,
                    nominal_probability,
                    support,
                    objectives,
                    model.objective_sense == MOI.MIN_SENSE,
                )
                terms = SDDP.Noise.(support, adjusted_probability)
                node_index, noise = SDDP.sample_noise(terms)
            end
        finally
            unlock(node.lock)
        end
    end
    return (
        scenario_path = scenario_path,
        sampled_states = sampled_states,
        objective_states = NTuple{0,Float64}[],
        belief_states = Tuple{Int,Dict{T,Float64}}[],
        cumulative_value = cumulative_value,
    )
end
