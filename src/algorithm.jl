#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

const SDDP_TIMER = TimerOutputs.TimerOutput()

# to_nodal_form is an internal helper function so users can pass arguments like:
# risk_measure = SDDP.Expectation(),
# risk_measure = Dict(1=>Expectation(), 2=>WorstCase())
# risk_measure = (node_index) -> node_index == 1 ? Expectation() : WorstCase()
# It will return a dictionary with a key for each node_index in the policy
# graph, and a corresponding value of whatever the user provided.
function to_nodal_form(model::PolicyGraph{T}, element) where {T}
    # Note: we don't copy element here, so if element is mutable, you should use
    # to_nodal_form(model, x -> new_element()) instead. A good example is
    # Vector{T}; use to_nodal_form(model, i -> T[]).
    store = Dict{T,typeof(element)}()
    for node_index in keys(model.nodes)
        store[node_index] = element
    end
    return store
end

function to_nodal_form(model::PolicyGraph{T}, builder::Function) where {T}
    store = Dict{T,Any}()
    for node_index in keys(model.nodes)
        store[node_index] = builder(node_index)
    end
    V = typeof(first(values(store)))
    for val in values(store)
        V = promote_type(V, typeof(val))
    end
    return Dict{T,V}(key => val for (key, val) in store)
end

function to_nodal_form(model::PolicyGraph{T}, dict::Dict{T,V}) where {T,V}
    for key in keys(model.nodes)
        if !haskey(dict, key)
            error("Missing key: $(key).")
        end
    end
    return dict
end

# Internal function: returns a dictionary with a key for each node, where the
# value is a list of other nodes that contain the same children. This is useful
# because on the backward pass we can add cuts to nodes with the same children
# without having to re-solve the children.
function get_same_children(model::PolicyGraph{T}) where {T}
    tmp = Dict{Set{T},Set{T}}()
    for (key, node) in model.nodes
        children = Set(child.term for child in node.children)
        if length(children) == 0
            continue
        elseif haskey(tmp, children)
            push!(tmp[children], key)
        else
            tmp[children] = Set{T}([key])
        end
    end
    same_children = Dict{T,Vector{T}}(key => T[] for key in keys(model.nodes))
    for set in values(tmp)
        for v in set
            same_children[v] = collect(setdiff(set, Ref(v)))
        end
    end
    return same_children
end

# Internal struct: storage for SDDP options and cached data. Users shouldn't
# interact with this directly.
struct Options{T}
    # The initial state to start from the root node.
    initial_state::Dict{Symbol,Float64}
    # The sampling scheme to use on the forward pass.
    sampling_scheme::AbstractSamplingScheme
    backward_sampling_scheme::AbstractBackwardSamplingScheme
    # Storage for the set of possible sampling states at each node. We only use
    # this if there is a cycle in the policy graph.
    starting_states::Dict{T,Vector{Dict{Symbol,Float64}}}
    # Risk measure to use at each node.
    risk_measures::Dict{T,AbstractRiskMeasure}
    # The delta by which to check if a state is close to a previously sampled
    # state.
    cycle_discretization_delta::Float64
    # Flag to add cuts to similar nodes.
    refine_at_similar_nodes::Bool
    # The node transition matrix.
    Φ::Dict{Tuple{T,T},Float64}
    # A list of nodes that contain a subset of the children of node i.
    similar_children::Dict{T,Vector{T}}
    stopping_rules::Vector{AbstractStoppingRule}
    dashboard_callback::Function
    print_level::Int
    start_time::Float64
    log::Vector{Log}
    log_file_handle::Any
    log_frequency::Int
    forward_pass::AbstractForwardPass
    duality_handler::AbstractDualityHandler
    # A callback called after the forward pass.
    forward_pass_callback::Any

    # Internal function: users should never construct this themselves.
    function Options(
        model::PolicyGraph{T},
        initial_state::Dict{Symbol,Float64},
        sampling_scheme::AbstractSamplingScheme,
        backward_sampling_scheme::AbstractBackwardSamplingScheme,
        risk_measures,
        cycle_discretization_delta::Float64,
        refine_at_similar_nodes::Bool,
        stopping_rules::Vector{AbstractStoppingRule},
        dashboard_callback::Function,
        print_level::Int,
        start_time::Float64,
        log::Vector{Log},
        log_file_handle,
        log_frequency::Int,
        forward_pass::AbstractForwardPass,
        duality_handler::AbstractDualityHandler,
        forward_pass_callback,
    ) where {T}
        return new{T}(
            initial_state,
            sampling_scheme,
            backward_sampling_scheme,
            to_nodal_form(model, x -> Dict{Symbol,Float64}[]),
            to_nodal_form(model, risk_measures),
            cycle_discretization_delta,
            refine_at_similar_nodes,
            build_Φ(model),
            get_same_children(model),
            stopping_rules,
            dashboard_callback,
            print_level,
            start_time,
            log,
            log_file_handle,
            log_frequency,
            forward_pass,
            duality_handler,
            forward_pass_callback,
        )
    end
end

# Internal function: set the incoming state variables of node to the values
# contained in state.
function set_incoming_state(node::Node, state::Dict{Symbol,Float64})
    for (state_name, value) in state
        JuMP.fix(node.states[state_name].in, value)
    end
    return
end

# Internal function: get the values of the outgoing state variables in node.
# Requires node.subproblem to have been solved with PrimalStatus ==
# FeasiblePoint.
function get_outgoing_state(node::Node)
    values = Dict{Symbol,Float64}()
    for (name, state) in node.states
        # To fix some cases of numerical infeasiblities, if the outgoing value
        # is outside its bounds, project the value back onto the bounds. There
        # is a pretty large (×5) penalty associated with this check because it
        # typically requires a call to the solver. It is worth reducing
        # infeasibilities though.
        outgoing_value = JuMP.value(state.out)
        if JuMP.has_upper_bound(state.out)
            current_bound = JuMP.upper_bound(state.out)
            if current_bound < outgoing_value
                outgoing_value = current_bound
            end
        end
        if JuMP.has_lower_bound(state.out)
            current_bound = JuMP.lower_bound(state.out)
            if current_bound > outgoing_value
                outgoing_value = current_bound
            end
        end
        values[name] = outgoing_value
    end
    return values
end

# Internal function: set the objective of node to the stage objective, plus the
# cost/value-to-go term.
function set_objective(node::Node{T}) where {T}
    objective_state_component = get_objective_state_component(node)
    belief_state_component = get_belief_state_component(node)
    if objective_state_component != JuMP.AffExpr(0.0) ||
       belief_state_component != JuMP.AffExpr(0.0)
        node.stage_objective_set = false
    end
    if !node.stage_objective_set
        JuMP.set_objective(
            node.subproblem,
            JuMP.objective_sense(node.subproblem),
            @expression(
                node.subproblem,
                node.stage_objective +
                objective_state_component +
                belief_state_component +
                bellman_term(node.bellman_function)
            )
        )
    end
    node.stage_objective_set = true
    return
end

# Internal function: overload for the case where JuMP.value fails on a
# Real number.
stage_objective_value(stage_objective::Real) = stage_objective
stage_objective_value(stage_objective) = JuMP.value(stage_objective)

"""
    write_subproblem_to_file(
        node::Node,
        filename::String;
        throw_error::Bool = false,
    )

Write the subproblem contained in `node` to the file `filename`.
"""
function write_subproblem_to_file(
    node::Node,
    filename::String;
    throw_error::Bool = false,
)
    model = MOI.FileFormats.Model(filename = filename)
    MOI.copy_to(model, JuMP.backend(node.subproblem))
    MOI.write_to_file(model, filename)
    if throw_error
        error(
            "Unable to retrieve solution from $(node.index).\n",
            "  Termination status: $(JuMP.termination_status(node.subproblem))\n",
            "  Primal status:      $(JuMP.primal_status(node.subproblem))\n",
            "  Dual status:        $(JuMP.dual_status(node.subproblem)).\n",
            "A MathOptFormat file was written to `$(filename)`.\n",
            "See https://odow.github.io/SDDP.jl/latest/tutorial/06_warnings/#Numerical-stability-1",
            "\nfor more information.",
        )
    end
    return
end

"""
    parameterize(node::Node, noise)

Parameterize node `node` with the noise `noise`.
"""
function parameterize(node::Node, noise)
    node.parameterize(noise)
    set_objective(node)
    return
end

function attempt_numerical_recovery(node::Node)
    @warn("Attempting to recover from serious numerical issues...")
    if JuMP.mode(node.subproblem) == JuMP.DIRECT
        @warn(
            "Unable to recover in direct mode! Remove `direct = true` when " *
            "creating the policy graph."
        )
    else
        MOI.Utilities.reset_optimizer(node.subproblem)
        optimize!(node.subproblem)
    end
    if JuMP.primal_status(node.subproblem) != JuMP.MOI.FEASIBLE_POINT
        write_subproblem_to_file(
            node,
            "subproblem_$(node.index).mof.json",
            throw_error = true,
        )
    end
    return
end

"""
    _initialize_solver(node::Node; throw_error::Bool)

After passing a model to a different process, we need to set the optimizer
again.

If `throw_error`, throw an error if the model is in direct mode.

See also: [`_uninitialize_solver`](@ref).
"""
function _initialize_solver(node::Node; throw_error::Bool)
    if mode(node.subproblem) == DIRECT
        if throw_error
            error(
                "Cannot use asynchronous solver with optimizers in direct mode.",
            )
        end
    elseif MOI.Utilities.state(backend(node.subproblem)) == MOIU.NO_OPTIMIZER
        if node.optimizer === nothing
            error(
                """
          You must supply an optimizer for the policy graph, either by passing
          one to the `optimizer` keyword argument to `PolicyGraph`, or by
          using `JuMP.set_optimizer(model, optimizer)`.
          """,
            )
        end
        set_optimizer(node.subproblem, node.optimizer)
    end
    return
end

"""
    _initialize_solver(model::PolicyGraph; throw_error::Bool)

After passing a model to a different process, we need to set the optimizer
again.

If `throw_error`, throw an error if the model is in direct mode.

See also: [`_uninitialize_solver`](@ref).
"""
function _initialize_solver(model::PolicyGraph; throw_error::Bool)
    for (_, node) in model.nodes
        _initialize_solver(node; throw_error = throw_error)
    end
    return
end

"""
    _uninitialize_solver(model; throw_error::Bool)

Before passing a model to a different process, we need to drop the inner solver
in case it has some C pointers that we cannot serialize (e.g., GLPK).

If `throw_error`, throw an error if the model is in direct mode.

See also: [`_initialize_solver`](@ref).
"""
function _uninitialize_solver(model::PolicyGraph; throw_error::Bool)
    for (_, node) in model.nodes
        if mode(node.subproblem) == DIRECT
            if throw_error
                error(
                    "Cannot use asynchronous solver with optimizers in direct mode.",
                )
            end
        elseif MOI.Utilities.state(backend(node.subproblem)) !=
               MOIU.NO_OPTIMIZER
            MOI.Utilities.drop_optimizer(node.subproblem)
        end
    end
    return
end

# Internal function: solve the subproblem associated with node given the
# incoming state variables state and realization of the stagewise-independent
# noise term noise.
function solve_subproblem(
    model::PolicyGraph{T},
    node::Node{T},
    state::Dict{Symbol,Float64},
    noise,
    scenario_path::Vector{Tuple{T,S}};
    duality_handler::Union{Nothing,AbstractDualityHandler},
) where {T,S}
    _initialize_solver(node; throw_error = false)
    # Parameterize the model. First, fix the value of the incoming state
    # variables. Then parameterize the model depending on `noise`. Finally,
    # set the objective.
    set_incoming_state(node, state)
    parameterize(node, noise)
    pre_optimize_ret = if node.pre_optimize_hook !== nothing
        node.pre_optimize_hook(
            model,
            node,
            state,
            noise,
            scenario_path,
            duality_handler,
        )
    else
        nothing
    end
    JuMP.optimize!(node.subproblem)
    if haskey(model.ext, :total_solves)
        model.ext[:total_solves] += 1
    else
        model.ext[:total_solves] = 1
    end
    if JuMP.primal_status(node.subproblem) != JuMP.MOI.FEASIBLE_POINT
        attempt_numerical_recovery(node)
    end
    state = get_outgoing_state(node)
    stage_objective = stage_objective_value(node.stage_objective)
    TimerOutputs.@timeit SDDP_TIMER "get_dual_solution" begin
        objective, dual_values = get_dual_solution(node, duality_handler)
    end
    if node.post_optimize_hook !== nothing
        node.post_optimize_hook(pre_optimize_ret)
    end
    return (
        state = state,
        duals = dual_values,
        objective = objective,
        stage_objective = stage_objective,
    )
end

# Internal function to get the objective state at the root node.
function initialize_objective_state(first_node::Node)
    objective_state = first_node.objective_state
    if objective_state !== nothing
        initial_objective_state = objective_state.initial_value
        return initial_objective_state, length(initial_objective_state)
    else
        return nothing, 0
    end
end

# Internal function: update the objective state given incoming `current_state`
# and `noise`.
update_objective_state(::Nothing, ::Any, ::Any) = nothing

function update_objective_state(obj_state, current_state, noise)
    if length(current_state) == 1
        obj_state.state = (obj_state.update(current_state[1], noise),)
    else
        obj_state.state = obj_state.update(current_state, noise)
    end
    return obj_state.state
end

# Internal function: calculate the initial belief state.
function initialize_belief(model::PolicyGraph{T}) where {T}
    current_belief = Dict{T,Float64}(keys(model.nodes) .=> 0.0)
    current_belief[model.root_node] = 1.0
    return current_belief
end

# Internal function: calculate the minimum distance between the state `state`
# and the list of states in `starting_states` using the distance measure `norm`.
function distance(
    starting_states::Vector{Dict{Symbol,Float64}},
    state::Dict{Symbol,Float64},
    norm::Function = inf_norm,
)
    if length(starting_states) == 0
        return Inf
    end
    return minimum(norm.(starting_states, Ref(state)))
end

# Internal function: the norm to use when checking the distance between two
# possible starting states. We're going to use: d(x, y) = |x - y| / (1 + |y|).
function inf_norm(x::Dict{Symbol,Float64}, y::Dict{Symbol,Float64})
    norm = 0.0
    for (key, value) in y
        if abs(x[key] - value) > norm
            norm = abs(x[key] - value) / (1 + abs(value))
        end
    end
    return norm
end

# Internal function: perform a backward pass of the SDDP algorithm along the
# scenario_path, refining the bellman function at sampled_states. Assumes that
# scenario_path does not end in a leaf node (i.e., the forward pass was solved
# with include_last_node = false)
function backward_pass(
    model::PolicyGraph{T},
    options::Options,
    scenario_path::Vector{Tuple{T,NoiseType}},
    sampled_states::Vector{Dict{Symbol,Float64}},
    objective_states::Vector{NTuple{N,Float64}},
    belief_states::Vector{Tuple{Int,Dict{T,Float64}}},
) where {T,NoiseType,N}
    TimerOutputs.@timeit SDDP_TIMER "prepare_backward_pass" begin
        restore_duality = prepare_backward_pass(
            model,
            options.duality_handler,
            options,
        )
    end
    # TODO(odow): improve storage type.
    cuts = Dict{T,Vector{Any}}(index => Any[] for index in keys(model.nodes))
    for index in length(scenario_path):-1:1
        outgoing_state = sampled_states[index]
        objective_state = get(objective_states, index, nothing)
        partition_index, belief_state = get(belief_states, index, (0, nothing))
        items = BackwardPassItems(T, Noise)
        if belief_state !== nothing
            # Update the cost-to-go function for partially observable model.
            for (node_index, belief) in belief_state
                belief == 0.0 && continue
                solve_all_children(
                    model,
                    model[node_index],
                    items,
                    belief,
                    belief_state,
                    objective_state,
                    outgoing_state,
                    options.backward_sampling_scheme,
                    scenario_path[1:index],
                    options.duality_handler,
                )
            end
            # We need to refine our estimate at all nodes in the partition.
            for node_index in model.belief_partition[partition_index]
                node = model[node_index]
                # Update belief state, etc.
                current_belief = node.belief_state::BeliefState{T}
                for (idx, belief) in belief_state
                    current_belief.belief[idx] = belief
                end
                new_cuts = refine_bellman_function(
                    model,
                    node,
                    node.bellman_function,
                    options.risk_measures[node_index],
                    outgoing_state,
                    items.duals,
                    items.supports,
                    items.probability .* items.belief,
                    items.objectives,
                )
                push!(cuts[node_index], new_cuts)
            end
        else
            node_index, _ = scenario_path[index]
            node = model[node_index]
            if length(node.children) == 0
                continue
            end
            solve_all_children(
                model,
                node,
                items,
                1.0,
                belief_state,
                objective_state,
                outgoing_state,
                options.backward_sampling_scheme,
                scenario_path[1:index],
                options.duality_handler,
            )
            new_cuts = refine_bellman_function(
                model,
                node,
                node.bellman_function,
                options.risk_measures[node_index],
                outgoing_state,
                items.duals,
                items.supports,
                items.probability,
                items.objectives,
            )
            push!(cuts[node_index], new_cuts)
            if options.refine_at_similar_nodes
                # Refine the bellman function at other nodes with the same
                # children, e.g., in the same stage of a Markovian policy graph.
                for other_index in options.similar_children[node_index]
                    copied_probability = similar(items.probability)
                    other_node = model[other_index]
                    for (idx, child_index) in enumerate(items.nodes)
                        copied_probability[idx] =
                            get(options.Φ, (other_index, child_index), 0.0) *
                            items.supports[idx].probability
                    end
                    new_cuts = refine_bellman_function(
                        model,
                        other_node,
                        other_node.bellman_function,
                        options.risk_measures[other_index],
                        outgoing_state,
                        items.duals,
                        items.supports,
                        copied_probability,
                        items.objectives,
                    )
                    push!(cuts[other_index], new_cuts)
                end
            end
        end
    end
    TimerOutputs.@timeit SDDP_TIMER "prepare_backward_pass" begin
        restore_duality()
    end
    return cuts
end

struct BackwardPassItems{T,U}
    "Given a (node, noise) tuple, index the element in the array."
    cached_solutions::Dict{Tuple{T,Any},Int}
    duals::Vector{Dict{Symbol,Float64}}
    supports::Vector{U}
    nodes::Vector{T}
    probability::Vector{Float64}
    objectives::Vector{Float64}
    belief::Vector{Float64}
    function BackwardPassItems(T, U)
        return new{T,U}(
            Dict{Tuple{T,Any},Int}(),
            Dict{Symbol,Float64}[],
            U[],
            T[],
            Float64[],
            Float64[],
            Float64[],
        )
    end
end

function solve_all_children(
    model::PolicyGraph{T},
    node::Node{T},
    items::BackwardPassItems,
    belief::Float64,
    belief_state,
    objective_state,
    outgoing_state::Dict{Symbol,Float64},
    backward_sampling_scheme::AbstractBackwardSamplingScheme,
    scenario_path,
    duality_handler::Union{Nothing,AbstractDualityHandler},
) where {T}
    length_scenario_path = length(scenario_path)
    for child in node.children
        if isapprox(child.probability, 0.0, atol = 1e-6)
            continue
        end
        child_node = model[child.term]
        for noise in
            sample_backward_noise_terms(backward_sampling_scheme, child_node)
            if length(scenario_path) == length_scenario_path
                push!(scenario_path, (child.term, noise.term))
            else
                scenario_path[end] = (child.term, noise.term)
            end
            if haskey(items.cached_solutions, (child.term, noise.term))
                sol_index = items.cached_solutions[(child.term, noise.term)]
                push!(items.duals, items.duals[sol_index])
                push!(items.supports, items.supports[sol_index])
                push!(items.nodes, child_node.index)
                push!(items.probability, items.probability[sol_index])
                push!(items.objectives, items.objectives[sol_index])
                push!(items.belief, belief)
            else
                # Update belief state, etc.
                if belief_state !== nothing
                    current_belief = child_node.belief_state::BeliefState{T}
                    current_belief.updater(
                        current_belief.belief,
                        belief_state,
                        current_belief.partition_index,
                        noise.term,
                    )
                end
                if objective_state !== nothing
                    update_objective_state(
                        child_node.objective_state,
                        objective_state,
                        noise.term,
                    )
                end
                TimerOutputs.@timeit SDDP_TIMER "solve_subproblem" begin
                    subproblem_results = solve_subproblem(
                        model,
                        child_node,
                        outgoing_state,
                        noise.term,
                        scenario_path,
                        duality_handler = duality_handler,
                    )
                end
                push!(items.duals, subproblem_results.duals)
                push!(items.supports, noise)
                push!(items.nodes, child_node.index)
                push!(items.probability, child.probability * noise.probability)
                push!(items.objectives, subproblem_results.objective)
                push!(items.belief, belief)
                items.cached_solutions[(child.term, noise.term)] =
                    length(items.duals)
            end
        end
    end
    if length(scenario_path) == length_scenario_path
        # No-op. There weren't any children to solve.
    else
        # Drop the last element (i.e., the one we added).
        pop!(scenario_path)
    end
    return
end

"""
    SDDP.calculate_bound(
        model::PolicyGraph,
        state::Dict{Symbol,Float64},
        risk_measure = Expectation(),
    )

Calculate the lower bound (if minimizing, otherwise upper bound) of the problem
model at the point state, assuming the risk measure at the root node is
risk_measure.
"""
function calculate_bound(
    model::PolicyGraph{T},
    root_state::Dict{Symbol,Float64} = model.initial_root_state;
    risk_measure = Expectation(),
) where {T}
    # Initialization.
    noise_supports = Any[]
    probabilities = Float64[]
    objectives = Float64[]
    current_belief = initialize_belief(model)
    # Solve all problems that are children of the root node.
    for child in model.root_children
        if isapprox(child.probability, 0.0, atol = 1e-6)
            continue
        end
        node = model[child.term]
        for noise in node.noise_terms
            if node.objective_state !== nothing
                update_objective_state(
                    node.objective_state,
                    node.objective_state.initial_value,
                    noise.term,
                )
            end
            # Update belief state, etc.
            if node.belief_state !== nothing
                belief = node.belief_state::BeliefState{T}
                partition_index = belief.partition_index
                belief.updater(
                    belief.belief,
                    current_belief,
                    partition_index,
                    noise.term,
                )
            end
            subproblem_results = solve_subproblem(
                model,
                node,
                root_state,
                noise.term,
                Tuple{T,Any}[(child.term, noise.term)],
                duality_handler = nothing,
            )
            push!(objectives, subproblem_results.objective)
            push!(probabilities, child.probability * noise.probability)
            push!(noise_supports, noise.term)
        end
    end
    # Now compute the risk-adjusted probability measure:
    risk_adjusted_probability = similar(probabilities)
    offset = adjust_probability(
        risk_measure,
        risk_adjusted_probability,
        probabilities,
        noise_supports,
        objectives,
        model.objective_sense == MOI.MIN_SENSE,
    )
    # Finally, calculate the risk-adjusted value.
    return sum(
        obj * prob for (obj, prob) in zip(objectives, risk_adjusted_probability)
    ) + offset
end

struct IterationResult{T}
    pid::Int
    bound::Float64
    cumulative_value::Float64
    has_converged::Bool
    status::Symbol
    cuts::Dict{T,Vector{Any}}
end

function iteration(model::PolicyGraph{T}, options::Options) where {T}
    TimerOutputs.@timeit SDDP_TIMER "forward_pass" begin
        forward_trajectory = forward_pass(model, options, options.forward_pass)
        options.forward_pass_callback(forward_trajectory)
    end
    TimerOutputs.@timeit SDDP_TIMER "backward_pass" begin
        cuts = backward_pass(
            model,
            options,
            forward_trajectory.scenario_path,
            forward_trajectory.sampled_states,
            forward_trajectory.objective_states,
            forward_trajectory.belief_states,
        )
    end
    TimerOutputs.@timeit SDDP_TIMER "calculate_bound" begin
        bound = calculate_bound(model)
    end
    push!(
        options.log,
        Log(
            length(options.log) + 1,
            bound,
            forward_trajectory.cumulative_value,
            time() - options.start_time,
            Distributed.myid(),
            model.ext[:total_solves],
        ),
    )
    has_converged, status =
        convergence_test(model, options.log, options.stopping_rules)
    return IterationResult(
        Distributed.myid(),
        bound,
        forward_trajectory.cumulative_value,
        has_converged,
        status,
        cuts,
    )
end

"""
    termination_status(model::PolicyGraph)

Query the reason why the training stopped.
"""
function termination_status(model::PolicyGraph)
    if model.most_recent_training_results === nothing
        return :model_not_solved
    end
    return model.most_recent_training_results.status
end

"""
    SDDP.train(model::PolicyGraph; kwargs...)

Train the policy for `model`.

## Keyword arguments

 - `iteration_limit::Int`: number of iterations to conduct before termination.

 - `time_limit::Float64`: number of seconds to train before termination.

 - `stoping_rules`: a vector of [`SDDP.AbstractStoppingRule`](@ref)s.

 - `print_level::Int`: control the level of printing to the screen. Defaults to
    `1`. Set to `0` to disable all printing.

 - `log_file::String`: filepath at which to write a log of the training
   progress. Defaults to `SDDP.log`.

 - `log_frequency::Int`: control the frequency with which the logging is
    outputted (iterations/log). Defaults to `1`.

 - `run_numerical_stability_report::Bool`: generate (and print) a numerical
   stability report prior to solve. Defaults to `true`.

 - `refine_at_similar_nodes::Bool`: if SDDP can detect that two nodes have the
    same children, it can cheaply add a cut discovered at one to the other. In
    almost all cases this should be set to `true`.

 - `cut_deletion_minimum::Int`: the minimum number of cuts to cache before
    deleting  cuts from the subproblem. The impact on performance is solver
    specific; however, smaller values result in smaller subproblems (and
    therefore quicker solves), at the expense of more time spent performing cut
    selection.

 - `risk_measure`: the risk measure to use at each node. Defaults to
   [`Expectation`](@ref).

 - `sampling_scheme`: a sampling scheme to use on the forward pass of the
    algorithm. Defaults to [`InSampleMonteCarlo`](@ref).

 - `backward_sampling_scheme`: a backward pass sampling scheme to use on the
    backward pass of the algorithm. Defaults to `CompleteSampler`.

 - `cut_type`: choose between `SDDP.SINGLE_CUT` and `SDDP.MULTI_CUT` versions of
   SDDP.

 - `dashboard::Bool`: open a visualization of the training over time. Defaults
    to `false`.

 - `parallel_scheme::AbstractParallelScheme`: specify a scheme for solving in
   parallel. Defaults to `Serial()`.

 - `forward_pass::AbstractForwardPass`: specify a scheme to use for the forward
   passes.

 - `forward_pass_resampling_probability::Union{Nothing,Float64}`: set to a value
   in `(0, 1)` to enable [`RiskAdjustedForwardPass`](@ref). Defaults to
   `nothing` (disabled).

 - `add_to_existing_cuts::Bool`: set to `true` to allow training a model that
   was previously trained. Defaults to `false`.

 - `duality_handler::AbstractDualityHandler`: specify a duality handler to use
   when creating cuts.

There is also a special option for infinite horizon problems

 - `cycle_discretization_delta`: the maximum distance between states allowed on
    the forward pass. This is for advanced users only and needs to be used in
    conjunction with a different `sampling_scheme`.
"""
function train(
    model::PolicyGraph;
    iteration_limit::Union{Int,Nothing} = nothing,
    time_limit::Union{Real,Nothing} = nothing,
    print_level::Int = 1,
    log_file::String = "SDDP.log",
    log_frequency::Int = 1,
    run_numerical_stability_report::Bool = true,
    stopping_rules = AbstractStoppingRule[],
    risk_measure = SDDP.Expectation(),
    sampling_scheme = SDDP.InSampleMonteCarlo(),
    cut_type = SDDP.SINGLE_CUT,
    cycle_discretization_delta::Float64 = 0.0,
    refine_at_similar_nodes::Bool = true,
    cut_deletion_minimum::Int = 1,
    backward_sampling_scheme::AbstractBackwardSamplingScheme = SDDP.CompleteSampler(),
    dashboard::Bool = false,
    parallel_scheme::AbstractParallelScheme = Serial(),
    forward_pass::AbstractForwardPass = DefaultForwardPass(),
    forward_pass_resampling_probability::Union{Nothing,Float64} = nothing,
    add_to_existing_cuts::Bool = false,
    duality_handler::AbstractDualityHandler = SDDP.ContinuousConicDuality(),
    forward_pass_callback::Function = (x) -> nothing,
)
    if !add_to_existing_cuts && model.most_recent_training_results !== nothing
        @warn("""
        Re-training a model with existing cuts!

        Are you sure you want to do this? The output from this training may be
        misleading because the policy is already partially trained.

        If you meant to train a new policy with different settings, you must
        build a new model.

        If you meant to refine a previously trained policy, turn off this
        warning by passing `add_to_existing_cuts = true` as a keyword argument
        to `SDDP.train`.

        In a future release, this warning may turn into an error.
        """)
    end
    if forward_pass_resampling_probability !== nothing
        forward_pass = RiskAdjustedForwardPass(
            forward_pass = forward_pass,
            risk_measure = risk_measure,
            resampling_probability = forward_pass_resampling_probability,
        )
    end
    # Reset the TimerOutput.
    TimerOutputs.reset_timer!(SDDP_TIMER)
    log_file_handle = open(log_file, "a")
    log = Log[]

    if print_level > 0
        print_helper(print_banner, log_file_handle)
        print_helper(
            print_problem_statistics,
            log_file_handle,
            model,
            model.most_recent_training_results !== nothing,
            parallel_scheme,
            risk_measure,
            sampling_scheme,
        )
    end
    if run_numerical_stability_report
        report = sprint(
            io -> numerical_stability_report(
                io,
                model,
                print = print_level > 0,
            ),
        )
        print_helper(print, log_file_handle, report)
    end
    if print_level > 0
        print_helper(print_iteration_header, log_file_handle)
    end
    # Convert the vector to an AbstractStoppingRule. Otherwise if the user gives
    # something like stopping_rules = [SDDP.IterationLimit(100)], the vector
    # will be concretely typed and we can't add a TimeLimit.
    stopping_rules = convert(Vector{AbstractStoppingRule}, stopping_rules)
    # Add the limits as stopping rules. An IterationLimit or TimeLimit may
    # already exist in stopping_rules, but that doesn't matter.
    if iteration_limit !== nothing
        push!(stopping_rules, IterationLimit(iteration_limit))
    end
    if time_limit !== nothing
        push!(stopping_rules, TimeLimit(time_limit))
    end
    if length(stopping_rules) == 0
        @warn(
            "You haven't specified a stopping rule! You can only terminate " *
            "the call to SDDP.train via a keyboard interrupt ([CTRL+C])."
        )
    end
    # Update the nodes with the selected cut type (SINGLE_CUT or MULTI_CUT)
    # and the cut deletion minimum.
    if cut_deletion_minimum < 0
        cut_deletion_minimum = typemax(Int)
    end
    for (key, node) in model.nodes
        node.bellman_function.cut_type = cut_type
        node.bellman_function.global_theta.cut_oracle.deletion_minimum =
            cut_deletion_minimum
        for oracle in node.bellman_function.local_thetas
            oracle.cut_oracle.deletion_minimum = cut_deletion_minimum
        end
    end
    dashboard_callback = if dashboard
        launch_dashboard()
    else
        (::Any, ::Any) -> nothing
    end
    options = Options(
        model,
        model.initial_root_state,
        sampling_scheme,
        backward_sampling_scheme,
        risk_measure,
        cycle_discretization_delta,
        refine_at_similar_nodes,
        stopping_rules,
        dashboard_callback,
        print_level,
        time(),
        log,
        log_file_handle,
        log_frequency,
        forward_pass,
        duality_handler,
        forward_pass_callback,
    )
    status = :not_solved
    try
        status = master_loop(parallel_scheme, model, options)
    catch ex
        if isa(ex, InterruptException)
            status = :interrupted
            interrupt(parallel_scheme)
        else
            close(log_file_handle)
            rethrow(ex)
        end
    finally
        # And close the dashboard callback if necessary.
        dashboard_callback(nothing, true)
    end
    training_results = TrainingResults(status, log)
    model.most_recent_training_results = training_results
    if print_level > 0
        print_helper(print_footer, log_file_handle, training_results)
        if print_level > 1
            print_helper(TimerOutputs.print_timer, log_file_handle, SDDP_TIMER)
            # Annoyingly, TimerOutputs doesn't end the print section with `\n`,
            # so we do it here.
            print_helper(println, log_file_handle)
        end
    end
    close(log_file_handle)
    return
end

# Internal function: helper to conduct a single simulation. Users should use the
# documented, user-facing function SDDP.simulate instead.
function _simulate(
    model::PolicyGraph{T},
    variables::Vector{Symbol};
    sampling_scheme::AbstractSamplingScheme,
    custom_recorders::Dict{Symbol,Function},
    duality_handler::Union{Nothing,AbstractDualityHandler},
    skip_undefined_variables::Bool,
    incoming_state::Dict{Symbol,Float64},
) where {T}
    # Sample a scenario path.
    scenario_path, _ = sample_scenario(model, sampling_scheme)

    # Storage for the simulation results.
    simulation = Dict{Symbol,Any}[]
    current_belief = initialize_belief(model)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0

    # Objective state interpolation.
    objective_state_vector, N =
        initialize_objective_state(model[scenario_path[1][1]])
    objective_states = NTuple{N,Float64}[]
    for (depth, (node_index, noise)) in enumerate(scenario_path)
        node = model[node_index]
        # Objective state interpolation.
        objective_state_vector = update_objective_state(
            node.objective_state,
            objective_state_vector,
            noise,
        )
        if objective_state_vector !== nothing
            push!(objective_states, objective_state_vector)
        end
        if node.belief_state !== nothing
            belief = node.belief_state::BeliefState{T}
            partition_index = belief.partition_index
            current_belief = belief.updater(
                belief.belief,
                current_belief,
                partition_index,
                noise,
            )
        else
            current_belief = Dict(node_index => 1.0)
        end
        # Solve the subproblem.
        subproblem_results = solve_subproblem(
            model,
            node,
            incoming_state,
            noise,
            scenario_path[1:depth],
            duality_handler = duality_handler,
        )
        # Add the stage-objective
        cumulative_value += subproblem_results.stage_objective
        # Record useful variables from the solve.
        store = Dict{Symbol,Any}(
            :node_index => node_index,
            :noise_term => noise,
            :stage_objective => subproblem_results.stage_objective,
            :bellman_term =>
                subproblem_results.objective -
                subproblem_results.stage_objective,
            :objective_state => objective_state_vector,
            :belief => copy(current_belief),
        )
        if objective_state_vector !== nothing && N == 1
            store[:objective_state] = store[:objective_state][1]
        end
        # Loop through the primal variable values that the user wants.
        for variable in variables
            if haskey(node.subproblem.obj_dict, variable)
                # Note: we broadcast the call to value for variables which are
                # containers (like Array, Containers.DenseAxisArray, etc). If
                # the variable is a scalar (e.g. just a plain VariableRef), the
                # broadcast preseves the scalar shape.
                # TODO: what if the variable container is a dictionary? They
                # should be using Containers.SparseAxisArray, but this might not
                # always be the case...
                store[variable] = JuMP.value.(node.subproblem[variable])
            elseif skip_undefined_variables
                store[variable] = NaN
            else
                error(
                    "No variable named $(variable) exists in the subproblem.",
                    " If you want to simulate the value of a variable, make ",
                    "sure it is defined in _all_ subproblems, or pass ",
                    "`skip_undefined_variables=true` to `simulate`.",
                )
            end
        end
        # Loop through any custom recorders that the user provided.
        for (sym, recorder) in custom_recorders
            store[sym] = recorder(node.subproblem)
        end
        # Add the store to our list.
        push!(simulation, store)
        # Set outgoing state as the incoming state for the next node.
        incoming_state = copy(subproblem_results.state)
    end
    return simulation
end

function _initial_state(model::PolicyGraph)
    return Dict(String(k) => v for (k, v) in model.initial_root_state)
end

"""
    simulate(
        model::PolicyGraph,
        number_replications::Int = 1,
        variables::Vector{Symbol} = Symbol[];
        sampling_scheme::AbstractSamplingScheme =
            InSampleMonteCarlo(),
        custom_recorders = Dict{Symbol, Function}(),
        duality_handler::Union{Nothing,AbstractDualityHandler} = nothing,
        skip_undefined_variables::Bool = false,
        parallel_scheme::AbstractParallelScheme = Serial(),
        incoming_state::Dict{String,Float64} = _intial_state(model),
     )::Vector{Vector{Dict{Symbol,Any}}}

Perform a simulation of the policy model with `number_replications` replications
using the sampling scheme `sampling_scheme`.

Use `incoming_state` to pass an initial value of the state variable, if it
differs from that at the root node. Each key should be the string name of the
state variable.

Returns a vector with one element for each replication. Each element is a vector
with one-element for each node in the scenario that was sampled. Each element in
that vector is a dictionary containing information about the subproblem that was
solved.

In that dictionary there are four special keys:
 - :node_index, which records the index of the sampled node in the policy model
 - :noise_term, which records the noise observed at the node
 - :stage_objective, which records the stage-objective of the subproblem
 - :bellman_term, which records the cost/value-to-go of the node.
The sum of :stage_objective + :bellman_term will equal the objective value of
the solved subproblem.

In addition to the special keys, the dictionary will contain the result of
`JuMP.value(subproblem[key])` for each `key` in `variables`. This is
useful to obtain the primal value of the state and control variables.

For more complicated data, the `custom_recorders` keyword argument can be used.

```julia
data = Dict{Symbol, Any}()
for (key, recorder) in custom_recorders
    data[key] = foo(subproblem)
end
```

For example, to record the dual of a constraint named `my_constraint`, pass the
following:

```julia
simulation_results = SDDP.simulate(model, 2;
    custom_recorders = Dict{Symbol, Function}(
        :constraint_dual => (sp) -> JuMP.dual(sp[:my_constraint])
    )
)
```

The value of the dual in the first stage of the second replication can be
accessed as:

```julia
simulation_results[2][1][:constraint_dual]
```

If you do not require dual variables (or if they are not available), pass
`duality_handler = nothing`.

If you attempt to simulate the value of a variable that is only defined in some
of the stage problems, an error will be thrown. To over-ride this (and return a
`NaN` instead), pass `skip_undefined_variables = true`.

Use `parallel_scheme::[AbstractParallelScheme](@ref)` to specify a scheme for
simulating in parallel. Defaults to [`Serial`](@ref).
"""
function simulate(
    model::PolicyGraph,
    number_replications::Int = 1,
    variables::Vector{Symbol} = Symbol[];
    sampling_scheme::AbstractSamplingScheme = InSampleMonteCarlo(),
    custom_recorders = Dict{Symbol,Function}(),
    duality_handler::Union{Nothing,AbstractDualityHandler} = nothing,
    skip_undefined_variables::Bool = false,
    parallel_scheme::AbstractParallelScheme = Serial(),
    incoming_state::Dict{String,Float64} = _initial_state(model),
)
    return _simulate(
        model,
        parallel_scheme,
        number_replications,
        variables;
        sampling_scheme = sampling_scheme,
        custom_recorders = custom_recorders,
        duality_handler = duality_handler,
        skip_undefined_variables = skip_undefined_variables,
        incoming_state = Dict(Symbol(k) => v for (k, v) in incoming_state),
    )
end

"""
    DecisionRule(model::PolicyGraph{T}; node::T)

Create a decision rule for node `node` in `model`.
"""
struct DecisionRule{T}
    model::PolicyGraph{T}
    node::Node{T}
    function DecisionRule(model::PolicyGraph{T}; node::T) where {T}
        return new{T}(model, model[node])
    end
end

function Base.show(io::IO, pi::DecisionRule)
    print(io, "A decision rule for node $(pi.node.index)")
    return
end

"""
    evaluate(
        rule::DecisionRule;
        incoming_state::Dict{Symbol,Float64},
        noise = nothing,
        controls_to_record = Symbol[],
    )

Evalute the decision rule `rule` at the point described by the `incoming_state`
and `noise`.

If the node is deterministic, omit the `noise` argument.

Pass a list of symbols to `controls_to_record` to save the optimal primal
solution corresponding to the names registered in the model.
"""
function evaluate(
    rule::DecisionRule{T};
    incoming_state::Dict{Symbol,Float64},
    noise = nothing,
    controls_to_record = Symbol[],
) where {T}
    ret = solve_subproblem(
        rule.model,
        rule.node,
        incoming_state,
        noise,
        Tuple{T,Any}[];
        duality_handler = nothing,
    )
    return (
        stage_objective = ret.stage_objective,
        outgoing_state = ret.state,
        controls = Dict(
            c => value.(rule.node.subproblem[c]) for c in controls_to_record
        ),
    )
end
