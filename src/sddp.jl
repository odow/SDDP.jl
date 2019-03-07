#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

const SDDP_TIMER = TimerOutputs.TimerOutput()

# to_nodal_form is an internal helper function so users can pass arguments like:
# risk_measure = Kokako.Expectation(),
# risk_measure = Dict(1=>Expectation(), 2=>WorstCase())
# risk_measure = (node_index) -> node_index == 1 ? Expectation() : WorstCase()
# It will return a dictionary with a key for each node_index in the policy
# graph, and a corresponding value of whatever the user provided.
function to_nodal_form(model::PolicyGraph{T}, element) where T
    # Note: we don't copy element here, so it element is mutable, you should use
    # to_nodal_form(model, x -> new_element()) instead. A good example is
    # Vector{T}; use to_nodal_form(model, i -> T[]).
    store = Dict{T, typeof(element)}()
    for node_index in keys(model.nodes)
        store[node_index] = element
    end
    return store
end
function to_nodal_form(model::PolicyGraph{T}, builder::Function) where T
    node = first(keys(model.nodes))
    element = builder(node)
    store = Dict{T, typeof(element)}()
    for node_index in keys(model.nodes)
        store[node_index] = builder(node_index)
    end
    return store
end
function to_nodal_form(model::PolicyGraph{T}, dict::Dict{T, V}) where {T, V}
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
#
# TODO(odow): this is inefficient as it is O(n²) in the number of nodes, but
# it's just a one-off hit so let's optimize later.
function get_same_children(model::PolicyGraph{T}) where T
    same_children = Dict{T, Vector{T}}()
    # For each node in the model
    for (node_index_1, node_1) in model.nodes
        same_children[node_index_1] = T[]
        # Get the set of child nodes.
        children_1 = Set(child.term for child in node_1.children)
        # Skip this one if there are no children.
        length(children_1) == 0 && continue
        # For each node in the model:
        for (node_index_2, node_2) in model.nodes
            node_index_1 == node_index_2 && continue
            # Get the set of child nodes.
            children_2 = Set(child.term for child in node_2.children)
            # Skip this one if there are no children.
            length(children_2) == 0 && continue
            # Record if node_1 has a superset of node_2's children.
            if children_2 ⊆ children_1
                push!(same_children[node_index_1], node_index_2)
            end
        end
    end
    return same_children
end

function build_Φ(model::PolicyGraph{T}) where T
    Φ = Dict{Tuple{T, T}, Float64}()
    for (node_index_1, node_1) in model.nodes
        for child in node_1.children
            Φ[(node_index_1, child.term)] = child.probability
        end
    end
    return Φ
end

# Internal struct: storage for SDDP options and cached data. Users shouldn't
# interact with this directly.
struct Options{T}
    # The initial state to start from the root node.
    initial_state::Dict{Symbol, Float64}
    # The sampling scheme to use on the forward pass.
    sampling_scheme::AbstractSamplingScheme
    # Storage for the set of possible sampling states at each node. We only use
    # this if there is a cycle in the policy graph.
    starting_states::Dict{T, Vector{Dict{Symbol, Float64}}}
    # Risk measure to use at each node.
    risk_measures::Dict{T, AbstractRiskMeasure}
    # The delta by which to check if a state is close to a previously sampled
    # state.
    cycle_discretization_delta::Float64
    # Flag to add cuts to similar nodes.
    refine_at_similar_nodes::Bool
    # The node transition matrix.
    Φ::Dict{Tuple{T, T}, Float64}
    # A list of nodes that contain a subset of the children of node i.
    similar_children::Dict{T, Vector{T}}
    # Internal function: users should never construct this themselves.
    function Options(policy_graph::PolicyGraph{T},
                     initial_state::Dict{Symbol, Float64},
                     sampling_scheme::AbstractSamplingScheme,
                     risk_measures,
                     cycle_discretization_delta::Float64,
                     refine_at_similar_nodes::Bool) where {T, S}
        return new{T}(
            initial_state,
            sampling_scheme,
            to_nodal_form(policy_graph, x -> Dict{Symbol, Float64}[]),
            to_nodal_form(policy_graph, risk_measures),
            cycle_discretization_delta,
            refine_at_similar_nodes,
            build_Φ(policy_graph),
            get_same_children(policy_graph)
        )
    end
end

# Internal function: set the incoming state variables of node to the values
# contained in state.
function set_incoming_state(node::Node, state::Dict{Symbol, Float64})
    for (state_name, value) in state
        JuMP.fix(node.states[state_name].in, value)
    end
    return
end

# Internal function: get the values of the outgoing state variables in node.
# Requires node.subproblem to have been solved with PrimalStatus ==
# FeasiblePoint.
function get_outgoing_state(node::Node)
    values = Dict{Symbol, Float64}()
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

# Internal function: get the values of the dual variables associated with the
# fixed incoming state variables. Requires node.subproblem to have been solved
# with DualStatus == FeasiblePoint.
function get_dual_variables(node::Node)
    # Note: due to JuMP's dual convention, we need to flip the sign for
    # maximization problems.
    dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1.0 : -1.0
    values = Dict{Symbol, Float64}()
    for (name, state) in node.states
        ref = JuMP.FixRef(state.in)
        values[name] = dual_sign * JuMP.dual(ref)
    end
    return values
end

# Internal function: set the objective of node to the stage objective, plus the
# cost/value-to-go term.
function set_objective(model::PolicyGraph{T}, node::Node{T}) where T
    objective_state_component = get_objective_state_component(node)
    if objective_state_component != JuMP.AffExpr(0.0)
        node.stage_objective_set = false
    end
    if !node.stage_objective_set
        JuMP.set_objective(
            node.subproblem,
            model.objective_sense,
            node.stage_objective + objective_state_component +
                bellman_term(node.bellman_function)
        )
    end
    node.stage_objective_set = true
end

# Internal function: overload for the case where JuMP.value fails on a
# Real number.
stage_objective_value(stage_objective::Real) = stage_objective
stage_objective_value(stage_objective) = JuMP.value(stage_objective)

function write_to_file(subproblem::JuMP.Model, filename::String)
    mps = MathOptFormat.MPS.Model()
    MOI.copy_to(mps, JuMP.backend(subproblem))
    MOI.write_to_file(mps, filename * ".mps")
    lp = MathOptFormat.LP.Model()
    MOI.copy_to(lp, JuMP.backend(subproblem))
    MOI.write_to_file(lp, filename * ".lp")
    return
end

# Internal function: solve the subproblem associated with node given the
# incoming state variables state and realization of the stagewise-independent
# noise term noise. If require_duals=true, also return the dual variables
# associated with the fixed constraint of the incoming state variables.
function solve_subproblem(model::PolicyGraph{T},
                          node::Node{T},
                          state::Dict{Symbol, Float64},
                          noise,
                          require_duals::Bool = true) where T
    # Parameterize the model. First, fix the value of the incoming state
    # variables. Then parameterize the model depending on `noise`. Finally,
    # set the objective. Note that we set the objective every time incase
    # the user calls set_stage_objective in the parameterize function.
    set_incoming_state(node, state)
    node.parameterize(noise)
    set_objective(model, node)
    JuMP.optimize!(node.subproblem)
    # Test for primal feasibility.
    primal_status = JuMP.primal_status(node.subproblem)
    if primal_status != JuMP.MOI.FEASIBLE_POINT
        write_to_file(node.subproblem, "subproblem")
        error("Unable to retrieve primal solution from ", node.index, ".",
              "\n  Termination status: ", JuMP.termination_status(node.subproblem),
              "\n  Primal status:      ", primal_status,
              "\n  Dual status:        ", JuMP.dual_status(node.subproblem),
              ".\n An MPS file was written to `subproblem.mps` and an LP file ",
              "written to `subproblem.lp`.")
    end
    # If require_duals = true, check for dual feasibility and return a dict with
    # the dual on the fixed constraint associated with each incoming state
    # variable. If require_duals=false, return an empty dictionary for
    # type-stability.
    dual_values = if require_duals
        dual_status = JuMP.dual_status(node.subproblem)
        if dual_status != JuMP.MOI.FEASIBLE_POINT
            write_to_file(node.subproblem, "subproblem.mps")
            error("Unable to retrieve dual solution from ", node.index, ".",
                  "\n  Termination status: ", JuMP.termination_status(node.subproblem),
                  "\n  Primal status:      ", primal_status,
                  "\n  Dual status:        ", dual_status,
                  ".\n An MPS file was written to `subproblem.mps`")
        end
        get_dual_variables(node)
    else
        Dict{Symbol, Float64}()
    end
    return get_outgoing_state(node),  # The outgoing state variable x'.
           dual_values,  # The dual variables on the incoming state variables.
           stage_objective_value(node.stage_objective),
           JuMP.objective_value(node.subproblem)  # C(x, u, ω) + θ
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
update_objective_state(obj_state::Nothing, current_state, noise) = nothing
function update_objective_state(obj_state, current_state, noise)
    if length(current_state) == 1
        obj_state.state = (obj_state.update(current_state[1], noise),)
    else
        obj_state.state = obj_state.update(current_state, noise)
    end
    return obj_state.state
end

# Internal function: perform a single forward pass of the SDDP algorithm given
# options.
function forward_pass(model::PolicyGraph{T}, options::Options) where T
    # First up, sample a scenario. Note that if a cycle is detected, this will
    # return the cycle node as well.
    TimerOutputs.@timeit SDDP_TIMER "sample_scenario" begin
        scenario_path, terminated_due_to_cycle = sample_scenario(
            model, options.sampling_scheme)
    end
    # Storage for the list of outgoing states that we visit on the forward pass.
    sampled_states = Dict{Symbol, Float64}[]
    # Our initial incoming state.
    incoming_state_value = copy(options.initial_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0
    # Objective state interpolation.
    objective_state_vector, N = initialize_objective_state(
        model[scenario_path[1][1]])
    objective_states = NTuple{N, Float64}[]
    # Iterate down the scenario.
    for (node_index, noise) in scenario_path
        node = model[node_index]
        # Objective state interpolation.
        objective_state_vector = update_objective_state(
            node.objective_state, objective_state_vector, noise)
        if objective_state_vector !== nothing
            push!(objective_states, objective_state_vector)
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
            incoming_state_value = splice!(
                starting_states, rand(1:length(starting_states))
            )
        end
        # ===== End: starting state for infinite horizon =====
        # Solve the subproblem, note that `require_duals = false`.
        TimerOutputs.@timeit SDDP_TIMER "solve_subproblem" begin
            (outgoing_state_value, duals, stage_objective, objective) =
                solve_subproblem(
                    model, node, incoming_state_value, noise, false)
        end
        # Cumulate the stage_objective.
        cumulative_value += stage_objective
        # Add the outgoing state variable to the list of states we have sampled
        # on this forward pass.
        push!(sampled_states, outgoing_state_value)
        # Set the outgoing state value as the incoming state value for the next
        # node.
        incoming_state_value = copy(outgoing_state_value)
    end
    if terminated_due_to_cycle
        # Get the last node in the scenario.
        final_node_index = scenario_path[end][1]
        # We terminated due to a cycle. Here is the list of possible starting
        # states for that node:
        starting_states = options.starting_states[final_node_index]
        # We also need the incoming state variable to the final node, which is
        # the outgoing state value of the last node:
        incoming_state_value = sampled_states[end]
        # If this incoming state value is more than δ away from another state,
        # add it to the list.
        if distance(starting_states, incoming_state_value) >
                options.cycle_discretization_delta
            push!(starting_states, incoming_state_value)
        end
    end
    # ===== End: drop off starting state if terminated due to cycle =====
    return scenario_path, sampled_states, objective_states, cumulative_value
end

# Internal function: calculate the minimum distance between the state `state`
# and the list of states in `starting_states` using the distance measure `norm`.
function distance(starting_states::Vector{Dict{Symbol, Float64}},
                  state::Dict{Symbol, Float64},
                  norm::Function = inf_norm)
    if length(starting_states) == 0
        return Inf
    else
        return minimum(norm.(starting_states, Ref(state)))
    end
end

# Internal function: the norm to use when checking the distance between two
# possible starting states. We're going to use: d(x, y) = |x - y| / (1 + |y|).
function inf_norm(x::Dict{Symbol, Float64}, y::Dict{Symbol, Float64})
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
function backward_pass(model::PolicyGraph{T},
                       options::Options,
                       scenario_path::Vector{Tuple{T, NoiseType}},
                       sampled_states::Vector{Dict{Symbol, Float64}},
                       objective_states::Vector{NTuple{N, Float64}}
                           ) where {T, NoiseType, N}
    for index in length(scenario_path):-1:1
        # Lookup node, noise realization, and outgoing state variables.
        node_index, noise = scenario_path[index]
        outgoing_state = sampled_states[index]
        node = model[node_index]
        # If our node has no children, it means that we terminated the forward
        # pass at a leaf node. In this case, we don't need to add any cuts so we
        # can skip back up the scenario path one node. This should only ever be
        # true on the last node, but it probably doesn't hurt to check every
        # time in case someone wants to implement a really weird sampling
        # scheme.
        if length(node.children) == 0
            continue
        end
        # Initialization.
        noise_supports = Noise[]
        child_indices = T[]
        original_probability = Float64[]
        dual_variables = Dict{Symbol, Float64}[]
        objective_realizations = Float64[]
        if length(node.children) == 0
            error("The `scenario_path` passed to the backward pass should not" *
                  " contain a leaf node.")
        end
        # Solve all children.
        for child in node.children
            child_node = model[child.term]
            for noise in child_node.noise_terms
                # There are multiple ways to check this. Here is one. Another
                # could be that |objective_state| = |scenario_path|.
                if child_node.objective_state !== nothing
                    update_objective_state(child_node.objective_state,
                        objective_states[index], noise.term)
                end
                TimerOutputs.@timeit SDDP_TIMER "solve_subproblem" begin
                    (new_outgoing_state, duals, stage_objective, obj) =
                        solve_subproblem(
                            model, child_node, outgoing_state, noise.term
                        )
                end
                push!(dual_variables, duals)
                push!(noise_supports, noise)
                push!(child_indices, child_node.index)
                push!(original_probability,
                    child.probability * noise.probability)
                push!(objective_realizations, obj)
            end
        end
        refine_bellman_function(
            model,
            node,
            node.bellman_function,
            options.risk_measures[node_index],
            outgoing_state,
            dual_variables,
            noise_supports,
            original_probability,
            objective_realizations
        )
        if options.refine_at_similar_nodes
            # Refine the bellman function at other nodes with the same children,
            # e.g., in the same stage of a Markovian policy graph.
            for other_index in options.similar_children[node_index]
                copied_probability = similar(original_probability)
                other_node = model[other_index]
                for (idx, child_index) in enumerate(child_indices)
                    copied_probability[idx] =
                        get(options.Φ, (other_index, child_index), 0.0) *
                        noise_supports[idx].probability
                end
                refine_bellman_function(
                    model,
                    other_node,
                    other_node.bellman_function,
                    options.risk_measures[other_index],
                    outgoing_state,
                    dual_variables,
                    noise_supports,
                    copied_probability,
                    objective_realizations
                )
            end
        end
    end
end

"""
    Kokako.calculate_bound(model::PolicyGraph, state::Dict{Symbol, Float64},
                           risk_measure=Expectation())

Calculate the lower bound (if minimizing, otherwise upper bound) of the problem
model at the point state, assuming the risk measure at the root node is
risk_measure.
"""
function calculate_bound(model::PolicyGraph,
                         root_state::Dict{Symbol, Float64} =
                            model.initial_root_state;
                         risk_measure = Expectation())
    # Initialization.
    noise_supports = Any[]
    probabilities = Float64[]
    objectives = Float64[]
    # Solve all problems that are children of the root node.
    for child in model.root_children
        node = model[child.term]
        for noise in node.noise_terms
            if node.objective_state !== nothing
                update_objective_state(node.objective_state,
                    node.objective_state.initial_value, noise.term)
            end
            (outgoing_state, duals, stage_objective, obj) =
                solve_subproblem(model, node, root_state, noise.term)
            push!(objectives, obj)
            push!(probabilities, child.probability * noise.probability)
            push!(noise_supports, noise.term)
        end
    end
    # Now compute the risk-adjusted probability measure:
    risk_adjusted_probability = similar(probabilities)
    adjust_probability(risk_measure,
                       risk_adjusted_probability,
                       probabilities,
                       noise_supports,
                       objectives,
                       model.objective_sense == MOI.MIN_SENSE)
    # Finally, calculate the risk-adjusted value.
    return sum(obj * prob for (obj, prob) in
        zip(objectives, risk_adjusted_probability))
end

struct TrainingResults
    status::Symbol
    log::Vector{Log}
end

"""
    termination_status(model::PolicyGraph)

Query the reason why the training stopped.
"""
function termination_status(model::PolicyGraph)
    if model.most_recent_training_results === nothing
        return :ModelNotSolved
    else
        return model.most_recent_training_results.status
    end
end

"""
    Kokako.train(model::PolicyGraph; kwargs...)

Train the policy of the model. Keyword arguments are
 - iteration_limit: number of iterations to conduct before termination
 - time_limit: number of seconds to train before termination
 - print_level: control the level of printing to the screen
 - log_file: filepath at which to write a log of the training progress
 - perform_numerical_stability_check: perform a numerical stability check prior
   to solve
 - risk_measure
 - stoping_rules
 - sampling_scheme: a sampling scheme to use on the forward pass of the
   algorithm. Defaults to InSampleMonteCarlo().
 - refine_at_similar_nodes

There is also a special option for infinite horizon problems
 - cycle_discretization_delta: the maximum distance between states allowed on
   the forward pass.
"""
function train(model::PolicyGraph;
               iteration_limit = nothing,
               time_limit = nothing,
               print_level = 1,
               log_file = "kokako.log",
               perform_numerical_stability_check::Bool = true,
               stopping_rules = AbstractStoppingRule[],
               risk_measure = Kokako.Expectation(),
               sampling_scheme = Kokako.InSampleMonteCarlo(),
               cycle_discretization_delta = 0.0,
               refine_at_similar_nodes = true
               )
    # Reset the TimerOutput.
    TimerOutputs.reset_timer!(SDDP_TIMER)
    log_file_handle = open(log_file, "a")

    if print_level > 0
        print_banner()
        print_banner(log_file_handle)
    end

    if perform_numerical_stability_check
        # TODO: don't do this check twice.
        numerical_stability_report(model, print = print_level > 0)
        numerical_stability_report(
            log_file_handle, model, print = print_level > 0)
    end

    if print_level > 0
        print_iteration_header()
        print_iteration_header(log_file_handle)
    end
    # Convert the vector to an AbstractStoppingRule. Otherwise if the user gives
    # something like stopping_rules = [Kokako.IterationLimit(100)], the vector
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
        @warn("You haven't specified a stopping rule! You can only terminate " *
              "the call to Kokako.train via a keyboard interrupt ([CTRL+C]).")
    end
    options = Options(
        model,
        model.initial_root_state,
        sampling_scheme,
        risk_measure,
        cycle_discretization_delta,
        refine_at_similar_nodes
    )
    # The default status. This should never be seen by the user.
    status = :not_solved
    log = Log[]
    try
        start_time = time()
        iteration_count = 1
        has_converged = false
        while !has_converged
            TimerOutputs.@timeit SDDP_TIMER "forward_pass" begin
                scenario_path, sampled_states, objective_states,
                    cumulative_value = forward_pass(model, options)
            end
            TimerOutputs.@timeit SDDP_TIMER "backward_pass" begin
                backward_pass(model, options, scenario_path, sampled_states,
                              objective_states)
            end
            TimerOutputs.@timeit SDDP_TIMER "calculate_bound" begin
                bound = calculate_bound(model)
            end
            push!(log, Log(iteration_count, bound, cumulative_value,
                time() - start_time)
            )
            has_converged, status = convergence_test(model, log, stopping_rules)
            if print_level > 0
                print_iteration(stdout, log[end])
                print_iteration(log_file_handle, log[end])
            end
            iteration_count += 1
        end
    catch ex
        if isa(ex, InterruptException)
            status = :interrupted
        else
            close(log_file_handle)
            rethrow(ex)
        end
    end
    training_results = TrainingResults(status, log)
    model.most_recent_training_results = training_results
    if print_level > 0
        print_footer(stdout, training_results)
        print_footer(log_file_handle, training_results)
        if print_level > 1
            TimerOutputs.print_timer(stdout, SDDP_TIMER)
            TimerOutputs.print_timer(log_file_handle, SDDP_TIMER)
        end
    end
    close(log_file_handle)
    return
end

# Internal function: helper to conduct a single simulation. Users should use the
# documented, user-facing function Kokako.simulate instead.
function _simulate(model::PolicyGraph,
                   variables::Vector{Symbol} = Symbol[];
                   sampling_scheme::AbstractSamplingScheme =
                       InSampleMonteCarlo(),
                   custom_recorders = Dict{Symbol, Function}())
    # Sample a scenario path.
    scenario_path, terminated_due_to_cycle = sample_scenario(
        model, sampling_scheme
    )
    # Storage for the simulation results.
    simulation = Dict{Symbol, Any}[]
    # The incoming state values.
    incoming_state = copy(model.initial_root_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0

    # Objective state interpolation.
    objective_state_vector, N = initialize_objective_state(
        model[scenario_path[1][1]])
    objective_states = NTuple{N, Float64}[]
    for (node_index, noise) in scenario_path
        node = model[node_index]
        # Objective state interpolation.
        objective_state_vector = update_objective_state(node.objective_state,
            objective_state_vector, noise)
        if objective_state_vector !== nothing
            push!(objective_states, objective_state_vector)
        end
        # Solve the subproblem.
        outgoing_state, duals, stage_objective, objective = solve_subproblem(
            model, node, incoming_state, noise)
        # Add the stage-objective
        cumulative_value += stage_objective
        # Record useful variables from the solve.
        store = Dict{Symbol, Any}(
            :node_index => node_index,
            :noise_term => noise,
            :stage_objective => stage_objective,
            :bellman_term => objective - stage_objective,
            :objective_state => objective_state_vector
        )
        if objective_state_vector !== nothing && N == 1
            store[:objective_state] = store[:objective_state][1]
        end
        # Loop through the primal variable values that the user wants.
        for variable in variables
            if haskey(node.subproblem.obj_dict, variable)
                store[variable] = JuMP.value(node.subproblem[variable])
            else
                error("No variable named $(variable) exists in the subproblem" *
                     ". If you want to simulate the value of a variable, make" *
                     " sure it is defined in _all_ subproblems.")
            end
        end
        # Loop through any custom recorders that the user provided.
        for (sym, foo) in custom_recorders
            store[sym] = foo(node.subproblem)
        end
        # Add the store to our list.
        push!(simulation, store)
        # Set outgoing state as the incoming state for the next node.
        incoming_state = copy(outgoing_state)
    end
    return simulation
end

"""
    simulate(model::PolicyGraph,
             number_replications::Int = 1,
             variables::Vector{Symbol} = Symbol[];
             sampling_scheme::AbstractSamplingScheme =
                 InSampleMonteCarlo(),
             custom_recorders = Dict{Symbol, Function}()
     )::Vector{Vector{Dict{Symbol, Any}}}

Perform a simulation of the policy model with `number_replications` replications
using the sampling scheme `sampling_scheme`.

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

For more complicated data, the `custom_recorders` keyword arguement can be used.

    data = Dict{Symbol, Any}()
    for (key, recorder) in custom_recorders
        data[key] = foo(subproblem)
    end

For example, to record the dual of a constraint named `my_constraint`, pass the
following:

    simulation_results = simulate(model, number_replications=2;
        custom_recorders = Dict(
            :constraint_dual = (sp) -> JuMP.dual(sp[:my_constraint])
        )
    )

The value of the dual in the first stage of the second replication can be
accessed as:

    simulation_results[2][1][:constraint_dual]
"""
function simulate(model::PolicyGraph,
                  number_replications::Int = 1,
                  variables::Vector{Symbol} = Symbol[];
                  sampling_scheme::AbstractSamplingScheme =
                      InSampleMonteCarlo(),
                  custom_recorders = Dict{Symbol, Function}())
    return [_simulate(
                model,
                variables;
                sampling_scheme = sampling_scheme,
                custom_recorders = custom_recorders)
                for i in 1:number_replications
            ]
end
