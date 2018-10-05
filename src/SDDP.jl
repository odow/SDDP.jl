const SDDP_TIMER = TimerOutputs.TimerOutput()

struct Options{SamplingScheme, T}
    # The initial state to start from the root node.
    initial_state::Dict{Symbol, Float64}
    # The sampling scheme to use on the forward pass.
    sampling_scheme::SamplingScheme
    # Storage for the set of possible sampling states at each node. We only use
    # this if there is a cycle in the policy graph.
    starting_states::Dict{T, Vector{Dict{Symbol, Float64}}}
    # The delta by which to check if a state is close to a previously sampled
    # state.
    cycle_discretization_delta::Float64
    # Internal function: users should never construct this themselves.
    function Options(policy_graph::PolicyGraph{T},
                     initial_state::Dict{Symbol, Float64},
                     sampling_scheme::S,
                     cycle_discretization_delta::Float64 = 0.0) where {T, S}
        # Initialize the storage for every node in the policy graph.
        starting_states = Dict{T, Vector{Dict{Symbol, Float64}}}()
        for node_index in keys(policy_graph.nodes)
            starting_states[node_index] = Dict{Symbol, Float64}[]
        end
        return new{S, T}(initial_state, sampling_scheme, starting_states,
            cycle_discretization_delta)
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
        # is outside its bounds, project the value back onto the bounds.
        outgoing_value = JuMP.result_value(state.out)
        if (JuMP.has_upper_bound(state.out) &&
                JuMP.upper_bound(state.out) < outgoing_value)
            outgoing_value = JuMP.upper_bound(state.out)
        elseif (JuMP.has_lower_bound(state.out) &&
                JuMP.lower_bound(state.out) > outgoing_value)
            outgoing_value = JuMP.lower_bound(state.out)
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
    dual_sign = JuMP.objective_sense(node.subproblem) == :Min ? 1.0 : -1.0
    values = Dict{Symbol, Float64}()
    for (name, state) in node.states
        ref = JuMP.FixRef(state.in)
        values[name] = dual_sign * JuMP.result_dual(ref)
    end
    return values
end


# Internal function: set the objective of node to the stage objective, plus the
# cost/value-to-go term.
function set_objective(graph::PolicyGraph{T}, node::Node{T}) where T
    JuMP.set_objective(
        node.subproblem,
        graph.objective_sense,
        node.stage_objective + bellman_term(node.bellman_function)
    )
end

# Internal function: overload for the case where JuMP.result_value fails on a
# Real number.
stage_objective_value(stage_objective::Real) = stage_objective
stage_objective_value(stage_objective) = JuMP.result_value(stage_objective)

# Internal function: solve the subproblem associated with node given the
# incoming state variables state and realization of the stagewise-independent
# noise term noise. If require_duals=true, also return the dual variables
# associated with the fixed constraint of the incoming state variables.
function solve_subproblem(graph::PolicyGraph{T},
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
    # TODO(odow): cache the need to call set_objective. Only call it if the
    # stage-objective changes.
    set_objective(graph, node)
    JuMP.optimize!(node.subproblem)
    # Test for primal feasibility.
    primal_status = JuMP.primal_status(node.subproblem)
    if primal_status != JuMP.MOI.FeasiblePoint
        error("Unable to solve node $(node.index). Primal status: " *
              "$(primal_status).")
    end
    # If require_duals = true, check for dual feasibility and return a dict with
    # the dual on the fixed constraint associated with each incoming state
    # variable. If require_duals=false, return an empty dictionary for
    # type-stability.
    dual_values = if require_duals
        dual_status = JuMP.dual_status(node.subproblem)
        if dual_status != JuMP.MOI.FeasiblePoint
            error("Unable to solve dual of node $(node.index). Dual status: " *
                  "$(dual_status).")
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

# Internal function: perform a single forward pass of the SDDP algorithm given
# options.
function forward_pass(graph::PolicyGraph{T}, options::Options) where T
    # First up, sample a scenario. Note that if a cycle is detected, this will
    # return the cycle node as well.
    TimerOutputs.@timeit SDDP_TIMER "sample_scenario" begin
        scenario_path = sample_scenario(graph, options.sampling_scheme,
            terminate_on_cycle = true)
    end
    # Storage for the list of outgoing states that we visit on the forward pass.
    sampled_states = Dict{Symbol, Float64}[]
    # Our initial incoming state.
    incoming_state_value = copy(options.initial_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0
    # Iterate down the scenario.
    for (node_index, noise) in scenario_path
        node = graph[node_index]
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
        # Solve the subproblem, note that `require_duals=false`.
        TimerOutputs.@timeit SDDP_TIMER "solve_subproblem" begin
            (outgoing_state_value, duals, stage_objective, objective) =
                solve_subproblem(
                    graph, node, incoming_state_value, noise, false)
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
    # Get the last node in the scenario.
    final_node_index = scenario_path[end][1]
    final_node = graph[final_node_index]
    if length(final_node.children) > 0
        # The last node in the scenario has children, so we must have terminated
        # due to a cycle. Here is the list of possible starting states for that
        # node:
        starting_states = options.starting_states[final_node_index]
        # We also need the incoming state variable to the final node, which is
        # the outgoing state value of the 2'nd to last node:
        incoming_state_value = sampled_states[end-1]
        # If this incoming state value is more than δ away from another state,
        # add it to the list.
        if distance(starting_states, incoming_state_value) >
                options.cycle_discretization_delta
            push!(starting_states, incoming_state_value)
        end
    end
    # ===== End: drop off starting state if terminated due to cycle =====
    return scenario_path, sampled_states, cumulative_value
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
# scenario_path, refining the bellman function at sampled_states.
function backward_pass(graph::PolicyGraph{T},
                       options::Options,
                       scenario_path::Vector{Tuple{T, NoiseType}},
                       sampled_states::Vector{Dict{Symbol, Float64}},
                       risk_measures
                           ) where {T, NoiseType}
    for index in (length(scenario_path)-1):-1:1
        # Lookup node, noise realization, and outgoing state variables.
        node_index, noise = scenario_path[index]
        outgoing_state = sampled_states[index]
        node = graph[node_index]
        # Initialization.
        noise_supports = NoiseType[]
        original_probability = Float64[]
        dual_variables = Dict{Symbol, Float64}[]
        objective_realizations = Float64[]
        # Solve all children.
        for child in node.children
            child_node = graph[child.term]
            for noise in child_node.noise_terms
                TimerOutputs.@timeit SDDP_TIMER "solve_subproblem" begin
                    (new_outgoing_state, duals, stage_objective, obj) =
                        solve_subproblem(
                            graph, child_node, outgoing_state, noise.term
                        )
                end
                push!(dual_variables, duals)
                push!(noise_supports, noise.term)
                push!(original_probability,
                    child.probability * noise.probability)
                push!(objective_realizations, obj)
            end
        end
        # TODO(odow): refine the bellman function at other nodes with the same
        # children, e.g., in the same stage of a Markovian policy graph.
        refine_bellman_function(
            graph,
            node,
            node.bellman_function,
            get_per_node(risk_measures, node_index),
            outgoing_state,
            dual_variables,
            noise_supports,
            original_probability,
            objective_realizations)
    end
end

"""
    Kokako.calculate_bound(graph::PolicyGraph, state::Dict{Symbol, Float64},
                           risk_measure=Expectation())

Calculate the lower bound (if minimizing, otherwise upper bound) of the problem
graph at the point state, assuming the risk measure at the root node is
risk_measure.
"""
function calculate_bound(graph::PolicyGraph,
                         root_state::Dict{Symbol, Float64} =
                            graph.initial_root_state;
                         risk_measure = Expectation())
    # Initialization.
    noise_supports = Any[]
    probabilities = Float64[]
    objectives = Float64[]
    # Solve all problems that are children of the root node.
    for child in graph.root_children
        node = graph[child.term]
        for noise in node.noise_terms
            (outgoing_state, duals, stage_objective, obj) =
                solve_subproblem(graph, node, root_state, noise.term)
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
                       graph.objective_sense == :Min)
    # Finally, calculate the risk-adjusted value.
    return sum(obj * prob for (obj, prob) in
        zip(objectives, risk_adjusted_probability))
end

# Internal functions: helpers so users can pass arguments like:
# risk_measure = Kokako.Expectation(),
# risk_measure = Dict(1=>Expectation(), 2=>WorstCase())
# risk_measure = (node_index) -> node_index == 1 ? Expectation() : WorstCase()
get_per_node(collection, node_index) = collection
function get_per_node(collection::Dict{T, <:Any}, node_index::T) where T
    return collection[node_index]
end
function get_per_node(collection::Function, node_index::T) where T
    return collection(node_index)
end

"""
    Kokako.train(graph::PolicyGraph; kwargs...)

Train the policy of the graph. Keyword arguments are
 - iteration_limit: number of iterations to conduct before termination. Defaults
   to 100_000.
 - time_limit: number of seconds to train before termination. Defaults to Inf.
 - print_level: control the level of printing to the screen.
 - sampling_scheme: a sampling scheme to use on the forward pass of the
   algorithm. Defaults to InSampleMonteCarlo().

There is also a special option for infinite horizon problems
 - cycle_discretization_delta: the maximum distance between states allowed on
   the forward pass.
"""
function train(graph::PolicyGraph;
               iteration_limit = 100_000,
               time_limit = Inf,
               stopping_rules = AbstractStoppingRule[],
               risk_measure = Kokako.Expectation(),
               sampling_scheme = Kokako.InSampleMonteCarlo(),
               print_level = 0,
               cycle_discretization_delta = 0.01
               )
    # Reset the TimerOutput.
    TimerOutputs.reset_timer!(SDDP_TIMER)
    if print_level > 0
        print_banner()
    end
    options = Options(graph, graph.initial_root_state, sampling_scheme,
        cycle_discretization_delta)
    status = :not_solved
    start_time = time()
    iteration_count = 1
    try
        while true
            @debug "Beginning iteration $(iteration_count)."
            should_stop, status = convergence_test(graph, stopping_rules)
            if should_stop
                break
            end
            if time() - start_time > time_limit
                status = :time_limit
                break
            end
            if iteration_count > iteration_limit
                status = :iteration_limit
                break
            end
            TimerOutputs.@timeit SDDP_TIMER "forward_pass" begin
                scenario_path, sampled_states, cumulative_value = forward_pass(
                    graph, options)
            end
            TimerOutputs.@timeit SDDP_TIMER "backward_pass" begin
                backward_pass(
                    graph, options, scenario_path, sampled_states, risk_measure)
            end
            TimerOutputs.@timeit SDDP_TIMER "calculate_bound" begin
                bound = calculate_bound(graph)
            end
            if print_level > 0
                print_iteration(iteration_count, cumulative_value, bound,
                    time() - start_time)
            end
            iteration_count += 1
        end
    catch ex
        if isa(ex, InterruptException)
            @warn("Terminating solve due to user interaction.")
            status = :interrupted
        else
            rethrow(ex)
        end
    end
    if print_level > 1
        TimerOutputs.print_timer(stdout, SDDP_TIMER)
    end
    return status
end

# Internal function: helper to conduct a single simulation. Users should use the
# documented, user-facing function Kokako.simulate instead.
function _simulate(graph::PolicyGraph,
                   variables::Vector{Symbol} = Symbol[];
                   sampling_scheme::AbstractSamplingScheme =
                       InSampleMonteCarlo(),
                   terminate_on_cycle::Bool = true,
                   max_depth::Int = 0,
                   custom_recorders = Dict{Symbol, Function}())
    # Sample a scenario path.
    scenario_path = sample_scenario(graph, sampling_scheme;
        terminate_on_cycle = terminate_on_cycle,
        max_depth = max_depth
    )
    # Storage for the simulation results.
    simulation = Dict{Symbol, Any}[]
    # The incoming state values.
    incoming_state = copy(graph.initial_root_state)
    # A cumulator for the stage-objectives.
    cumulative_value = 0.0
    for (node_index, noise) in scenario_path
        node = graph[node_index]
        # Solve the subproblem.
        outgoing_state, duals, stage_objective, objective = solve_subproblem(
            graph, node, incoming_state, noise)
        # Add the stage-objective
        cumulative_value += stage_objective
        # Record useful variables from the solve.
        store = Dict{Symbol, Any}(
            :node_index => node_index,
            :noise_term => noise,
            :stage_objective => stage_objective,
            :bellman_term => objective - stage_objective
        )
        # Loop through the primal variable values that the user wants.
        for variable in variables
            store[variable] = JuMP.result_value(node.subproblem[variable])
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
    simulate(graph::PolicyGraph,
             number_replications::Int = 1,
             variables::Vector{Symbol} = Symbol[];
             sampling_scheme::AbstractSamplingScheme =
                 InSampleMonteCarlo(),
             terminate_on_cycle::Bool = true,
             max_depth::Int = 0,
             custom_recorders = Dict{Symbol, Function}()
     )::Vector{Vector{Dict{Symbol, Any}}}

Perform a simulation of the policy graph with `number_replications` replications
using the sampling scheme `sampling_scheme`. If `terminate_on_cycle=false` then
`max_depth` must be set >0.

Returns a vector with one element for each replication. Each element is a vector
with one-element for each node in the scenario that was sampled. Each element in
that vector is a dictionary containing information about the subproblem that was
solved.

In that dictionary there are four special keys:
 - :node_index, which records the index of the sampled node in the policy graph
 - :noise_term, which records the noise observed at the node
 - :stage_objective, which records the stage-objective of the subproblem
 - :bellman_term, which records the cost/value-to-go of the node.
The sum of :stage_objective + :bellman_term will equal the objective value of
the solved subproblem.

In addition to the special keys, the dictionary will contain the result of
`JuMP.result_value(subproblem[key])` for each `key` in `variables`. This is
useful to obtain the primal value of the state and control variables.

For more complicated data, the `custom_recorders` keyword arguement can be used.

    data = Dict{Symbol, Any}()
    for (key, recorder) in custom_recorders
        data[key] = foo(subproblem)
    end

For example, to record the dual of a constraint named `my_constraint`, pass the
following:

    simulation_results = simulate(graph, number_replications=2;
        custom_recorders = Dict(
            :constraint_dual = (sp) -> JuMP.result_dual(sp[:my_constraint])
        )
    )

The value of the dual in the first stage of the second replication can be
accessed as:

    simulation_results[2][1][:constraint_dual]
"""
function simulate(graph::PolicyGraph,
                  number_replications::Int = 1,
                  variables::Vector{Symbol} = Symbol[];
                  sampling_scheme::AbstractSamplingScheme =
                      InSampleMonteCarlo(),
                  terminate_on_cycle::Bool = true,
                  max_depth::Int = 0,
                  custom_recorders = Dict{Symbol, Function}())
    return [_simulate(
                graph,
                variables;
                sampling_scheme = sampling_scheme,
                terminate_on_cycle = terminate_on_cycle,
                max_depth = max_depth,
                custom_recorders = custom_recorders)
                for i in 1:number_replications
            ]
end
