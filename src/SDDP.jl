struct Options{SamplingScheme, T}
    initial_state::Dict{Symbol, Float64}
    sampling_scheme::SamplingScheme
    starting_states::Dict{T, Vector{Dict{Symbol, Float64}}}
    cycle_discretization_delta::Float64
    function Options(policy_graph::PolicyGraph{T},
                     initial_state::Dict{Symbol, Float64},
                     sampling_scheme::S,
                     cycle_discretization_delta::Float64) where {T, S}
        starting_states = Dict{T, Vector{Dict{Symbol, Float64}}}()
        for node_index in keys(policy_graph.nodes)
            starting_states[node_index] = Dict{Symbol, Float64}[]
        end
        return new{S, T}(initial_state,
                         sampling_scheme,
                         starting_states,
                         cycle_discretization_delta)
    end
end

# Internal function: set the incoming state variables of node to the values
# contained in state.
function set_incoming_state(node::Node, state::Dict{Symbol, Float64})
    for (state_name, value) in state
        JuMP.fix(node.states[state_name].incoming, value)
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
        outgoing_value = JuMP.result_value(state.outgoing)
        if JuMP.has_upper_bound(state.outgoing) && JuMP.upper_bound(state.outgoing) < outgoing_value
            outgoing_value = JuMP.upper_bound(state.outgoing)
        elseif JuMP.has_lower_bound(state.outgoing) && JuMP.lower_bound(state.outgoing) > outgoing_value
            outgoing_value = JuMP.lower_bound(state.outgoing)
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
        ref = JuMP.FixRef(state.incoming)
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
           stage_objective_value(node.stage_objective),  # C(x, u, ω)
           JuMP.objective_value(node.subproblem)  # C(x, u, ω) + θ
end

# Internal function: perform a single forward pass of the SDDP algorithm given
# options.
function forward_pass(graph::PolicyGraph{T}, options::Options) where T
    @debug "Beginning forward pass."
    scenario_path = sample_scenario(graph, options.sampling_scheme)
    @debug "Forward pass = $(scenario_path)."
    sampled_states = Dict{Symbol, Float64}[]
    state = copy(options.initial_state)
    cumulative_value = 0.0
    for (node_index, noise) in scenario_path
        @debug "Solving $(node_index) with noise=$(noise) on forward pass."
        node = graph[node_index]
        # ===== Begin: starting state for infinite horizon =====
        starting_states = options.starting_states[node_index]
        if distance(starting_states, state) > options.cycle_discretization_delta
            push!(starting_states, state)
        end
        num_starting_states = length(starting_states)
        state = splice!(starting_states, rand(1:num_starting_states))
        # ===== End: starting state for infinite horizon =====
        new_state, duals, stage_objective, objective = solve_subproblem(graph, node, state, noise)
        @debug "Solved $(node_index): $(new_state), $(duals), $(stage_objective), $(objective)"
        cumulative_value += stage_objective
        push!(sampled_states, new_state)
        state = copy(new_state)
    end
    # ===== Begin: drop off starting state if terminated due to cycle =====
    final_node_index = scenario_path[end][1]
    final_node = graph[final_node_index]
    if length(final_node.children) > 0  # Terminated due to cycle.
        if distance(
                options.starting_states[final_node_index],
                sampled_states[end-1]) > options.cycle_discretization_delta
            push!(options.starting_states[final_node_index], sampled_states[end-1])
        end
    end
    # ===== End: drop off starting state if terminated due to cycle =====
    return scenario_path, sampled_states, cumulative_value
end

function distance(starting_states, state)
    if length(starting_states) == 0
        return Inf
    else
        return maximum(inf_norm.(starting_states, Ref(state)))
    end
end

function inf_norm(x::Dict{Symbol, Float64}, y::Dict{Symbol, Float64})
    norm = 0.0
    for (key, value) in x
        if abs(value - y[key]) > norm
            norm = abs(value - y[key])
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
                new_outgoing_state, duals, stage_obj, obj = solve_subproblem(
                    graph, child_node, outgoing_state, noise.term)
                push!(dual_variables, duals)
                push!(noise_supports, noise.term)
                push!(original_probability, child.probability * noise.probability)
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
function calculate_bound(graph::PolicyGraph, state::Dict{Symbol, Float64},
                         risk_measure=Expectation())
    noise_supports = Any[]
    probabilities = Float64[]
    objectives = Float64[]
    for child in graph.root_children
        node = graph[child.term]
        for noise in node.noise_terms
            outgoing_state, duals, stage_obj, obj = solve_subproblem(
                graph, node, state, noise.term)
            push!(objectives, obj)
            push!(probabilities, child.probability * noise.probability)
            push!(noise_supports, noise.term)
        end
    end
    risk_adjusted_probability = similar(probabilities)
    adjust_probability(risk_measure,
                       risk_adjusted_probability,
                       probabilities,
                       noise_supports,
                       objectives,
                       graph.objective_sense == :Min)
    return sum(obj * prob for (obj, prob) in
        zip(objectives, risk_adjusted_probability))
end

function calculate_bound(graph::PolicyGraph, risk_measure::AbstractRiskMeasure=Expectation())
    calculate_bound(graph, graph.initial_root_state, risk_measure)
end

# Internal functions: helpers so users can pass arguments like:
# risk_measure = Kokako.Expectation(),
# risk_measure = Dict(1=>Expectation(), 2=>WorstCase())
# risk_measure = (node_index) -> node_index == 1 ? Expectation() : WorstCase()
get_per_node(collection, node_index) = collection
get_per_node(collection::Dict{T, <:Any}, node_index::T) where T = collection[node_index]
get_per_node(collection::Function, node_index::T) where T = collection(node_index)

function train(graph::PolicyGraph;
               iteration_limit = 100_000,
               time_limit = Inf,
               stopping_rules = AbstractStoppingRules[],
               # initial_state = nothing,
               # cut_selection = nothing,
               risk_measure = Kokako.Expectation(),
               sampling_scheme = Kokako.MonteCarlo(),
               print_level = 0,
               cycle_discretization_delta = 0.01
               # log_file = "log.txt",
               # cut_file = "cuts.csv"
               )
    if print_level > 0
        print_banner()
    end
    options = Options(graph, graph.initial_root_state, sampling_scheme, cycle_discretization_delta)
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
            scenario_path, sampled_states, cumulative_value = forward_pass(graph, options)
            backward_pass(graph, options, scenario_path, sampled_states, risk_measure)
            # bound = calculate_bound(graph, initial_state)
            bound = calculate_bound(graph)
            if print_level > 0
                print_iteration(iteration_count, cumulative_value, bound)
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
    return status
end

"""
Special keys:
 - :node_index
 - :noise_term
 - :stage_objective
 - :bellman_term
"""
function simulate(graph::PolicyGraph{T}, variables;
                  sampling_scheme = MonteCarlo(),
                  terminate_on_cycle::Bool = true,
                  max_depth::Int = 0) where T
    scenario_path = sample_scenario(graph, sampling_scheme;
        terminate_on_cycle = terminate_on_cycle, max_depth = max_depth)
    simulation = Dict{Symbol, Any}[]
    state = copy(graph.initial_root_state)
    cumulative_value = 0.0
    for (node_index, noise) in scenario_path
        node = graph[node_index]
        new_state, duals, stage_objective, objective = solve_subproblem(
            graph, node, state, noise)
        cumulative_value += stage_objective
        store = Dict{Symbol, Any}(
            :node_index => node_index,
            :noise_term => noise,
            :stage_objective => stage_objective,
            :bellman_term => objective - stage_objective
        )
        for variable in variables
            store[variable] = JuMP.result_value(node.subproblem[variable])
        end
        push!(simulation, store)
        state = copy(new_state)
    end
    return simulation
end

function simulate(graph, N::Int, variables::Vector{Symbol};
                  terminate_on_cycle::Bool = true, max_depth::Int = 0)
    return [
        simulate(graph, variables;
            sampling_scheme = MonteCarlo(),
            terminate_on_cycle = terminate_on_cycle,
            max_depth = max_depth)
        for i in 1:N]
end
