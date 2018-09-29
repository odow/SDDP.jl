struct Options{SamplingScheme}
    initial_state::Dict{Symbol, Float64}
    sampling_scheme::SamplingScheme
end

function set_incoming_state(node::Node, state::Dict{Symbol, Float64})
    for (state_name, value) in state
        JuMP.fix(node.states[state_name].incoming, value)
    end
    return
end

function get_outgoing_state(node::Node)
    values = Dict{Symbol, Float64}()
    for (name, state) in node.states
        values[name] = JuMP.result_value(state.outgoing)
    end
    return values
end

function get_dual_variables(node::Node)
    values = Dict{Symbol, Float64}()
    for (name, state) in node.states
        values[name] = JuMP.result_dual(state.incoming)
    end
    return values
end

function set_objective(node::Node)
    JuMP.set_objective(node.subproblem, node.optimization_sense,
                      node.stage_objective)  # TODO(odow) + bellman_term(node.bellman_function))
end

function solve_subproblem(graph::PolicyGraph{T}, node, state, noise,
                          require_duals::Bool = false) where T
    # Parameterize the model. First, fix the value of the incoming state
    # variables. Then parameterize the model depending on `noise`. Finally,
    # set the objective. Note that we set the objective every time incase
    # the user calls set_stage_objective in the parameterize function.
    #
    # TODO(odow): cache the need to call set_objective. Only call it if the
    # objective changes.
    set_incoming_state(node, state)
    node.parameterize(noise)
    set_objective(node)
    JuMP.optimize!(node.subproblem)
    # Test for primal feasibility.
    primal_status = JuMP.primal_status(node.subproblem)
    if primal_status != JuMP.MOI.FeasiblePoint
        error("Unable to solve node $(node.index). Primal status: " *
              "$(primal_status).")
    end
    state = get_outgoing_state(node)
    stage_objective = JuMP.result_value(node.stage_objective)
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
    return state, dual_values, stage_objective
end

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
        state, duals, stage_objective = solve_subproblem(graph, node, state, noise)
        cumulative_value += stage_objective
        push!(sampled_states, state)
    end
    return scenario_path, sampled_states, cumulative_value
end

function backward_pass(graph::PolicyGraph{T}, options::Options,
                       scenario_path,
                       sampled_states::Vector{Dict{Symbol, Float64}})
    for index in (length(scenario_path)-1):-1:1
        node_index, noise = scenario_path[index]
        state = sampled_states[index]
        node = graph[node_index]

        noise_supports = []
        original_probability = []
        dual_variables = Dict{Symbol, Float64}[]
        objective_realizations = Float64[]
        for child in node.children
            child_node = graph[child.term]
            for noise in child_node.noise_terms
                outgoing_state, duals, obj = solve_subproblem(graph, child_node,
                                                              state, noise.term)
                push!(dual_variables, duals)
                push!(noise_supports, noise.term)
                push!(original_probability, child.probabiltiy * noise.probability)
                push!(objective_realizations, obj)
            end
        end
        # refine_bellman_function(graph,
        #                         node,
        #                         state,
        #                         dual_variables,
        #                         noise_supports,
        #                         original_probability,
        #                         objective_realizations)
    end
end
# function solve(graph::PolicyGraph;
#                iteration_limit = 100_000,
#                time_limit = Inf,
#                stopping_rules = [],
#                cut_selection = nothing,
#                risk_measure = Kokako.Expectation(),
#                sampling_scheme = Kokako.MonteCarlo(),
#                print_level = 0,
#                log_file = "log.txt",
#                cut_file = "cuts.csv")
#     status = :not_solved
#     # while !convergence_test(graph)
#     #     iteration(graph, options)
#     # end
#     return status
# end
