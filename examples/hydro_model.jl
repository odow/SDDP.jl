using Kokako

# ========== Linear policy graph ==========
graph = Kokako.LinearGraph(5)
# For infinite horizon:
Kokako.add_edge(graph, (5 => 1, 0.9)

indices = [1, 2, 3, 4, 5]

# ========== Markovian policy graph ==========
graph = Kokako.MarkovianGraph([
        reshape([1.0], 1, 1),
        [0.5, 0.5],
        [0.5 0.5; 0.3 0.4],
        [0.5 0.5; 0.3 0.4],
        [0.5 0.5; 0.3 0.4]
    ]
)
# For infinite horizon:
Kokako.add_edge(graph, ((5, 1) => (1, 1), 0.9)
Kokako.add_edge(graph, ((5, 2) => (1, 1), 0.9)

indices = [(1,1), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2)]

graph = Kokako.MarkovianGraph(
    stages = 5,
    transition_matrix = [0.5 0.5; 0.3 0.4]
    root_node_transition = [0.5, 0.5]
)
# For infinite horizon:
Kokako.add_edge(graph, ((5, 1) => (1, 1), 0.4)
Kokako.add_edge(graph, ((5, 2) => (1, 1), 0.4)
Kokako.add_edge(graph, ((5, 1) => (1, 2), 0.4)
Kokako.add_edge(graph, ((5, 2) => (1, 2), 0.4)

indices = [(1,1), (1,2), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2)]

# ========== General graph ==========
graph = Kokako.Graph(
    :root,
    [:stage_1, :stage_2, :stage_3],
    [
        (:root => :stage_1, 1.0),
        (:stage_1 => :stage_2, 1.0),
        (:stage_2 => :stage_3, 1.0),
        (:stage_3 => :stage_1, 0.9)
    ]
)

# ========== The model ==========
model = Kokako.PolicyGraph(graph,
            bellman_function = AverageCut(),
            optimizer = with_optimizer(GLPK.Optimizer),
            direct_mode = true
                ) do subproblem, index
    @variables(subproblem, begin
        0 <= x′ <= 1,
        x
        y′
        y
        u >= 0
    end)
    @constraints(subproblem, begin
        con, x′ == x + u
        y′ == y - u
    end)

    Kokako.add_state_variable(subproblem, :x, x, x′)
    Kokako.add_state_variable(subproblem, :y, y, y′)

    # If not probability given in 3'rd argument, defaults to uniform.
    Kokako.parameterize(subproblem, [1, 2, 3], [0.5, 0.2, 0.3]) do ω
        # We parameterize the JuMP model using JuMP syntax.
        # JuMP.setRHS(subproblem, con = ω)
        Kokako.set_stage_objective(subproblem, :Min, ω * u)
    end
end

# Kokako.solve(model,
#     # Stop after K iterations.
#     iteration_limit = 10,
#     # Stop after T seconds.
#     time_limt = 20,
#     # A list of stopping rules.
#     stopping_rules = [Kokako.BoundStalling(), Kokako.Statistical()],
#     # Cut selection techniques to be used at each node.
#     cut_selection = Kokako.LevelOne(),
#     # A risk measure to use at the root node for returning the lower bound.
#     root_node_risk_measure = Kokako.Expectation(),
#     # A risk measure for each node.
#     risk_measure = (index) -> index < 5 ? Kokako.Expectation() : Kokako.WorstCase(),
#     # How to sample on the forward pass.
#     sampling_scheme = Kokako.MonteCarlo(),
#     # Control the level of logging to screen.
#     print_level = 0,
#     # Pipe the log to a file.
#     log_file = "log.txt",
#     # Write the cuts to a file.
#     cut_output_file = "cuts.json"
# )
#
# # Query the termination status. One of:
# #    IterationLimit
# #    TimeLimit
# #    BoundStall
# #    Statistical
# Kokako.termination_status(model) == Kokako.IterationLimit
#
# # Single realization of historical (potentially out-of-sample) dataset.
# Kokako.simulate(model,
#     sampling_scheme = Kokako.Historical(
#         nodes = [1, 2, 3, 4, 5],
#         noise = [1, 2, 4, 3, 2]
#     )
# )
#
# # 100 Monte Carlo replications.
# Kokako.simulate(model,
#     sampling_scheme = Kokako.MonteCarlo(100)
# )
#
# # 100 Monte Carlo replications using a different distribution.
# Kokako.simulate(model,
#     sampling_scheme = Kokako.MonteCarlo(100,
#         (index) -> Kokako.DiscreteDistribution(
#             [4, 5, 6], [0.4, 0.3, 0.3]
#         )
#     )
# )
