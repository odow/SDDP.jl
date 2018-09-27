using Kokako

# ========== Linear policy graph ==========
graph = Kokako.LinearGraph(stages = 5)
# For infinite horizon:
Kokako.add_edge(graph, (5, 1) => 0.9)

indices = [1, 2, 3, 4, 5]

# ========== Markovian policy graph ==========
graph = Kokako.MarkovianGraph(stages = 5,
    transition_matrix = [
        [1.0],
        [0.5, 0.5],
        [0.5 0.5; 0.3 0.4],
        [0.5 0.5; 0.3 0.4],
        [0.5 0.5; 0.3 0.4]
    ],
    root_node_transition = [1.0],  # Argument not needed.
)
# For infinite horizon:
Kokako.add_edge(graph, ((5, 1), (1, 1)) => 0.9)
Kokako.add_edge(graph, ((5, 2), (1, 1)) => 0.9)

indices = [(1,1), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2)]

graph = Kokako.MarkovianGraph(stages = 5,
    transition_matrix = [0.5 0.5; 0.3 0.4]
    root_node_transition = [0.5, 0.5]
)
# For infinite horizon:
Kokako.add_edge(graph, ((5, 1), (1, 1)) => 0.4)
Kokako.add_edge(graph, ((5, 2), (1, 1)) => 0.4)
Kokako.add_edge(graph, ((5, 1), (1, 2)) => 0.4)
Kokako.add_edge(graph, ((5, 2), (1, 2)) => 0.4)

indices = [(1,1), (1,2), (2,1), (2,2), (3,1), (3,2), (4,1), (4,2), (5,1), (5,2)]

# ========== General graph ==========
graph = Kokako.Graph(
    root_node = :root,
    nodes = [:stage_1, :stage_2, :stage_3],
    edges = [
        (:root, :stage_1) => 1.0,
        (:stage_1, :stage_2) => 1.0,
        (:stage_2, :stage_3) => 1.0,
        (:stage_3, :stage_1) => 0.9
    ]
)

# ========== The model ==========
model = Kokako.PolicyGraph(
            graph = graph,
            root_node_risk_measure = Kokako.Expectation(),
            optimizer = with_optimizer(GLPK.Optimizer),
            # jump_mode = Direct / Manual / Automatic?
            sense = :Min
                ) do subproblem, index

    @state(subproblem, 0 <= x′ <= 1, x==1)
    @state(subproblem, y′, y==2)

    # JuMP definitions
    @variable(subproblem, u >= 0)
    @constraints(subproblem, begin
        con, x′ == x + u
        y′ == y - u
    end)

    # Set a risk measure for the subproblem.
    Kokako.set_risk_measure(subproblem, Kokako.Expectation())

    # If not probability given in 3'rd argument, defaults to uniform.
    Kokako.parameterize(subproblem, Kokako.DiscreteDistribution(
            [1, 2, 3], [0.5, 0.2, 0.3])
            ) do ω
        # We parameterize the JuMP model using JuMP syntax. The tricky thing is
        # what to do with the stage objective...
        JuMP.setRHS(subproblem, con = ω)
        @stageobjective(subproblem, ω * u)
    end
end

Kokako.solve(model,
    iteration_limit = 10,
    time_limt = 20,
    stopping_rules = [Kokako.BoundStalling(), Kokako.Statistical()],
    cut_selection = Kokako.LevelOne(),
    risk_measure = (index) -> index < 5 ? Kokako.Expectation() : Kokako.WorstCase(),
    sampling_scheme = Kokako.Basic(),
    print_level = 0,
    log_file = "log.txt",
    cut_output_file = "cuts.json"
)

Kokako.termination_status(model) == Kokako.IterationLimit

# Single realization of historical (potentially out-of-sample) dataset.
Kokako.simulate(model,
    sampling_scheme = Kokako.Historical(
        nodes = [1, 2, 3, 4, 5],
        noise = [1, 2, 4, 3, 2]
    )
)

# 100 Monte Carlo replications.
Kokako.simulate(model,
    sampling_scheme = Kokako.MonteCarlo(100)
)

# 100 Monte Carlo replications using a different distribution.
Kokako.simulate(model,
    sampling_scheme = Kokako.MonteCarlo(100,
        (index) -> Kokako.DiscreteDistribution(
            [4, 5, 6], [0.4, 0.3, 0.3]
        )
    )
)
