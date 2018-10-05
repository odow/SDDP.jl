using Kokako, GLPK, Test

function infinite_trivial()
    graph = Kokako.Graph(
        :root_node,
        [:week],
        [(:root_node => :week, 1.0), (:week => :week, 0.9)]
    )
    model = Kokako.PolicyGraph(graph,
                bellman_function = Kokako.AverageCut(lower_bound = 0),
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do subproblem, node
        @variable(subproblem, state, Kokako.State, initial_value = 0)
        @constraint(subproblem, state.in == state.out)
        @stageobjective(subproblem, 2.0)
    end
    Kokako.train(model, iteration_limit = 100)
    @test Kokako.calculate_bound(model) ≈ 2.0 / (1 - 0.9) atol = 1e-3
end

infinite_trivial()
