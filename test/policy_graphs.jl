using Kokako, Test, GLPK

@testset "Linear" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode=false) do subproblem, stage
    end

    @test_throws Exception Kokako.PolicyGraph(Kokako.LinearGraph(2)
                                                  ) do subproblem, stage
    end

    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               optimizer = with_optimizer(GLPK.Optimizer)
                                   ) do subproblem, stage
    end
end

@testset "Markovian" begin
    graph = Kokako.MarkovianGraph([
            ones(Float64, (1, 1)),
            [0.5 0.5],
            [0.5 0.5; 0.3 0.4],
            [0.5 0.5; 0.3 0.4],
            [0.5 0.5; 0.3 0.4]
        ]
    )
    model = Kokako.PolicyGraph(graph, direct_mode=false) do subproblem, stage
    end
end

@testset "General" begin
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
    model = Kokako.PolicyGraph(graph, direct_mode=false) do subproblem, stage
    end
end
