using Kokako, Test

@testset "MonteCarlo" begin
    @testset "Acyclic" begin
        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   direct_mode=false) do node, stage
            @variable(node, 0 <= x <= 1)
            Kokako.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        scenario = Kokako.sample_scenario(Kokako.MonteCarlo(), model)
        @test length(scenario) == 2
        for (stage, (node, noise)) in enumerate(scenario)
            @test stage == node
            @test noise in stage * [1, 3]
        end
    end

    @testset "Cyclic" begin
        graph = Kokako.LinearGraph(2)
        Kokako.add_edge(graph, 2=>1, 0.9)
        model = Kokako.PolicyGraph(graph, direct_mode=false) do node, stage
            @variable(node, 0 <= x <= 1)
            Kokako.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        scenario = Kokako.sample_scenario(Kokako.MonteCarlo(), model,
                                          max_cycles = 2)
        @test length(scenario) == 4
        for (index, (node, noise)) in enumerate(scenario)
            stage = (index - 1) % 2 + 1
            @test stage == node
            @test noise in stage * [1, 3]
        end
    end
end
