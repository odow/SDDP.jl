using Kokako, Test

@testset "TimeLimit" begin
    graph = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode = false) do node, stage
        @variable(node, x, Kokako.State, initial_value = 0)
    end
    rule = Kokako.TimeLimit(0.5)
    @test Kokako.stopping_rule_status(rule) == :time_limit
    @test Kokako.convergence_test(graph, [Kokako.Log(1, 0.0, 0.0, 1.0)], rule)
    @test !Kokako.convergence_test(graph, [Kokako.Log(1, 0.0, 0.0, 0.1)], rule)
end

@testset "IterationLimit" begin
    graph = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode = false) do node, stage
        @variable(node, x, Kokako.State, initial_value = 0)
    end
    rule = Kokako.IterationLimit(2)
    @test Kokako.stopping_rule_status(rule) == :iteration_limit
    @test Kokako.convergence_test(graph,
        [Kokako.Log(1, 0.0, 0.0, 1.0), Kokako.Log(2, 0.0, 0.0, 1.0)], rule)
    @test !Kokako.convergence_test(graph, [Kokako.Log(1, 0.0, 0.0, 0.1)], rule)
end

@testset "BoundStalling" begin
    graph = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode = false) do node, stage
        @variable(node, x, Kokako.State, initial_value = 0)
    end
    rule = Kokako.BoundStalling(3, 1.0)
    @test Kokako.stopping_rule_status(rule) == :bound_stalling
    @test Kokako.convergence_test(graph, [
        Kokako.Log(1, 0.0, 0.0, 1.0),
        Kokako.Log(2, 0.9, 0.0, 1.0),
        Kokako.Log(3, 1.5, 0.0, 1.0),
        Kokako.Log(4, 1.5, 0.0, 1.0)
    ], rule)
    @test !Kokako.convergence_test(graph, [
        Kokako.Log(1, 0.0, 0.0, 1.0),
        Kokako.Log(2, 1.9, 0.0, 1.0),
        Kokako.Log(3, 2.0, 0.0, 1.0),
        Kokako.Log(4, 2.0, 0.0, 1.0)
    ], rule)
end
