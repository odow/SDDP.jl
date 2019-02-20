#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

using Kokako, GLPK, Random, Test

@testset "TimeLimit" begin
    graph = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               lower_bound = 0.0,
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
                               lower_bound = 0.0,
                               direct_mode = false) do node, stage
        @variable(node, x, Kokako.State, initial_value = 0)
    end
    rule = Kokako.IterationLimit(2)
    @test Kokako.stopping_rule_status(rule) == :iteration_limit
    @test Kokako.convergence_test(graph,
        [Kokako.Log(1, 0.0, 0.0, 1.0), Kokako.Log(2, 0.0, 0.0, 1.0)], rule)
    @test !Kokako.convergence_test(graph, [Kokako.Log(1, 0.0, 0.0, 0.1)], rule)
end

@testset "Statistical" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               bellman_function = Kokako.AverageCut(lower_bound = 0.0),
                               optimizer = with_optimizer(GLPK.Optimizer),
                               sense = :Min
                               ) do node, stage
        @variable(node, x >= 0, Kokako.State, initial_value = 0.0)
        Kokako.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_lower_bound(x.out, ω)
        end
        @stageobjective(node, x.out)
    end
    Kokako.train(model, iteration_limit = 1)
    rule = Kokako.Statistical(num_replications = 20)
    @test Kokako.stopping_rule_status(rule) == :statistical
    Random.seed!(123)
    @test Kokako.convergence_test(model, [Kokako.Log(1, 6.0, 9.0, 1.0)], rule)
    @test !Kokako.convergence_test(model, [Kokako.Log(1, 0.0, 9.0, 1.0)], rule)
    @test !Kokako.convergence_test(model, [Kokako.Log(1, 12.0, 9.0, 1.0)], rule)
end

@testset "BoundStalling" begin
    graph = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               lower_bound = 0.0,
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
