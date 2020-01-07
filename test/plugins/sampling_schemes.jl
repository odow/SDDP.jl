#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
using Test

@testset "InSampleMonteCarlo" begin
    @testset "Acyclic" begin
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
            @variable(node, 0 <= x <= 1)
            SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        @test_throws ErrorException SDDP.InSampleMonteCarlo(
            max_depth = 0,
            terminate_on_dummy_leaf = false,
            terminate_on_cycle = false,
        )
        scenario, terminated_due_to_cycle =
            SDDP.sample_scenario(model, SDDP.InSampleMonteCarlo())
        @test length(scenario) == 2
        @test !terminated_due_to_cycle
        for (stage, (node, noise)) in enumerate(scenario)
            @test stage == node
            @test noise in stage * [1, 3]
        end
    end

    @testset "Cyclic" begin
        graph = SDDP.LinearGraph(2)
        SDDP.add_edge(graph, 2 => 1, 0.9)
        model =
            SDDP.PolicyGraph(graph, lower_bound = 0.0, direct_mode = false) do node, stage
                @variable(node, 0 <= x <= 1)
                SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                    JuMP.set_upper_bound(x, ω)
                end
            end
        scenario, terminated_due_to_cycle = SDDP.sample_scenario(
            model,
            SDDP.InSampleMonteCarlo(terminate_on_dummy_leaf = false, max_depth = 4),
        )
        @test length(scenario) == 4
        @test !terminated_due_to_cycle  # Terminated due to max depth.
        for (index, (node, noise)) in enumerate(scenario)
            stage = (index - 1) % 2 + 1
            @test stage == node
            @test noise in stage * [1, 3]
        end
    end
end

@testset "OutOfSampleMonteCarlo" begin
    @testset "Acyclic" begin
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
            @variable(node, 0 <= x <= 1)
            SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        @test_throws ErrorException SDDP.OutOfSampleMonteCarlo(
            (node) -> nothing,
            model,
            max_depth = 0,
            terminate_on_dummy_leaf = false,
            terminate_on_cycle = false,
        )
        sampler = SDDP.OutOfSampleMonteCarlo(model, use_insample_transition = true) do stage
            return [SDDP.Noise(2 * stage, 0.4), SDDP.Noise(4 * stage, 0.6)]
        end
        scenario, terminated_due_to_cycle = SDDP.sample_scenario(model, sampler)
        @test length(scenario) == 2
        @test !terminated_due_to_cycle
        for (stage, (node, noise)) in enumerate(scenario)
            @test stage == node
            @test noise in stage * [2, 4]
        end
        sampler =
            SDDP.OutOfSampleMonteCarlo(model, use_insample_transition = false) do stage
                if stage == 0
                    return [SDDP.Noise(2, 1.0)]
                else
                    return SDDP.Noise{Int}[],
                    [SDDP.Noise(2 * stage, 0.4), SDDP.Noise(4 * stage, 0.6)]
                end
            end
        scenario, terminated_due_to_cycle = SDDP.sample_scenario(model, sampler)
        @test length(scenario) == 1
        @test !terminated_due_to_cycle
        node, noise = scenario[1]
        @test node == 2
        @test noise in [4, 8]
    end

    @testset "Cyclic" begin
        graph = SDDP.LinearGraph(2)
        SDDP.add_edge(graph, 2 => 1, 0.9)
        model =
            SDDP.PolicyGraph(graph, lower_bound = 0.0, direct_mode = false) do node, stage
                @variable(node, 0 <= x <= 1)
                SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                    JuMP.set_upper_bound(x, ω)
                end
            end
        sampler = SDDP.OutOfSampleMonteCarlo(
            model,
            use_insample_transition = true,
            terminate_on_dummy_leaf = false,
            max_depth = 4,
        ) do stage
            return [SDDP.Noise(2 * stage, 0.4), SDDP.Noise(4 * stage, 0.6)]
        end
        scenario, terminated_due_to_cycle = SDDP.sample_scenario(model, sampler)
        @test length(scenario) == 4
        @test !terminated_due_to_cycle  # Terminated due to max depth.
        for (index, (node, noise)) in enumerate(scenario)
            stage = (index - 1) % 2 + 1
            @test stage == node
            @test noise in stage * [2, 4]
        end
    end
end

@testset "Historical" begin
    @test_throws Exception SDDP.Historical([[1, 2], [3, 4]], [0.6, 0.6])
    @testset "Single trajectory" begin
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
            @variable(node, 0 <= x <= 1)
            SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        scenario, terminated_due_to_cycle =
            SDDP.sample_scenario(model, SDDP.Historical([(1, 0.1), (2, 0.2), (1, 0.3)]))
        @test length(scenario) == 3
        @test !terminated_due_to_cycle
        @test scenario == [(1, 0.1), (2, 0.2), (1, 0.3)]
    end
    @testset "Multiple historical" begin
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
            @variable(node, 0 <= x <= 1)
            SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        scenario_A = [(1, 0.1), (2, 0.2), (1, 0.3)]
        scenario_B = [(1, 0.4), (2, 0.5)]
        for i = 1:10
            scenario, terminated_due_to_cycle = SDDP.sample_scenario(
                model,
                SDDP.Historical([scenario_A, scenario_B], [0.2, 0.8]),
            )
            if length(scenario) == 3
                @test scenario == scenario_A
            else
                @test length(scenario) == 2
                @test scenario == scenario_B
            end
            @test !terminated_due_to_cycle
        end
    end
end
