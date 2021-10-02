#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestSamplingSchemes

using SDDP
using Test

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_InSampleMonteCarlo_Acyclic()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
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
    return
end

function test_InSampleMonteCarlo_Cyclic()
    graph = SDDP.LinearGraph(2)
    SDDP.add_edge(graph, 2 => 1, 0.9)
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
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
    return
end

function test_OutOfSampleMonteCarlo_Acyclic()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    @test_throws ErrorException SDDP.OutOfSampleMonteCarlo(
        (node) -> nothing,
        model,
        max_depth = 0,
        terminate_on_dummy_leaf = false,
        terminate_on_cycle = false,
    )
    sampler = SDDP.OutOfSampleMonteCarlo(
        model,
        use_insample_transition = true,
    ) do stage
        return [SDDP.Noise(2 * stage, 0.4), SDDP.Noise(4 * stage, 0.6)]
    end
    scenario, terminated_due_to_cycle = SDDP.sample_scenario(model, sampler)
    @test length(scenario) == 2
    @test !terminated_due_to_cycle
    for (stage, (node, noise)) in enumerate(scenario)
        @test stage == node
        @test noise in stage * [2, 4]
    end
    sampler = SDDP.OutOfSampleMonteCarlo(
        model,
        use_insample_transition = false,
    ) do stage
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
    return
end

function test_OutOfSampleMonteCarlo_Cyclic()
    graph = SDDP.LinearGraph(2)
    SDDP.add_edge(graph, 2 => 1, 0.9)
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
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

function test_Historical()
    @test_throws Exception SDDP.Historical([[1, 2], [3, 4]], [0.6, 0.6])
    return
end

function test_Historical_SingleTrajectory()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    scenario, terminated_due_to_cycle = SDDP.sample_scenario(
        model,
        SDDP.Historical([(1, 0.1), (2, 0.2), (1, 0.3)]),
    )
    @test length(scenario) == 3
    @test !terminated_due_to_cycle
    @test scenario == [(1, 0.1), (2, 0.2), (1, 0.3)]
    return
end

function test_Historical_multiple()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    scenario_A = [(1, 0.1), (2, 0.2), (1, 0.3)]
    scenario_B = [(1, 0.4), (2, 0.5)]
    for i in 1:10
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
    return
end

function test_PSR()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x, ω)
        end
    end
    scheme = SDDP.PSRSamplingScheme(2)
    scenario_1, term_1 = SDDP.sample_scenario(model, scheme)
    @test length(scenario_1) == 2
    @test !term_1
    @test length(scheme.scenarios) == 1
    scenario_2, term_2 = SDDP.sample_scenario(model, scheme)
    @test length(scenario_2) == 2
    @test !term_2
    @test length(scheme.scenarios) == 2
    scenario_3, _ = SDDP.sample_scenario(model, scheme)
    @test scenario_1 == scenario_3
    @test length(scheme.scenarios) == 2
    scenario_4, _ = SDDP.sample_scenario(model, scheme)
    @test scenario_2 == scenario_4
    @test length(scheme.scenarios) == 2
    return
end

end  # module

TestSamplingSchemes.runtests()
