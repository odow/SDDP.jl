#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import GLPK
using SDDP
using Test

@testset "Forward Pass" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = GLPK.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    forward_trajectory = SDDP.forward_pass(
        model,
        SDDP.Options(
            model,
            Dict(:x => 1.0),
            SDDP.InSampleMonteCarlo(),
            SDDP.CompleteSampler(),
            SDDP.Expectation(),
            0.0,
            true,
            SDDP.AbstractStoppingRule[],
            (a, b) -> nothing,
            0,
            0.0,
            SDDP.Log[],
            IOBuffer(),
            1,
            SDDP.DefaultForwardPass(),
            SDDP.ContinuousConicDuality(),
            x -> nothing,
        ),
        SDDP.DefaultForwardPass(),
    )
    simulated_value = 0.0
    for ((node_index, noise), state) in
        zip(forward_trajectory.scenario_path, forward_trajectory.sampled_states)
        @test state[:x] == noise
        simulated_value += noise
    end
    @test simulated_value == forward_trajectory.cumulative_value
end

@testset "RevisitingForwardPass" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = GLPK.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    fp = SDDP.RevisitingForwardPass(2; sub_pass = SDDP.DefaultForwardPass())
    @test length(fp.archive) == 0
    for i in 1:5
        pass = SDDP.forward_pass(
            model,
            SDDP.Options(
                model,
                Dict(:x => 1.0),
                SDDP.InSampleMonteCarlo(),
                SDDP.CompleteSampler(),
                SDDP.Expectation(),
                0.0,
                true,
                SDDP.AbstractStoppingRule[],
                (a, b) -> nothing,
                0,
                0.0,
                SDDP.Log[],
                IOBuffer(),
                1,
                fp,
                SDDP.ContinuousConicDuality(),
                x -> nothing,
            ),
            fp,
        )
        if i <= 2
            @test length(fp.archive) == i
        elseif i == 3
            @test length(fp.archive) == 2
            @test pass.cumulative_value == fp.archive[1].cumulative_value
        elseif i == 4
            @test length(fp.archive) == 2
            @test pass.cumulative_value == fp.archive[2].cumulative_value
        elseif i == 5
            @test length(fp.archive) == 3
        end
    end
end

@testset "RiskAdjustedForwardPass" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = GLPK.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    @test_throws ArgumentError SDDP.train(
        model;
        iteration_limit = 5,
        forward_pass_resampling_probability = 0.0,
    )
    @test_throws ArgumentError SDDP.train(
        model;
        iteration_limit = 5,
        forward_pass_resampling_probability = 1.0,
    )

    forward_pass = SDDP.RiskAdjustedForwardPass(
        forward_pass = SDDP.DefaultForwardPass(),
        risk_measure = SDDP.WorstCase(),
        resampling_probability = 0.9,
    )
    SDDP.train(
        model;
        print_level = 0,
        iteration_limit = 20,
        forward_pass = forward_pass,
    )
    @test length(forward_pass.archive) < 10
    SDDP.train(
        model;
        iteration_limit = 10,
        print_level = 0,
        forward_pass_resampling_probability = 0.9,
    )
    @test SDDP.termination_status(model) == :iteration_limit
end
