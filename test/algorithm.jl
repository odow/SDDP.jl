#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Test, GLPK

@testset "Forward Pass" begin
    model = SDDP.PolicyGraph(SDDP.LinearGraph(2);
                sense = :Max,
                bellman_function = SDDP.BellmanFunction(upper_bound = 100.0),
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_upper_bound(x.out, ω)
        end
    end
    forward_trajectory = SDDP.forward_pass(
        model,
        SDDP.Options(
            model,
            Dict(:x => 1.0),
            SDDP.InSampleMonteCarlo(),
            SDDP.Expectation(),
            0.0,
            true
        )
    )
    simulated_value = 0.0
    for ((node_index, noise), state) in zip(forward_trajectory.scenario_path, forward_trajectory.sampled_states)
        @test state[:x] == noise
        simulated_value += noise
    end
    @test simulated_value == forward_trajectory.cumulative_value
end

@testset "to nodal forms" begin
    model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_lower_bound(x.out, ω)
        end
    end
    SDDP.train(model; iteration_limit = 1, print_level = 0,
        risk_measure = SDDP.Expectation())
    @test SDDP.termination_status(model) == :iteration_limit
    SDDP.train(model; iteration_limit = 1, print_level = 0,
        risk_measure = Dict(1 => SDDP.Expectation(), 2 => SDDP.WorstCase()))
    @test SDDP.termination_status(model) == :iteration_limit
    SDDP.train(model; iteration_limit = 1, print_level = 0,
        risk_measure = (idx) -> idx == 1 ? SDDP.Expectation() : SDDP.WorstCase()
    )
    @test SDDP.termination_status(model) == :iteration_limit
end

@testset "solve" begin
    model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_lower_bound(x.out, ω)
        end
    end
    SDDP.train(model; iteration_limit = 4, print_level = 0)
    @test SDDP.termination_status(model) == :iteration_limit
    @testset "simulate" begin
        simulations = SDDP.simulate(model, 11, [:x])
        @test length(simulations) == 11
        @test all(length.(simulations) .== 2)

        simulation = simulations[1][1]
        @test length(keys(simulation)) == 7
        @test sort(collect(keys(simulation))) ==
            [:belief, :bellman_term, :node_index, :noise_term, :objective_state,
            :stage_objective, :x]
        @test typeof(simulation[:x]) == SDDP.State{Float64}
    end
end

@testset "simulate" begin
    model = SDDP.LinearPolicyGraph(
            stages=2, lower_bound=0.0, optimizer=with_optimizer(GLPK.Optimizer)
            ) do sp, t
        @variable(sp, x[i=1:2] >= i, SDDP.State, initial_value = 2i)
        @stageobjective(sp, x[1].out + x[2].out)
    end
    simulations = SDDP.simulate(model, 1, [:x])
    @test simulations[1][1][:x] == [SDDP.State(2.0, 1.0), SDDP.State(4.0, 2.0)]
end

@testset "infeasible model" begin
    model = SDDP.LinearPolicyGraph(
                stages = 2,
                lower_bound = 0.0,
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(node, x.out <= -1)
        @stageobjective(node, x.out)
    end
    @test_throws Exception SDDP.train(model; iteration_limit = 1, print_level = 0)
    @test isfile("subproblem.mps")
    rm("subproblem.mps")
    @test isfile("subproblem.lp")
    rm("subproblem.lp")
end
