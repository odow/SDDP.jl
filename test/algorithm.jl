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
            SDDP.CompleteSampler(),
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

@testset "refine_at_similar_nodes" begin
    model = SDDP.MarkovianPolicyGraph(
        transition_matrices = [ [0.5 0.5], [0.2 0.8; 0.8 0.2]],
        optimizer = with_optimizer(GLPK.Optimizer),
        lower_bound = 0.0
    ) do sp, index
        stage, markov_state = index
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.out >= stage)
        @stageobjective(sp, (stage + markov_state) * x.out)
    end
    SDDP.train(model, iteration_limit = 1, refine_at_similar_nodes = false, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 5.7 || SDDP.calculate_bound(model) ≈ 6.3
    mi1 = length(model[(1,1)].bellman_function.global_theta.cut_oracle.cuts)
    mi2 = length(model[(1,2)].bellman_function.global_theta.cut_oracle.cuts)
    @test mi1 + mi2 == 1

    model = SDDP.MarkovianPolicyGraph(
        transition_matrices = [ [0.5 0.5], [0.2 0.8; 0.8 0.2]],
        optimizer = with_optimizer(GLPK.Optimizer),
        lower_bound = 0.0
    ) do sp, index
        stage, markov_state = index
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.out >= stage)
        @stageobjective(sp, (stage + markov_state) * x.out)
    end
    SDDP.train(model, iteration_limit = 1, refine_at_similar_nodes = true, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 9.5
    @test length(model[(1, 1)].bellman_function.global_theta.cut_oracle.cuts) == 1
    @test length(model[(1, 2)].bellman_function.global_theta.cut_oracle.cuts) == 1
end

@testset "optimize_hook" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        optimizer = with_optimizer(GLPK.Optimizer),
        lower_bound = 0.0
    ) do  sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0)
        @stageobjective(sp, x.out)
    end
    pre_optimize_called = 0
    post_optimize_called = 0
    node = model[1]
    SDDP.pre_optimize_hook(node) do model, node, state, noise, scenario_path, require_duals
        pre_optimize_called = 1
        return pre_optimize_called
    end
    SDDP.post_optimize_hook(node) do ret
        post_optimize_called = ret + 2
        return
    end
    SDDP.solve_subproblem(
        model, node, Dict(:x => 0.0), nothing, Tuple{Int, Any}[(1, nothing)];
        require_duals = false
    )
    @test pre_optimize_called == 1
    @test post_optimize_called == 3
end

@testset "relax_integrality and enforce_integrality" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2, lower_bound = 0.0,
        direct_mode = false
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 2.0)
        @variable(sp, b1, Bin)
        @variable(sp, 0.2 <= b2, Bin)
        @variable(sp, 0.5 <= b3 <= 1.2, Bin)
        @variable(sp, i1, Int)
        @variable(sp, 6.2 >= i2, Int)
        @variable(sp, -8 <= i3 <= 2, Int)
        @stageobjective(sp, b1 + b2 + b2 + i3 + i1)
    end

    for node in [model[1], model[2]]
        @test JuMP.is_binary(node.subproblem[:b1])
        @test !JuMP.has_lower_bound(node.subproblem[:b1])
        @test !JuMP.has_upper_bound(node.subproblem[:b1])

        @test JuMP.is_binary(node.subproblem[:b2])
        @test JuMP.lower_bound(node.subproblem[:b2]) == 0.2
        @test !JuMP.has_upper_bound(node.subproblem[:b2])

        @test JuMP.is_binary(node.subproblem[:b3])
        @test JuMP.lower_bound(node.subproblem[:b3]) == 0.5
        @test JuMP.upper_bound(node.subproblem[:b3]) == 1.2

        @test JuMP.is_integer(node.subproblem[:i1])
        @test !JuMP.has_lower_bound(node.subproblem[:i1])
        @test !JuMP.has_upper_bound(node.subproblem[:i1])

        @test JuMP.is_integer(node.subproblem[:i2])
        @test JuMP.upper_bound(node.subproblem[:i2]) == 6.2
        @test !JuMP.has_lower_bound(node.subproblem[:i2])

        @test JuMP.is_integer(node.subproblem[:i3])
        @test JuMP.lower_bound(node.subproblem[:i3]) == -8
        @test JuMP.upper_bound(node.subproblem[:i3]) == 2
    end
    binaries, integers = SDDP.relax_integrality(model)
    for node in [model[1], model[2]]
        @test !JuMP.is_binary(node.subproblem[:b1])
        @test JuMP.lower_bound(node.subproblem[:b1]) == 0.0
        @test JuMP.upper_bound(node.subproblem[:b1]) == 1.0

        @test !JuMP.is_binary(node.subproblem[:b2])
        @test JuMP.lower_bound(node.subproblem[:b2]) == 0.2
        @test JuMP.upper_bound(node.subproblem[:b2]) == 1.0

        @test !JuMP.is_binary(node.subproblem[:b3])
        @test JuMP.lower_bound(node.subproblem[:b3]) == 0.5
        @test JuMP.upper_bound(node.subproblem[:b3]) == 1.0

        @test !JuMP.is_integer(node.subproblem[:i1])
        @test !JuMP.has_lower_bound(node.subproblem[:i1])
        @test !JuMP.has_upper_bound(node.subproblem[:i1])

        @test !JuMP.is_integer(node.subproblem[:i2])
        @test JuMP.upper_bound(node.subproblem[:i2]) == 6.2
        @test !JuMP.has_lower_bound(node.subproblem[:i2])

        @test !JuMP.is_integer(node.subproblem[:i3])
        @test JuMP.lower_bound(node.subproblem[:i3]) == -8
        @test JuMP.upper_bound(node.subproblem[:i3]) == 2
    end
    SDDP.enforce_integrality(binaries, integers)
    for node in [model[1], model[2]]
        @test JuMP.is_binary(node.subproblem[:b1])
        @test !JuMP.has_lower_bound(node.subproblem[:b1])
        @test !JuMP.has_upper_bound(node.subproblem[:b1])

        @test JuMP.is_binary(node.subproblem[:b2])
        @test JuMP.lower_bound(node.subproblem[:b2]) == 0.2
        @test !JuMP.has_upper_bound(node.subproblem[:b2])

        @test JuMP.is_binary(node.subproblem[:b3])
        @test JuMP.lower_bound(node.subproblem[:b3]) == 0.5
        @test JuMP.upper_bound(node.subproblem[:b3]) == 1.2

        @test JuMP.is_integer(node.subproblem[:i1])
        @test !JuMP.has_lower_bound(node.subproblem[:i1])
        @test !JuMP.has_upper_bound(node.subproblem[:i1])

        @test JuMP.is_integer(node.subproblem[:i2])
        @test JuMP.upper_bound(node.subproblem[:i2]) == 6.2
        @test !JuMP.has_lower_bound(node.subproblem[:i2])

        @test JuMP.is_integer(node.subproblem[:i3])
        @test JuMP.lower_bound(node.subproblem[:i3]) == -8
        @test JuMP.upper_bound(node.subproblem[:i3]) == 2
    end
end
