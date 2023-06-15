#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestForwardPasses

using SDDP
using Test
import HiGHS
import Random

function runtests()
    for name in names(@__MODULE__, all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_DefaultForwardPass()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    forward_trajectory = SDDP.forward_pass(
        model,
        SDDP.Options(model, Dict(:x => 1.0)),
        SDDP.DefaultForwardPass(),
    )
    simulated_value = 0.0
    for ((node_index, noise), state) in
        zip(forward_trajectory.scenario_path, forward_trajectory.sampled_states)
        @test state[:x] == noise
        simulated_value += noise
    end
    @test simulated_value == forward_trajectory.cumulative_value
    return
end

function test_RevisitingForwardPass()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
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
            SDDP.Options(model, Dict(:x => 1.0); forward_pass = fp),
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
    return
end

function test_RiskAdjustedForwardPass()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
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
    return
end

function test_DefaultForwardPass_cyclic()
    graph = SDDP.LinearGraph(3)
    SDDP.add_edge(graph, 3 => 1, 0.9)
    model = SDDP.PolicyGraph(
        graph;
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    pass = SDDP.DefaultForwardPass()
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        sampling_scheme = SDDP.InSampleMonteCarlo(terminate_on_cycle = true),
        forward_pass = pass,
    )
    forward_trajectory = SDDP.forward_pass(model, options, pass)
    @test length(forward_trajectory.scenario_path) == 4
    @test length(forward_trajectory.sampled_states) == 4
    @test options.starting_states[1] == [Dict(:x => 4.5)]
    @test isempty(options.starting_states[2])
    @test isempty(options.starting_states[3])
    return
end

function test_DefaultForwardPass_cyclic_include_last_node()
    graph = SDDP.LinearGraph(3)
    SDDP.add_edge(graph, 3 => 1, 0.9)
    model = SDDP.PolicyGraph(
        graph;
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    pass = SDDP.DefaultForwardPass(include_last_node = false)
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        sampling_scheme = SDDP.InSampleMonteCarlo(terminate_on_cycle = true),
        forward_pass = pass,
    )
    forward_trajectory = SDDP.forward_pass(model, options, pass)
    @test length(forward_trajectory.scenario_path) == 3
    @test length(forward_trajectory.sampled_states) == 3
    @test options.starting_states[1] == [Dict(:x => 4.5)]
    @test isempty(options.starting_states[2])
    @test isempty(options.starting_states[3])
    return
end

function test_DefaultForwardPass_acyclic_include_last_node()
    graph = SDDP.LinearGraph(3)
    model = SDDP.PolicyGraph(
        graph;
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
    end
    pass = SDDP.DefaultForwardPass(include_last_node = false)
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        sampling_scheme = SDDP.InSampleMonteCarlo(terminate_on_cycle = true),
        forward_pass = pass,
    )
    forward_trajectory = SDDP.forward_pass(model, options, pass)
    @test length(forward_trajectory.scenario_path) == 3
    @test length(forward_trajectory.sampled_states) == 3
    @test isempty(options.starting_states[1])
    @test isempty(options.starting_states[2])
    @test isempty(options.starting_states[3])
    return
end

function test_RegularizedForwardPass()
    function main(capacity_cost, forward_pass, hint)
        Random.seed!(1245)
        graph = SDDP.LinearGraph(2)
        SDDP.add_edge(graph, 2 => 2, 0.95)
        model = SDDP.PolicyGraph(
            graph;
            sense = :Min,
            lower_bound = 0.0,
            optimizer = HiGHS.Optimizer,
        ) do sp, node
            @variable(sp, 0 <= x <= 400, SDDP.State, initial_value = hint)
            @variable(sp, 0 <= y, SDDP.State, initial_value = 0)
            if node == 1
                @stageobjective(sp, capacity_cost * x.out)
                @constraint(sp, y.out == y.in)
            else
                @variable(sp, 0 <= u_prod <= 200)
                @variable(sp, u_overtime >= 0)
                @stageobjective(sp, 100u_prod + 300u_overtime + 50y.out)
                @constraint(sp, x.out == x.in)
                @constraint(sp, y.out <= x.in)
                @constraint(sp, c_bal, y.out == y.in + u_prod + u_overtime)
                SDDP.parameterize(sp, [100, 300]) do ω
                    set_normalized_rhs(c_bal, -ω)
                    return
                end
            end
            return
        end
        SDDP.train(model; print_level = 0, forward_pass = forward_pass)
        results = SDDP.simulate(model, 1, [:x])
        log = model.most_recent_training_results.log
        return results[1][1][:x].out, length(log)
    end
    for (cost, hint) in [(0, 400), (200, 100), (400, 0)]
        fp = SDDP.RegularizedForwardPass()
        reg_capacity, reg_num_iterations = main(cost, fp, hint)
        capacity, num_iterations = main(cost, SDDP.DefaultForwardPass(), hint)
        @test isapprox(reg_capacity, capacity; atol = 1e-2)
        @test reg_num_iterations <= num_iterations
    end
    fp = SDDP.RegularizedForwardPass()
    reg_capacity, reg_num_iterations = main(0, fp, 0)
    capacity, num_iterations = main(0, SDDP.DefaultForwardPass(), 0)
    @test isapprox(reg_capacity, capacity; atol = 1e-2)
    @test reg_num_iterations > num_iterations
    return
end

end  # module

TestForwardPasses.runtests()
