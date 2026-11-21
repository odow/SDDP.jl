# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module TestForwardPasses

using SDDP
using Test
import HiGHS
import Random

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

function test_DefaultForwardPass()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_upper_bound(x.out, ω)
            return
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
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_upper_bound(x.out, ω)
            return
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
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_upper_bound(x.out, ω)
            return
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

    forward_pass = SDDP.RiskAdjustedForwardPass(;
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
            JuMP.set_upper_bound(x.out, ω)
            return
        end
    end
    pass = SDDP.DefaultForwardPass()
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        sampling_scheme = SDDP.InSampleMonteCarlo(; terminate_on_cycle = true),
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
            JuMP.set_upper_bound(x.out, ω)
            return
        end
    end
    pass = SDDP.DefaultForwardPass(; include_last_node = false)
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        sampling_scheme = SDDP.InSampleMonteCarlo(; terminate_on_cycle = true),
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
            JuMP.set_upper_bound(x.out, ω)
            return
        end
    end
    pass = SDDP.DefaultForwardPass(; include_last_node = false)
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        sampling_scheme = SDDP.InSampleMonteCarlo(; terminate_on_cycle = true),
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
        SDDP.train(
            model;
            print_level = 0,
            forward_pass = forward_pass,
            iteration_limit = 10,
            parallel_scheme = SDDP.Serial(),
        )
        return SDDP.calculate_bound(model)
    end
    for (cost, hint) in [(0, 400), (200, 100), (400, 0)]
        fp = SDDP.RegularizedForwardPass()
        reg_bound = main(cost, fp, hint)
        bound = main(cost, SDDP.DefaultForwardPass(), hint)
        @test reg_bound >= bound - 1.0
    end
    # Test that initializing with a bad guess performs poorly
    fp = SDDP.RegularizedForwardPass()
    reg_bound = main(400, fp, 400)
    bound = main(400, SDDP.DefaultForwardPass(), 0)
    @test reg_bound < bound
    return
end

function test_ImportanceSamplingForwardPass()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_upper_bound(x.out, ω)
            return
        end
    end
    forward_trajectory = SDDP.forward_pass(
        model,
        SDDP.Options(model, Dict(:x => 1.0)),
        SDDP.ImportanceSamplingForwardPass(),
    )
    simulated_value = 0.0
    for ((node_index, noise), state) in
        zip(forward_trajectory.scenario_path, forward_trajectory.sampled_states)
        @test state[:x] == noise
        simulated_value += noise
    end
    @test simulated_value == forward_trajectory.cumulative_value
    SDDP.train(model; cut_type = SDDP.MULTI_CUT, iteration_limit = 1)
    forward_trajectory = SDDP.forward_pass(
        model,
        SDDP.Options(model, Dict(:x => 1.0)),
        SDDP.ImportanceSamplingForwardPass(),
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

function _risk_biased_model(sense)
    return SDDP.LinearPolicyGraph(;
        stages = 4,
        sense = sense,
        lower_bound = sense == :Min ? 0.0 : -100.0,
        upper_bound = sense == :Min ? 100.0 : 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.out == x.in)
        SDDP.parameterize(sp, [1.0, 2.0, 3.0], [0.4, 0.4, 0.2]) do ω
            # Cost of ω when minimizing, profit of -ω when maximizing. In both
            # cases, ω = 3 is the unfavorable outcome.
            @stageobjective(sp, sense == :Min ? ω : -ω)
            return
        end
    end
end

function test_RiskBiasedForwardPass()
    for sense in (:Min, :Max)
        model = _risk_biased_model(sense)
        options = SDDP.Options(model, Dict(:x => 0.0))
        pass = SDDP.RiskBiasedForwardPass()
        # Before training, no risk-adjusted probabilities exist, so we sample
        # from the nominal distribution.
        trajectory = SDDP.forward_pass(model, options, pass)
        @test length(trajectory.scenario_path) == 4
        @test trajectory.cumulative_value ≈
              sum(ω for (_, ω) in trajectory.scenario_path) *
              (sense == :Min ? 1 : -1)
        n_iterations = 5
        SDDP.train(
            model;
            risk_measure = SDDP.WorstCase(),
            iteration_limit = n_iterations,
            print_level = 0,
        )
        # Each iteration adds one cut to each of nodes 1, 2, and 3. WorstCase
        # puts all of the risk-adjusted probability on ω = 3.
        q = [0.4, 0.4, 0.2]
        r = [0.0, 0.0, 1.0]
        p = copy(q)
        for _ in 1:n_iterations
            p = (p .+ r .+ q) ./ 3
        end
        for t in 1:3
            terms = model[t].risk_adjusted_transitions
            @test [x.term for x in terms] == [(t + 1, ω) for ω in 1.0:3.0]
            @test [x.probability for x in terms] ≈ p
            @test all(x.probability > 0 for x in terms)
        end
        @test isempty(model[4].risk_adjusted_transitions)
        # The forward pass should sample ω = 3 with probability p[3], instead of
        # the nominal probability q[3].
        n_worst, n_total = 0, 0
        for _ in 1:2_000
            trajectory = SDDP.forward_pass(model, options, pass)
            @test length(trajectory.scenario_path) == 4
            for (_, ω) in trajectory.scenario_path[2:end]
                n_worst += ω == 3.0
                n_total += 1
            end
        end
        @test isapprox(n_worst / n_total, p[3]; atol = 0.03)
    end
    return
end

function test_RiskBiasedForwardPass_Expectation()
    model = _risk_biased_model(:Min)
    SDDP.train(
        model;
        forward_pass = SDDP.RiskBiasedForwardPass(),
        risk_measure = SDDP.Expectation(),
        iteration_limit = 3,
        print_level = 0,
    )
    node = model[1]
    terms = node.risk_adjusted_transitions
    @test [t.term for t in terms] == [(2, 1.0), (2, 2.0), (2, 3.0)]
    @test [t.probability for t in terms] ≈ [0.4, 0.4, 0.2]
    @test isempty(model[4].risk_adjusted_transitions)
    return
end

function test_RiskBiasedForwardPass_cyclic()
    model = SDDP.PolicyGraph(
        SDDP.UnicyclicGraph(0.5; num_nodes = 2);
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 0.0)
        @constraint(sp, x.out == x.in)
        SDDP.parameterize(sp, [1.0, 2.0]) do ω
            @stageobjective(sp, ω)
            return
        end
    end
    SDDP.train(
        model;
        forward_pass = SDDP.RiskBiasedForwardPass(),
        risk_measure = SDDP.WorstCase(),
        iteration_limit = 20,
        print_level = 0,
    )
    # Nodes 1 and 2 each get one cut per iteration. WorstCase puts all of the
    # risk-adjusted probability on ω = 2, which is rescaled to sum to the
    # probability of the edge to the child, 0.5 for node 2 and 1.0 for node 1.
    for (t, child, edge) in ((1, 2, 1.0), (2, 1, 0.5))
        q = edge .* [0.5, 0.5]
        r = edge .* [0.0, 1.0]
        p = copy(q)
        for _ in 1:20
            p = (p .+ r .+ q) ./ 3
        end
        terms = model[t].risk_adjusted_transitions
        @test [x.term for x in terms] == [(child, 1.0), (child, 2.0)]
        @test [x.probability for x in terms] ≈ p
    end
    # Also check that we respect terminate_on_cycle and max_depth.
    for sampling_scheme in (
        SDDP.InSampleMonteCarlo(; terminate_on_cycle = true),
        SDDP.InSampleMonteCarlo(;
            max_depth = 5,
            terminate_on_dummy_leaf = false,
        ),
    )
        options = SDDP.Options(model, Dict(:x => 0.0); sampling_scheme)
        trajectory =
            SDDP.forward_pass(model, options, SDDP.RiskBiasedForwardPass())
        @test all(ω in (1.0, 2.0) for (_, ω) in trajectory.scenario_path)
        if sampling_scheme.max_depth == 5
            @test length(trajectory.scenario_path) == 5
        else
            @test length(trajectory.scenario_path) <= 3
        end
    end
    return
end

function test_RiskBiasedForwardPass_MonteCarloSampler()
    model = _risk_biased_model(:Min)
    SDDP.train(
        model;
        forward_pass = SDDP.RiskBiasedForwardPass(),
        backward_sampling_scheme = SDDP.MonteCarloSampler(1),
        risk_measure = SDDP.WorstCase(),
        iteration_limit = 5,
        print_level = 0,
    )
    # Even though the backward pass evaluated only one noise per child, every
    # nominal (child, noise) pair has a positive sampling probability.
    for t in 1:3
        terms = model[t].risk_adjusted_transitions
        @test sort([x.term for x in terms]) == [(t + 1, ω) for ω in 1.0:3.0]
        @test all(x.probability > 0 for x in terms)
    end
    return
end

function test_RiskBiasedForwardPass_errors()
    model = _risk_biased_model(:Min)
    sampling_scheme =
        SDDP.Historical([[(1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)]])
    options = SDDP.Options(model, Dict(:x => 0.0); sampling_scheme)
    @test_throws(
        ErrorException(
            "RiskBiasedForwardPass requires the sampling scheme to be an " *
            "`InSampleMonteCarlo`. Got: $(typeof(sampling_scheme))",
        ),
        SDDP.forward_pass(model, options, SDDP.RiskBiasedForwardPass()),
    )
    return
end

function test_LoggingForwardPass()
    model = SDDP.PolicyGraph(
        SDDP.UnicyclicGraph(0.5; num_nodes = 5);
        sense = :Max,
        upper_bound = 100.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 3.0)
        @variable(node, y, SDDP.State, initial_value = 3.0)
        @stageobjective(node, x.out + 2 * y.out)
        @constraint(node, x.out + y.out <= x.in + y.in)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            set_upper_bound(x.out, ω)
            set_upper_bound(y.out, 4 - ω)
            return
        end
    end
    filename = joinpath(mktempdir(), "states.csv")
    SDDP.train(model; forward_pass = SDDP.LoggingForwardPass(; filename))
    contents = read(filename, String)
    @test occursin("iteration,index,x,y\n", contents)
    for iteration in 1:length(model.most_recent_training_results.log)
        for t in 1:5
            @test occursin("$iteration,$t,", contents)
        end
    end
    return
end

end  # module

TestForwardPasses.runtests()
