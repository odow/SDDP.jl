#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

module TestStoppingRules

using Random
using SDDP
using Test
import HiGHS

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

function test_TimeLimit()
    graph = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0)
    end
    rule = SDDP.TimeLimit(0.5)
    @test SDDP.stopping_rule_status(rule) == :time_limit
    @test SDDP.convergence_test(
        graph,
        [SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    @test !SDDP.convergence_test(
        graph,
        [SDDP.Log(1, 0.0, 0.0, 0.1, 1, 1, " ", false)],
        rule,
    )
    return
end

function test_IterationLimit()
    graph = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0)
    end
    rule = SDDP.IterationLimit(2)
    @test SDDP.stopping_rule_status(rule) == :iteration_limit
    @test SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 0.0, 0.0, 1.0, 1, 1, " ", false),
        ],
        rule,
    )
    @test !SDDP.convergence_test(
        graph,
        [SDDP.Log(1, 0.0, 0.0, 0.1, 1, 1, " ", false)],
        rule,
    )
    return
end

function test_Statistical()
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
        optimizer = HiGHS.Optimizer,
        sense = :Min,
    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0.0)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_lower_bound(x.out, ω)
        end
        @stageobjective(node, x.out)
    end
    SDDP.train(model, iteration_limit = 1, print_level = 0)
    rule = SDDP.Statistical(num_replications = 20)
    @test SDDP.stopping_rule_status(rule) == :statistical
    Random.seed!(123)
    @test SDDP.convergence_test(
        model,
        [SDDP.Log(1, 6.0, 9.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    @test !SDDP.convergence_test(
        model,
        [SDDP.Log(1, 0.0, 9.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    @test SDDP.convergence_test(
        model,
        [SDDP.Log(1, 12.0, 9.0, 1.0, 1, 1, " ", false)],
        rule,
    )

    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        bellman_function = SDDP.BellmanFunction(upper_bound = 6.0),
        optimizer = HiGHS.Optimizer,
        sense = :Max,
    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0.0)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
        end
        @stageobjective(node, x.out)
    end
    SDDP.train(model, iteration_limit = 1, print_level = 0)
    rule = SDDP.Statistical(num_replications = 20)
    @test SDDP.stopping_rule_status(rule) == :statistical
    Random.seed!(123)
    @test SDDP.convergence_test(
        model,
        [SDDP.Log(1, 6.0, 9.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    @test SDDP.convergence_test(
        model,
        [SDDP.Log(1, 0.0, 9.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    @test !SDDP.convergence_test(
        model,
        [SDDP.Log(1, 12.0, 9.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    return
end

function test_BoundStalling()
    graph = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0)
    end
    rule = SDDP.BoundStalling(3, 1.0)
    @test SDDP.stopping_rule_status(rule) == :bound_stalling
    # Not enough iterations to terminate.
    @test !SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 1.9, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(3, 2.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(4, 2.0, 0.0, 1.0, 1, 1, " ", false),
        ],
        rule,
    )
    # Now there is. But only just...
    @test SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 1.9, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(3, 2.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(4, 2.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(5, 2.9, 0.0, 1.0, 1, 1, " ", false),
        ],
        rule,
    )
    # This also meets the test, but we don't terminate because it hasn't
    # differed from the initial bound.
    @test !SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.1, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 0.0, 0.2, 1.1, 1, 2, " ", false),
            SDDP.Log(3, 0.0, 0.1, 1.2, 1, 3, " ", false),
            SDDP.Log(4, 0.0, 0.0, 1.3, 1, 4, " ", false),
        ],
        rule,
    )
    # This also meets the test, because it looks like a deterministic
    # policy
    @test SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(3, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(4, 0.0, 0.0, 1.0, 1, 1, " ", false),
        ],
        rule,
    )
    return
end

function test_StoppingChain()
    graph = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0)
    end
    rule = SDDP.StoppingChain(SDDP.IterationLimit(2), SDDP.TimeLimit(60.0))
    @test SDDP.stopping_rule_status(rule) ==
          Symbol("iteration_limit ∧ time_limit")
    # Not enough iterations to terminate.
    @test !SDDP.convergence_test(
        graph,
        [SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false)],
        rule,
    )
    # How there is. But not enough time.
    @test !SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 0.0, 0.0, 59.0, 1, 1, " ", false),
        ],
        rule,
    )
    # Both satisfied.
    @test SDDP.convergence_test(
        graph,
        [
            SDDP.Log(1, 0.0, 0.0, 1.0, 1, 1, " ", false),
            SDDP.Log(2, 0.0, 0.0, 59.0, 1, 1, " ", false),
            SDDP.Log(3, 0.0, 0.0, 60.1, 1, 1, " ", false),
        ],
        rule,
    )
    return
end

function test_SimulationStoppingRule()
    graph = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0)
        @stageobjective(node, x.out)
    end
    rule = SDDP.SimulationStoppingRule()
    @test rule.replications == -1
    @test SDDP.stopping_rule_status(rule) == :simulation_stopping
    log = [
        SDDP.Log(1, 0.000000e+00, 8.316000e+03, 1.559195, 1, 14, "", false),
        SDDP.Log(2, 3.171195e+03, 8.767171e+03, 1.801409, 1, 136, "", false),
        SDDP.Log(3, 4.057980e+03, 4.500000e+03, 1.807249, 1, 150, "", false),
        SDDP.Log(4, 4.074139e+03, 2.314272e+03, 1.813528, 1, 164, "", false),
        SDDP.Log(5, 4.074139e+03, 4.716000e+03, 1.819679, 1, 178, "", false),
        SDDP.Log(6, 4.074139e+03, 2.308500e+03, 1.824431, 1, 192, "", false),
        SDDP.Log(7, 4.074139e+03, 2.308500e+03, 1.830817, 1, 206, "", false),
        SDDP.Log(8, 4.074139e+03, 2.308500e+03, 1.837420, 1, 220, "", false),
        SDDP.Log(9, 4.074139e+03, 5.132230e+03, 1.843861, 1, 234, "", false),
        SDDP.Log(10, 4.074139e+03, 5.197500e+03, 1.850351, 1, 248, "", false),
        SDDP.Log(11, 4.074139e+03, 4.716000e+03, 1.856620, 1, 262, "", false),
        SDDP.Log(12, 4.074139e+03, 2.308500e+03, 1.862838, 1, 276, "", false),
        SDDP.Log(13, 4.074139e+03, 2.308500e+03, 1.869224, 1, 290, "", false),
        SDDP.Log(14, 4.074139e+03, 2.308500e+03, 1.875853, 1, 304, "", false),
        SDDP.Log(15, 4.074139e+03, 2.308500e+03, 1.882504, 1, 318, "", false),
        SDDP.Log(16, 4.074139e+03, 5.197500e+03, 1.889759, 1, 332, "", false),
        SDDP.Log(17, 4.074139e+03, 5.132230e+03, 1.896462, 1, 346, "", false),
        SDDP.Log(18, 4.074139e+03, 8.086500e+03, 1.903102, 1, 360, "", false),
        SDDP.Log(19, 4.074139e+03, 2.308500e+03, 1.910075, 1, 374, "", false),
        SDDP.Log(20, 4.074139e+03, 5.132230e+03, 1.917460, 1, 388, "", false),
    ]
    @test !SDDP.convergence_test(graph, log[1:1], rule)
    @test rule.replications == 1
    @test !SDDP.convergence_test(graph, log[1:4], rule)
    @test !SDDP.convergence_test(graph, log[1:10], rule)
    @test !SDDP.convergence_test(graph, log[1:19], rule)
    @test SDDP.convergence_test(graph, log[1:20], rule)
    return
end

function test_FirstStageStoppingRule()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0)
        @stageobjective(node, x.out)
    end
    rule = SDDP.FirstStageStoppingRule()
    SDDP.train(model; stopping_rules = [rule])
    @test length(rule.data) == 50
    log = model.most_recent_training_results.log
    set_lower_bound(model[1].subproblem[:x].out, 1.0)
    @test !SDDP.convergence_test(model, log, rule)
    @test length(rule.data) == 51
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0)
        SDDP.parameterize(node, 1:2) do w
            return @stageobjective(node, w * x.out)
        end
        return
    end
    @test_throws(
        ErrorException(
            "FirstStageStoppingRule cannot be applied because first-stage is " *
            "not deterministic",
        ),
        SDDP.train(
            model;
            print_level = 0,
            stopping_rules = [SDDP.FirstStageStoppingRule()],
        ),
    )
    graph = SDDP.Graph(0, [1, 2], [(0 => 1, 0.5), (0 => 2, 0.5)])
    model = SDDP.PolicyGraph(
        graph;
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x >= 0, SDDP.State, initial_value = 0)
        @stageobjective(node, x.out)
        return
    end
    @test_throws(
        ErrorException(
            "FirstStageStoppingRule cannot be applied because first-stage is " *
            "not deterministic",
        ),
        SDDP.train(
            model;
            print_level = 0,
            stopping_rules = [SDDP.FirstStageStoppingRule()],
        ),
    )
    return
end

end  # module

TestStoppingRules.runtests()
