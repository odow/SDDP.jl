#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestValueFunctions

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

function test_ValueFunction_Min()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @constraint(sp, x.out == x.in)
        @stageobjective(sp, 2 * x.out)
    end
    V1 = SDDP.ValueFunction(model[1])
    @test SDDP.evaluate(V1, Dict(:x => 1.0)) == (0.0, Dict(:x => 0.0))
    SDDP.train(model; iteration_limit = 2, print_level = 0)
    V1 = SDDP.ValueFunction(model[1])
    @test sprint(show, V1) == "A value function for node 1"
    for (xhat, yhat, pihat) in
        [(0.0, 0.0, 0.0), (1.0, 2.0, 2.0), (2.0, 4.0, 2.0)]
        @test SDDP.evaluate(V1, Dict(:x => xhat)) == (yhat, Dict(:x => pihat))
    end
    return
end

function test_ValueFunction_Max()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Max,
        upper_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @constraint(sp, x.out == x.in)
        @stageobjective(sp, -2 * x.out)
    end
    SDDP.train(model; iteration_limit = 2, print_level = 0)
    V1 = SDDP.ValueFunction(model[1])
    for (xhat, yhat, pihat) in
        [(0.0, 0.0, 0.0), (1.0, 2.0, 2.0), (2.0, 4.0, 2.0)]
        (y, duals) = SDDP.evaluate(V1, Dict(:x => xhat))
        @test y == -yhat
        @test duals == Dict(:x => -pihat)
    end
    return
end

function test_ValueFunction_optimizer()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
        direct_mode = true,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @constraint(sp, x.out == x.in)
        @stageobjective(sp, 2 * x.out)
    end
    SDDP.train(model; iteration_limit = 2, print_level = 0)
    V1 = SDDP.ValueFunction(model[1])
    @test_throws JuMP.NoOptimizer() SDDP.evaluate(V1, Dict(:x => 1.0))
    JuMP.set_optimizer(V1, HiGHS.Optimizer)
    (y, _) = SDDP.evaluate(V1, Dict(:x => 1.0))
    @test y == 2.0
    return
end

function test_ValueFunction_objective_state()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        SDDP.add_objective_state(
            sp;
            initial_value = 0.0,
            lipschitz = 10.0,
        ) do p, ω
            return p + ω
        end
        @constraint(sp, x.out == x.in)
        SDDP.parameterize(sp, [1, 2]) do ω
            price = SDDP.objective_state(sp)
            @stageobjective(sp, price * x.out)
        end
    end
    SDDP.train(model; iteration_limit = 10, print_level = 0)
    V1 = SDDP.ValueFunction(model[1])
    @test_throws AssertionError SDDP.evaluate(V1, Dict(:x => 1.0))
    @test SDDP.evaluate(V1, Dict(:x => 1.0); objective_state = 1) ==
          (2.5, Dict(:x => 2.5))
    @test SDDP.evaluate(V1, Dict(:x => 0.0); objective_state = 2) ==
          (0.0, Dict(:x => 3.5))
    return
end

function test_ValueFunction_belief_state()
    graph = SDDP.MarkovianGraph(Matrix{Float64}[[0.5 0.5], [1.0 0.0; 0.0 1.0]])
    SDDP.add_ambiguity_set(graph, [(1, 1), (1, 2)])
    SDDP.add_ambiguity_set(graph, [(2, 1), (2, 2)])
    model = SDDP.PolicyGraph(
        graph;
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, node
        (t, i) = node
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @constraint(sp, x.out == x.in)
        P = [[0.2, 0.8], [0.8, 0.2]]
        SDDP.parameterize(sp, [1, 2], P[i]) do ω
            @stageobjective(sp, ω * x.out)
        end
    end
    SDDP.train(model; iteration_limit = 10, print_level = 0)
    V11 = SDDP.ValueFunction(model[(1, 1)])
    @test_throws AssertionError SDDP.evaluate(V11, Dict(:x => 1.0))
    b = Dict((1, 1) => 0.8, (1, 2) => 0.2)
    (y, duals) = SDDP.evaluate(V11, Dict(:x => 1.0); belief_state = b)
    @test duals[:x] ≈ y ≈ 1.68
    V12 = SDDP.ValueFunction(model[(1, 2)])
    (y, duals) = SDDP.evaluate(V12, Dict(:x => 1.0); belief_state = b)
    @test duals[:x] ≈ y ≈ 1.68
    return
end

function test_ValuaeFunction_plot()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @variable(sp, y >= 0, SDDP.State, initial_value = 0)
        @constraint(sp, x.out >= x.in)
        @constraint(sp, x.out >= 2 * x.in - 1)
        @constraint(sp, y.out == y.in)
        @stageobjective(sp, x.out + y.out)
    end
    SDDP.train(model; iteration_limit = 3, print_level = 0)
    V1 = SDDP.ValueFunction(model[1])
    SDDP.plot(V1; x = 0:0.1:2, y = 0, open = false)
    SDDP.plot(V1; x = 0:0.1:2, y = 0:0.1:2, open = false)
    return
end

function test_ValueFunction_risk_measure_min()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @constraint(sp, x.out == x.in)
        SDDP.parameterize(sp, [1, 2]) do w
            @stageobjective(sp, w * x.out)
            return
        end
        return
    end
    SDDP.train(
        model;
        iteration_limit = 2,
        print_level = 0,
        risk_measure = SDDP.CVaR(0.25),
        cut_type = SDDP.MULTI_CUT,
    )
    V1 = SDDP.ValueFunction(model[1])
    for (xhat, yhat, pihat) in
        [(0.0, 0.0, 0.0), (1.0, 2.0, 2.0), (2.0, 4.0, 2.0)]
        @test SDDP.evaluate(V1, Dict(:x => xhat)) == (yhat, Dict(:x => pihat))
    end
    return
end

function test_ValueFunction_risk_measure_max()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        upper_bound = 10.0,
        sense = :Max,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 1.5)
        @constraint(sp, x.out == x.in)
        SDDP.parameterize(sp, [1, 2]) do w
            @stageobjective(sp, w * x.out)
            return
        end
        return
    end
    SDDP.train(
        model;
        iteration_limit = 2,
        print_level = 0,
        risk_measure = SDDP.CVaR(0.25),
        cut_type = SDDP.MULTI_CUT,
    )
    V1 = SDDP.ValueFunction(model[1])
    for (xhat, yhat, pihat) in
        [(0.0, 0.0, 1.0), (1.0, 1.0, 1.0), (2.0, 2.0, 1.0)]
        @test SDDP.evaluate(V1, Dict(:x => xhat)) == (yhat, Dict(:x => pihat))
        @test SDDP.evaluate(V1, Dict("x" => xhat)) == (yhat, Dict(:x => pihat))
    end
    return
end

end  # module

TestValueFunctions.runtests()
