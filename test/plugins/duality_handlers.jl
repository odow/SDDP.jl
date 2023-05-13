#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestDualityHandlers

using SDDP
using Test
import HiGHS

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

# Single-stage model helps set up a node and subproblem to test dual
# calculations
function easy_single_stage(duality_handler)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = 0,
        optimizer = HiGHS.Optimizer,
    ) do sp, stage
        @variable(sp, x[1:2], Bin, SDDP.State, initial_value = 0)
        @variable(sp, y)
        if stage == 1
            @stageobjective(sp, 0)
        else
            @constraint(sp, y >= x[1].in + x[2].in)
            fix(x[1].in, 0)
            fix(x[2].in, 0)
            @stageobjective(sp, y)
        end
    end
    node = model.nodes[2]
    options =
        SDDP.Options(model, Dict(:x => 1.0); duality_handler = duality_handler)
    _ = SDDP.prepare_backward_pass(model, duality_handler, options)
    SDDP._initialize_solver(node; throw_error = false)
    optimize!(node.subproblem)
    obj, dual_vars = SDDP.get_dual_solution(node, duality_handler)

    if duality_handler == SDDP.ContinuousConicDuality()
        @test all(values(dual_vars) .<= ones(2))
    else
        @test all(values(dual_vars) .>= -ones(2))
    end
    return
end

# 'Exclusive or' function, no obvious choice of "tightest" dual
function xor_single_stage(duality_handler)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = 0,
        optimizer = HiGHS.Optimizer,
    ) do sp, stage
        @variable(sp, x[1:2], Bin, SDDP.State, initial_value = 1)
        @variable(sp, y)
        if stage == 1
            @stageobjective(sp, 0)
        else
            @constraints(sp, begin
                y >= x[1].in - x[2].in
                y >= x[2].in - x[1].in
                y <= x[1].in + x[2].in
                y <= 2 - x[1].in - x[2].in
            end)
            fix(x[1].in, 0)
            fix(x[2].in, 0)
            @stageobjective(sp, y)
        end
    end
    node = model.nodes[2]
    options =
        SDDP.Options(model, Dict(:x => 1.0); duality_handler = duality_handler)
    _ = SDDP.prepare_backward_pass(model, duality_handler, options)
    SDDP._initialize_solver(node; throw_error = false)
    optimize!(node.subproblem)
    obj, dual_vars = SDDP.get_dual_solution(node, duality_handler)
    if duality_handler == SDDP.ContinuousConicDuality()
        @test sum(values(dual_vars)) >= -1
    else
        @test sum(values(dual_vars)) <= 1
    end
    return
end

function test_easy_continuous()
    easy_single_stage(SDDP.ContinuousConicDuality())
    return
end

function test_easy_LagrangianDuality()
    easy_single_stage(SDDP.LagrangianDuality())
    return
end

function test_xor_continuous()
    xor_single_stage(SDDP.ContinuousConicDuality())
    return
end

function test_xor_LagrangianDuality()
    xor_single_stage(SDDP.LagrangianDuality())
    return
end

function test_prepare_backward_pass()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        direct_mode = false,
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
    options = SDDP.Options(
        model,
        Dict(:x => 1.0);
        duality_handler = SDDP.ContinuousConicDuality(),
    )
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
    undo_relax = SDDP.prepare_backward_pass(
        model,
        SDDP.ContinuousConicDuality(),
        options,
    )
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
    undo_relax()
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
    return
end

function test_kelleys_min()
    model = SDDP.LinearPolicyGraph(
        stages = 10,
        sense = :Min,
        lower_bound = -1000,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 1.1)
        @stageobjective(sp, (-5 + t) * x.out)
        @constraint(sp, x.out == x.in)
    end
    set_optimizer(model, HiGHS.Optimizer)
    for t in 1:10
        SDDP.parameterize(model[t], nothing)
        SDDP.set_incoming_state(model[t], Dict(:x => 1.1))
        JuMP.optimize!(model[t].subproblem)
        lobj, lagrange =
            SDDP.get_dual_solution(model[t], SDDP.LagrangianDuality())
        JuMP.optimize!(model[t].subproblem)
        cobj, conic =
            SDDP.get_dual_solution(model[t], SDDP.ContinuousConicDuality())
        @test isapprox(lobj, cobj, atol = 1e-5)
        csc, scd =
            SDDP.get_dual_solution(model[t], SDDP.StrengthenedConicDuality())
        @test csc == cobj
        for (k, v) in lagrange
            @test isapprox(v, conic[k], atol = 1e-5)
            @test conic[k] == scd[k]
        end
    end
    return
end

function test_kelleys_max()
    model = SDDP.LinearPolicyGraph(
        stages = 10,
        sense = :Max,
        upper_bound = 1000,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 1.1)
        @stageobjective(sp, (-5 + t) * x.out)
        @constraint(sp, x.out == x.in)
    end
    set_optimizer(model, HiGHS.Optimizer)
    for t in 1:10
        SDDP.parameterize(model[t], nothing)
        SDDP.set_incoming_state(model[t], Dict(:x => 1.1))
        JuMP.optimize!(model[t].subproblem)
        lobj, lagrange =
            SDDP.get_dual_solution(model[t], SDDP.LagrangianDuality())
        JuMP.optimize!(model[t].subproblem)
        cobj, conic =
            SDDP.get_dual_solution(model[t], SDDP.ContinuousConicDuality())
        @test isapprox(lobj, cobj, atol = 1e-5)
        csc, scd =
            SDDP.get_dual_solution(model[t], SDDP.StrengthenedConicDuality())
        @test csc == cobj
        for (k, v) in lagrange
            @test isapprox(v, conic[k], atol = 1e-5)
            @test conic[k] == scd[k]
        end
    end
    return
end

function test_kelleys_abs_function()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = -10.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 1.0)
        @constraint(sp, x.out >= 1.2(x.in - 1))
        @constraint(sp, x.out >= 0.1(x.in - 1))
        @constraint(sp, x.out >= -x.in)
        @stageobjective(sp, x.out)
    end
    set_optimizer(model, HiGHS.Optimizer)
    SDDP.parameterize(model[1], nothing)
    SDDP.set_incoming_state(model[1], Dict(:x => 0.5))
    JuMP.optimize!(model[1].subproblem)
    lobj, lagrange = SDDP.get_dual_solution(model[1], SDDP.LagrangianDuality())
    @test isapprox(lobj, -10.05, atol = 1e-5)
    @test isapprox(lagrange[:x], 0.1, atol = 1e-5)
    SDDP.set_incoming_state(model[1], Dict(:x => 1.5))
    JuMP.optimize!(model[1].subproblem)
    lobj, lagrange = SDDP.get_dual_solution(model[1], SDDP.LagrangianDuality())
    @test isapprox(lobj, -9.4, atol = 1e-5)
    @test isapprox(lagrange[:x], 1.2, atol = 1e-5)
    return
end

function test_kelleys_abs_function_max()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 10.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 1.0)
        @constraint(sp, x.out <= 1.2(x.in - 1))
        @constraint(sp, x.out <= 0.1(x.in - 1))
        @stageobjective(sp, x.out)
    end
    set_optimizer(model, HiGHS.Optimizer)
    SDDP.parameterize(model[1], nothing)
    SDDP.set_incoming_state(model[1], Dict(:x => 0.5))
    JuMP.optimize!(model[1].subproblem)
    lobj, lagrange = SDDP.get_dual_solution(model[1], SDDP.LagrangianDuality())
    @test isapprox(lobj, 9.4, atol = 1e-5)
    @test isapprox(lagrange[:x], 1.2, atol = 1e-5)
    SDDP.set_incoming_state(model[1], Dict(:x => 1.5))
    JuMP.optimize!(model[1].subproblem)
    lobj, lagrange = SDDP.get_dual_solution(model[1], SDDP.LagrangianDuality())
    @test isapprox(lobj, 10.05, atol = 1e-5)
    @test isapprox(lagrange[:x], 0.1, atol = 1e-5)
    return
end

"""
Test duality in a naturally integer problem
"""
function test_kelleys_ip_min()
    model = SDDP.LinearPolicyGraph(
        stages = 10,
        sense = :Min,
        lower_bound = -1000,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, Int, SDDP.State, initial_value = 1.0)
        @stageobjective(sp, (-5 + t) * x.out)
        @constraint(sp, x.out == x.in)
    end
    set_optimizer(model, HiGHS.Optimizer)
    for t in 1:10
        SDDP.parameterize(model[t], nothing)
        SDDP.set_incoming_state(model[t], Dict(:x => 1.0))
        JuMP.optimize!(model[t].subproblem)
        lobj, lagrange =
            SDDP.get_dual_solution(model[t], SDDP.LagrangianDuality())
        csc, scd =
            SDDP.get_dual_solution(model[t], SDDP.StrengthenedConicDuality())
        @test isapprox(lobj, csc, atol = 1e-5)
        for (k, v) in lagrange
            @test isapprox(v, scd[k], atol = 1e-5)
        end
    end
    return
end

function test_kelleys_ip_max()
    model = SDDP.LinearPolicyGraph(
        stages = 10,
        sense = :Max,
        upper_bound = 1000,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x, Int, SDDP.State, initial_value = 2.0)
        @stageobjective(sp, (-5 + t) * x.out)
        @constraint(sp, x.out == x.in)
    end
    set_optimizer(model, HiGHS.Optimizer)
    l = SDDP.LagrangianDuality()
    for t in 1:10
        SDDP.parameterize(model[t], nothing)
        SDDP.set_incoming_state(model[t], Dict(:x => 2.0))
        JuMP.optimize!(model[t].subproblem)
        lobj, lagrange = SDDP.get_dual_solution(model[t], l)
        csc, scd =
            SDDP.get_dual_solution(model[t], SDDP.StrengthenedConicDuality())
        @test isapprox(lobj, csc, atol = 1e-5)
        for (k, v) in lagrange
            @test isapprox(v, scd[k], atol = 1e-5)
        end
    end
    return
end

function test_LagrangianDuality_warn()
    @test_logs (:warn,) SDDP.LagrangianDuality(atol = 1e-6)
    return
end

function test_BanditDuality_show()
    @test sprint(show, SDDP.BanditDuality()) ==
          "BanditDuality with arms:\n * SDDP.ContinuousConicDuality()\n * SDDP.StrengthenedConicDuality()"
    return
end

function test_BanditDuality_show()
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = -100.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x[1:2] <= 5, SDDP.State, initial_value = 0.0)
        if t == 1
            @stageobjective(sp, -1.5 * x[1].out - 4 * x[2].out)
        else
            @variable(sp, 0 <= y[1:4] <= 1, Bin)
            @variable(sp, ω[1:2])
            @stageobjective(sp, -16 * y[1] - 19 * y[2] - 23 * y[3] - 28 * y[4])
            @constraint(
                sp,
                2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= ω[1] - x[1].in
            )
            @constraint(
                sp,
                6 * y[1] + 1 * y[2] + 3 * y[3] + 2 * y[4] <= ω[2] - x[2].in
            )
            steps = range(5, stop = 15, length = 10)
            SDDP.parameterize(sp, [[i, j] for i in steps for j in steps]) do φ
                return JuMP.fix.(ω, φ)
            end
        end
    end
    handler = SDDP.BanditDuality()
    SDDP.train(model, duality_handler = handler, iteration_limit = 20)
    @test sum(
        l.duality_key == " " for l in model.most_recent_training_results.log
    ) > 10
    @test sum(
        l.duality_key == "S" for l in model.most_recent_training_results.log
    ) > 0
    return
end

end

TestDualityHandlers.runtests()
