#  Copyright 2017-19, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, JuMP, GLPK, Test

# Single-stage model helps set up a node and subproblem to test dual
# calculations
function easy_single_stage(integrality_handler)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = 0,
        integrality_handler = integrality_handler,
        optimizer = with_optimizer(GLPK.Optimizer)) do sp, stage
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
    integrality_handler == SDDP.ContinuousRelaxation() && SDDP.relax_integrality(model)
    JuMP.optimize!(node.subproblem)
    dual_vars = SDDP.get_dual_variables(node, integrality_handler)

    if integrality_handler == SDDP.ContinuousRelaxation()
        @test all(values(dual_vars) .<= ones(2))
    else
        @test all(values(dual_vars) .>= -ones(2))
        # Cannot add a variable without an upper bound
        @test_throws Exception @variable(node.subproblem, w, SDDP.State, initial_value = 0)
        @test_throws Exception @variable(node.subproblem, w, Int, SDDP.State, initial_value = 0)
    end
end

# 'Exclusive or' function, no obvious choice of "tightest" dual
function xor_single_stage(integrality_handler)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = 0,
        integrality_handler = integrality_handler,
        optimizer = with_optimizer(GLPK.Optimizer)) do sp, stage
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
    integrality_handler == SDDP.ContinuousRelaxation() && SDDP.relax_integrality(model)
    JuMP.optimize!(node.subproblem)
    dual_vars = SDDP.get_dual_variables(node, integrality_handler)

    if integrality_handler == SDDP.ContinuousRelaxation()
        @test sum(values(dual_vars)) >= -1
    else
        @test sum(values(dual_vars)) <= 1
    end
end

for integrality_handler in [SDDP.SDDiP(), SDDP.ContinuousRelaxation()]
    easy_single_stage(integrality_handler)
    xor_single_stage(integrality_handler)
end
