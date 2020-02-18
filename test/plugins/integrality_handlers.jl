#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using GLPK
using SDDP
using Test

@testset "get_dual_variables" begin
    # Single-stage model helps set up a node and subproblem to test dual
    # calculations
    function easy_single_stage(integrality_handler)
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            sense = :Min,
            lower_bound = 0,
            integrality_handler = integrality_handler,
            optimizer = GLPK.Optimizer,
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
        if integrality_handler == SDDP.ContinuousRelaxation()
            SDDP.relax_integrality(model, integrality_handler)
        end
        optimize!(node.subproblem)
        dual_vars = SDDP.get_dual_variables(node, integrality_handler)

        if integrality_handler == SDDP.ContinuousRelaxation()
            @test all(values(dual_vars) .<= ones(2))
        else
            @test all(values(dual_vars) .>= -ones(2))
            # Cannot add a variable without an upper bound
            @test_throws Exception @variable(
                node.subproblem,
                w,
                SDDP.State,
                initial_value = 0
            )
            @test_throws Exception @variable(
                node.subproblem,
                w,
                Int,
                SDDP.State,
                initial_value = 0
            )
        end
    end

    # 'Exclusive or' function, no obvious choice of "tightest" dual
    function xor_single_stage(integrality_handler)
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            sense = :Min,
            lower_bound = 0,
            integrality_handler = integrality_handler,
            optimizer = GLPK.Optimizer,
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
        if integrality_handler == SDDP.ContinuousRelaxation()
            SDDP.relax_integrality(model, integrality_handler)
        end
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
end

@testset "relax_integrality and enforce_integrality" begin
    model =
        SDDP.LinearPolicyGraph(stages = 2, lower_bound = 0.0, direct_mode = false) do sp, t
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
    @testset "ContinuousRelaxation" begin
        binaries, integers = SDDP.relax_integrality(model, SDDP.ContinuousRelaxation())
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
    @testset "SDDiP" begin
        binaries, integers = SDDP.relax_integrality(model, SDDP.SDDiP())
        @test isempty(binaries)
        @test isempty(integers)
    end
end

@testset "update_integrality_handler!" begin
    @testset "SDDiP" begin
        integrality_handler = SDDP.SDDiP()
        SDDP.update_integrality_handler!(integrality_handler, GLPK.Optimizer, 3)
        @test length(integrality_handler.subgradients) ==
        length(integrality_handler.old_rhs) ==
        length(integrality_handler.best_mult) ==
        length(integrality_handler.slacks) ==
        3
        @test integrality_handler.optimizer == GLPK.Optimizer
    end
    @testset "ContinuousRelaxation" begin
        integrality_handler = SDDP.ContinuousRelaxation()
        @test SDDP.update_integrality_handler!(integrality_handler, GLPK.Optimizer, 3) ==
              integrality_handler
    end
end

@testset "setup_state" begin
    function new_model(add_state)
        model = SDDP.PolicyGraph(
            SDDP.LinearGraph(2),
            lower_bound = 0.0,
            integrality_handler = SDDP.SDDiP(),
            direct_mode = false,
        ) do node, stage
            add_state(node)
        end
        return model
    end
    bin_state(node) = @variable(node, x, SDDP.State, Bin, initial_value = 0)
    model = new_model(bin_state)
    for stage = 1:2
        node = model[stage]
        @test haskey(node.states, :x)
        @test length(keys(node.states)) == 1
        @test node.states[:x] == node.subproblem[:x]
    end

    int_noupper(node) = @variable(node, x, SDDP.State, Int, initial_value = 0)
    @test_throws Exception new_model(int_noupper)
    int_state(node) = @variable(node, x <= 5, SDDP.State, Int, initial_value = 0)
    model = new_model(int_state)
    first_state = Symbol("_bin_x[1]")
    for stage = 1:2
        node = model[stage]
        @test haskey(node.states, first_state)
        @test length(keys(node.states)) == 3
    end
    cont_noupper(node) = @variable(node, x, SDDP.State, initial_value = 0)
    @test_throws Exception new_model(cont_noupper)
    cont_state(node) = @variable(node, x <= 5, SDDP.State, initial_value = 0)
    model = new_model(cont_state)
    for stage = 1:2
        node = model[stage]
        @test haskey(node.states, first_state)
        @test length(keys(node.states)) == 6
    end
    cont_state2(node) =
        @variable(node, x <= 5, SDDP.State, initial_value = 0, epsilon = 0.01)
    model = new_model(cont_state2)
    for stage = 1:2
        node = model[stage]
        @test haskey(node.states, first_state)
        @test length(keys(node.states)) == 9
    end
end
