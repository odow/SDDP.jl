#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Distributed

# !!! IMPORTANT !!!
#
# Workers **DON'T** inherit their parent's Pkg environment!
# Here's the relevant Julia issue: https://github.com/JuliaLang/julia/issues/28781
#
# This can cause reeeeeaally hard to track down bugs because
# a) workers may have different versions of packages on them
# b) you will run into a lot of complilation errors depending on the order of
#    code loading.
#
# As hack, run the following script:
@everywhere begin
    import Pkg
    Pkg.activate(".")
end

@everywhere begin
    using Test
    using GLPK
    using SDDP
end

@testset "Asynchronous" begin
    a = SDDP.Asynchronous()
    @test a.slave_ids == Distributed.workers()
    b = SDDP.Asynchronous([1, 2])
    @test b.slave_ids == [1, 2]
end

@testset "slave_update" begin
    model =
        SDDP.LinearPolicyGraph(stages = 2, sense = :Min, lower_bound = 0.0) do node, stage
            @variable(node, x, SDDP.State, initial_value = 0.0)
            @stageobjective(node, x.out)
            SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x.out, ω)
            end
        end

    result = SDDP.IterationResult(
        1,
        0.0,
        0.0,
        false,
        :not_converged,
        Dict(
            1 => Any[(
                theta = 1.0,
                pi = Dict(:x => 2.0),
                x = Dict(:x => 3.0),
                obj_y = nothing,
                belief_y = nothing,
            )],
            2 => Any[],
        ),
    )
    SDDP.slave_update(model, result)
    cons = JuMP.all_constraints(
        model[1].subproblem,
        GenericAffExpr{Float64,VariableRef},
        MOI.GreaterThan{Float64},
    )
    @test length(cons) == 1
    @test replace(sprint(print, cons[1]), "≥" => ">=") == "noname - 2 x_out >= -5.0"

    result = SDDP.IterationResult(
        1,
        0.0,
        0.0,
        false,
        :not_converged,
        Dict(
            1 => Any[(
                theta = 1.0,
                pi = Dict(:x => 2.0),
                x = Dict(:x => 3.0),
                obj_y = nothing,
                belief_y = nothing,
            ), nothing],
            2 => Any[],
        ),
    )
    @test_throws ErrorException SDDP.slave_update(model, result)
end

@testset "Async solve" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            JuMP.set_lower_bound(x.out, ω)
        end
    end
    SDDP.train(
        model,
        iteration_limit = 20,
        parallel_scheme = SDDP.Asynchronous(; use_master = false) do m
            for (key, node) in m.nodes
                JuMP.set_optimizer(node.subproblem, GLPK.Optimizer)
            end
        end,
    )
    @test SDDP.termination_status(model) == :iteration_limit
    @test all(l -> l.pid != 1, model.most_recent_training_results.log)
    @test SDDP.calculate_bound(model) == 6.0
    SDDP.train(
        model,
        iteration_limit = 20,
        parallel_scheme = SDDP.Asynchronous(; use_master = true) do m
            for (key, node) in m.nodes
                JuMP.set_optimizer(node.subproblem, GLPK.Optimizer)
            end
        end,
    )
    @test SDDP.termination_status(model) == :iteration_limit
    @test any(l -> l.pid == 1, model.most_recent_training_results.log)
    @test SDDP.calculate_bound(model) == 6.0
end

@testset "simulate parallel" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0.0,
        sense = :Min,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x[i = 1:2] >= i, SDDP.State, initial_value = 2i)
        @stageobjective(sp, x[1].out + x[2].out)
    end
    simulations = SDDP.simulate(
        model,
        20;
        custom_recorders = Dict{Symbol,Function}(:myid => (args...) -> Distributed.myid()),
        parallel_scheme = SDDP.Asynchronous(use_master = false),
    )
    @test all([s[1][:myid] != 1 for s in simulations])
end

@testset "trap_error" begin
    @test SDDP.trap_error(InvalidStateException("a", :a)) === nothing
    @test SDDP.trap_error(InterruptException()) === nothing
    ex = DomainError(-1.0)
    @test_throws ex SDDP.trap_error(ex)
    flag = true
    ex = try
        throw(InterruptException())
        flag = false
    catch ex
        Distributed.RemoteException(CapturedException(ex, catch_backtrace()))
    end
    @test SDDP.trap_error(ex) === nothing
    @test flag == true
end
