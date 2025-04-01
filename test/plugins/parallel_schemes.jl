#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Distributed

procs = Distributed.addprocs(4)

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
    # @__FILE__         := test/plugins/parallel_schemes.jl
    # @__DIR__          := test/plugins
    # dirname(@__DIR__) := test
    Pkg.activate(joinpath(dirname(@__DIR__), "Project.toml"))
end

@everywhere begin
    using Test
    using HiGHS
    using SDDP
end

function test_Asynchronous()
    a = SDDP.Asynchronous()
    @test a.slave_ids == Distributed.workers()
    @test length(a.slave_ids) > 1
    b = SDDP.Asynchronous([1, 2])
    @test b.slave_ids == [1, 2]
    return
end

function test_Asynchronous_optimizer()
    model = SDDP.LinearPolicyGraph(; stages = 2, lower_bound = 0.0) do sp, _
        @variable(sp, x, SDDP.State, initial_value = 0.0)
    end
    a = SDDP.Asynchronous(HiGHS.Optimizer)
    a.init_callback(model)
    @test solver_name(model[2].subproblem) == "HiGHS"
    model = SDDP.LinearPolicyGraph(;
        stages = 3,
        lower_bound = 0.0,
        direct_mode = true,
        optimizer = HiGHS.Optimizer,
    ) do sp, stage
        @variable(sp, 0 <= x <= 100, SDDP.State, initial_value = 0)
        @stageobjective(sp, x.in)
    end
    scheme = SDDP.Asynchronous()
    @test_throws(
        ErrorException(
            "Cannot use asynchronous solver with optimizers in direct mode.",
        ),
        scheme.init_callback(model),
    )
    @test_throws(
        ErrorException(
            "Cannot use asynchronous solver with optimizers in direct mode.",
        ),
        SDDP._uninitialize_solver(model; throw_error = true),
    )
    model = SDDP.LinearPolicyGraph(; stages = 3, lower_bound = 0.0) do sp, stage
        @variable(sp, 0 <= x <= 100, SDDP.State, initial_value = 0)
        @stageobjective(sp, x.in)
    end
    scheme = SDDP.Asynchronous()
    @test_throws(
        ErrorException(
            "You must supply an optimizer for the policy graph, either by passing\n" *
            "one to the `optimizer` keyword argument to `PolicyGraph`, or by\n" *
            "using `JuMP.set_optimizer(model, optimizer)`.\n",
        ),
        scheme.init_callback(model),
    )
    return
end

function test_slave_update()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Min,
        lower_bound = 0.0,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_upper_bound(x.out, ω)
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
        false,
    )
    SDDP.slave_update(model, result)
    cons = JuMP.all_constraints(
        model[1].subproblem,
        GenericAffExpr{Float64,VariableRef},
        MOI.GreaterThan{Float64},
    )
    @test length(cons) == 1
    obj = JuMP.constraint_object(cons[1])
    @test obj.set == MOI.GreaterThan(-5.0)
    @test length(obj.func.terms) == 2
    result = SDDP.IterationResult(
        1,
        0.0,
        0.0,
        false,
        :not_converged,
        Dict(
            1 => Any[
                (
                    theta = 1.0,
                    pi = Dict(:x => 2.0),
                    x = Dict(:x => 3.0),
                    obj_y = nothing,
                    belief_y = nothing,
                ),
                nothing,
            ],
            2 => Any[],
        ),
        false,
    )
    @test_throws ErrorException SDDP.slave_update(model, result)
    return
end

function test_async_solve()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0.0)
        @stageobjective(node, x.out)
        SDDP.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
            return JuMP.set_lower_bound(x.out, ω)
        end
    end
    solver =
        JuMP.optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
    SDDP.train(
        model;
        stopping_rules = [SDDP.IterationLimit(20)],
        parallel_scheme = SDDP.Asynchronous(solver; use_master = false),
    )
    @test SDDP.termination_status(model) == :iteration_limit
    @test all(l -> l.pid != 1, model.most_recent_training_results.log)
    @test SDDP.calculate_bound(model) == 6.0
    SDDP.train(
        model;
        stopping_rules = [SDDP.IterationLimit(20)],
        parallel_scheme = SDDP.Asynchronous(; use_master = true) do m
            for (key, node) in m.nodes
                JuMP.set_optimizer(node.subproblem, HiGHS.Optimizer)
                JuMP.set_silent(node.subproblem)
            end
        end,
    )
    @test SDDP.termination_status(model) == :iteration_limit
    @test any(l -> l.pid == 1, model.most_recent_training_results.log)
    @test SDDP.calculate_bound(model) == 6.0
    return
end

function test_simulate_parallel()
    model = SDDP.LinearPolicyGraph(;
        stages = 2,
        lower_bound = 0.0,
        sense = :Min,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, x[i = 1:2] >= i, SDDP.State, initial_value = 2i)
        @stageobjective(sp, x[1].out + x[2].out)
    end
    simulations = SDDP.simulate(
        model,
        20;
        custom_recorders = Dict{Symbol,Function}(
            :myid => (args...) -> Distributed.myid(),
        ),
        parallel_scheme = SDDP.Asynchronous(; use_master = false),
    )
    @test all([s[1][:myid] != 1 for s in simulations])
    return
end

function test_trap_error()
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
    return
end

test_Asynchronous()
test_slave_update()
test_async_solve()
test_simulate_parallel()
test_trap_error()

Distributed.rmprocs(procs)
