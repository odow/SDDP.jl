#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestThreads

using SDDP
using Test
import HiGHS

function runtests()
    if Threads.nthreads() == 1
        return  # Skip tests if running in serial
    end
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_threaded()
    if haskey(ENV, "JULIA_NUM_THREADS")
        num_threads = get(ENV, "JULIA_NUM_THREADS", "0")
        @test parse(Int, num_threads) == Threads.nthreads()
        @test Threads.nthreads() > 1
    end
    c_eta, c_pt = [0.8, 0.5], [2, 5, 8, 11, 14]
    df_demand = rand(10:10:60, 24)
    model = SDDP.LinearPolicyGraph(;
        stages = 24,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do subproblem, t
        @variable(
            subproblem,
            0 <= x_volume[1:2] <= 8,
            SDDP.State,
            initial_value = 1,
        )
        @variable(subproblem, u_u[1:2] >= 0)
        @variable(subproblem, u_v[1:2] >= 0)
        @variable(subproblem, 0 <= u_x[1:5] <= 5, Int)
        @variable(subproblem, u_slack >= 0)
        @variable(subproblem, w_noise)
        @constraint(
            subproblem,
            [j in 1:2],
            x_volume[j].out == x_volume[j].in + c_eta[j] * u_v[j] - u_u[j],
        )
        @constraint(
            subproblem,
            sum(u_x) + sum(u_u) - sum(u_v) + u_slack == df_demand[t] + w_noise,
        )
        SDDP.parameterize(subproblem, [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]) do w
            return JuMP.fix(w_noise, w)
        end
        @stageobjective(subproblem, c_pt' * u_x + 35 * u_slack)
        return
    end
    SDDP.train(model; iteration_limit = 100, parallel_scheme = SDDP.Threaded())
    thread_ids_seen =
        Set{Int}(log.pid for log in model.most_recent_training_results.log)
    @test 2 <= length(thread_ids_seen) <= Threads.nthreads()
    recorder = Dict{Symbol,Function}(:thread_id => sp -> Threads.threadid())
    simulations = SDDP.simulate(
        model,
        100;
        parallel_scheme = SDDP.Threaded(),
        custom_recorders = recorder,
    )
    thread_ids_seen =
        Set{Int}(data[:thread_id] for sim in simulations for data in sim)
    @test length(thread_ids_seen) == Threads.nthreads()
    return
end

function test_threaded_warning()
    model = SDDP.PolicyGraph(
        SDDP.UnicyclicGraph(0.95);
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 1, SDDP.State, initial_value = 1)
        @stageobjective(sp, x.out)
        return
    end
    @test_logs((:warn,), SDDP.train(model; parallel_scheme = SDDP.Threaded()))
    return
end

end  # module

TestThreads.runtests()
