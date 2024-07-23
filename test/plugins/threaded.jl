#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function test_threaded()
    num_periods = 24
    num_thermal = 5
    num_battery = 2
    c_capacity = 14
    c_VOLL = 35
    c_eta = [0.8, 0.5]
    c_pt = [2, 5, 8, 11, 14]
    c_q = fill(c_capacity, 5)
    df_demand = rand(10:10:60, num_periods)
    model = SDDP.LinearPolicyGraph(;
        stages = num_periods,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do subproblem, t
        @variables(subproblem, begin
            0 <= x_volume[1:num_battery] <= 8, SDDP.State, (initial_value = 1)
            u_u[1:num_battery] >= 0
            u_v[1:num_battery] >= 0
            0 <= u_x[i in 1:num_thermal] <= c_q[i]
            u_slack >= 0
            w_noise
        end)
        @constraints(subproblem, begin
            battery[j in 1:num_battery],
                x_volume[j].out == x_volume[j].in + c_eta[j] * u_v[j] - u_u[j]
            sum(u_x) + sum(u_u) - sum(u_v) + u_slack == df_demand[t] + w_noise
        end)
        O = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
        SDDP.parameterize(subproblem, O) do omega
            JuMP.fix(w_noise, omega)
            return
        end
        @stageobjective(subproblem, c_pt' * u_x + u_slack * c_VOLL)
        return
    end
    SDDP.train(
        model;
        iteration_limit = 100,
        parallel_scheme = SDDP.Threaded(),
    )
    thread_ids_seen =
        Set{Int}(log.pid for log in model.most_recent_training_results.log)
    if Threads.nthreads() == 1
        @test length(thread_ids_seen) == 1
    else
        @test length(thread_ids_seen) > 1
    end
    simulations = SDDP.simulate(
        model,
        100;
        parallel_scheme = SDDP.Threaded(),
        custom_recorders = Dict{Symbol, Function}(
            :thread_id => sp -> Threads.threadid(),
        ),
    )
    thread_ids_seen = Set{Int}()
    for sim in simulations
        thread_ids = unique(data[:thread_id] for data in sim)
        @test length(thread_ids) == 1
        push!(thread_ids_seen, only(thread_ids))
    end
    if Threads.nthreads() == 1
        @test length(thread_ids_seen) == 1
    else
        @test length(thread_ids_seen) > 1
    end
    return
end

test_threaded()
