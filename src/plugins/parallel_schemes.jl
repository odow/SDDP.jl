#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    Serial()

Run SDDP in serial mode.
"""
struct Serial <: AbstractParallelScheme end

function master_loop(::Serial, model::PolicyGraph{T}, options::Options) where {T}
    while true
        result = iteration(model, options)
        if result.has_converged
            return result.status
        end
        options.dashboard_callback(options.log[end], false)
        if options.print_level > 0
            print_helper(print_iteration, options.log_file_handle, options.log[end])
        end
    end
end

"""
    Asynchronous(worker_ids::Vector{Int} = workers())

Run SDDP in asynchronous mode, using all available workers.
"""
struct Asynchronous <: AbstractParallelScheme
    slave_ids::Vector{Int}
    function Asynchronous(slave_ids::Vector{Int} = Distributed.workers())
        return new(slave_ids)
    end
end

"""
    slave_update(model::PolicyGraph, result::IterationResult)

A callback called on a slave whenever a new result is available.
"""
function slave_update(model::PolicyGraph, result::IterationResult)
    for (node, cuts) in result.cuts
        for cut in cuts
            if cut === nothing
                error(
                    "This model uses features that are not suppored in async " *
                    "mode. Use `parallel_scheme = Serial()` instead.",
                )
            end
            _add_cut(
                node.bellman_function.global_theta,
                cut.theta,
                cut.pi,
                cut.x,
                JuMP.AffExpr(0.0);
                cut_selection = false,
            )
        end
    end
    return
end

# Use the "function-barrier" technique to avoid type-instabilities in slave_loop.
function init_slave_loop(args...)
    return slave_loop(SDDP.__async_model__::PolicyGraph, args...)
end

function slave_loop(
    model::PolicyGraph{T},
    jobs::Distributed.RemoteChannel{Channel{Options}},
    updates::Distributed.RemoteChannel{Channel{IterationResult{T}}},
    results::Distributed.RemoteChannel{Channel{IterationResult{T}}},
) where {T}
    while isopen(jobs)
        job = take!(jobs)
        while isready(updates)
            update = take!(updates)
            slave_update(model, update)
        end
        result = iteration(model, job)
        put!(results, result)
    end
end

function master_loop(async::Asynchronous, model::PolicyGraph{T}, options::Options) where {T}
    # Initialize the remote channels. There are three types:
    # 1) jobs: master -> slaves: which stores a list of jobs that the slaves
    #       collectively pull from.
    # 2) updates: master -> slaves[i]: a unique channel for each slave, which
    #       is used to distribute results found by other slaves.
    # 3) results: slaves -> master: a channel which slaves collectively push to
    #       to feed the master new results.
    #
    # When initializing the slaves, we copy `model` to `SDDP.__async_model__`
    # on each slave. This is a global variable (bad), but we use the
    # function-barrier technique to improve performance when using model on each
    # slave.
    jobs = Distributed.RemoteChannel(() -> Channel{Options}(Inf))
    updates = Dict(
        pid => Distributed.RemoteChannel(() -> Channel{IterationResult{T}}(Inf))
        for pid in slave_ids
    )
    results = Distributed.RemoteChannel(() -> Channel{IterationResult{T}}(Inf))
    for pid in async.slave_ids
        Distributed.remotecall_fetch(
            Core.eval,
            pid,
            SDDP,
            Expr(:(=), :__async_model__, model),
        )
        put!(jobs, options)
        Distributed.remote_do(init_slave_loop, pid, jobs, updates[pid], results)
    end

    while isopen(jobs)
        wait(results)
        while isready(results)
            result = take!(results)
            options.dashboard_callback(options.log[end], false)
            if options.print_level > 0
                print_helper(print_iteration, options.log_file_handle, options.log[end])
            end
            if result.converged
                close(jobs)
                return result.status
            end
            for pid in async.slave_ids
                pid == result.pid && continue
                put!(updates[pid], result)
            end
            put!(jobs, options)
        end
    end
    return :error
end
