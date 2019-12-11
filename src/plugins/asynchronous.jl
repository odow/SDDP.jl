#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    Serial()

Run SDDP in serial mode.
"""
struct Serial <: AbstractParallelScheme end

function master_loop(
    ::Serial, model::PolicyGraph{T}, options::Options
) where {T}
    while true
        result = iteration(model, options)
        if result.has_converged
            return result.status
        end
        options.dashboard_callback(options.log[end], false)
        if options.print_level > 0
            print_helper(
                print_iteration, options.log_file_handle, options.log[end]
            )
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

struct AsyncResult{T}
    result::IterationResult
end

struct AsyncJob
    options::Options
end

"""
    slave_update(model::PolicyGraph{T}, result::AsyncResult{T}) where {T}

A callback called on a slave whenever a new result is available.
"""
function slave_update(model::PolicyGraph{T}, result::AsyncResult{T}) where {T}
    for (node, cuts) in result.result.cuts
        for cut in cuts
            if cut === nothing
                error(
                    "This model uses features that are not suppored in async " *
                    "mode. Use `parallel_scheme = Serial()` instead."
                )
            end
            _add_cut(
                node.bellman_function.global_theta,
                cut.theta,
                cut.pi,
                cut.x,
                JuMP.AffExpr(0.0);
                cut_selection = false
            )
        end
    end
    return
end

"""
    slave_work(model::PolicyGraph{T}, job::AsyncJob) where {T}

Process `job` on a slave.
"""
function slave_work(model::PolicyGraph{T}, job::AsyncJob) where {T}
    result = iteration(model, job.options)
    return AsyncResult(result)
end

# Use the "function-barrier" technique to avoid type-instabilities in slave_loop.
function init_slave_loop(args...)
    return slave_loop(SDDP.__async_model__::PolicyGraph, args...)
end

function slave_loop(
    model::PolicyGraph{T},
    jobs::Distributed.RemoteChannel{Channel{AsyncJob}},
    updates::Distributed.RemoteChannel{Channel{AsyncResult{T}}},
    results::Distributed.RemoteChannel{Channel{AsyncResult{T}}}
) where {T}
    while isopen(jobs)
        job = take!(jobs)
        while isready(updates)
            update = take!(updates)
            slave_update(model, update)
        end
        result = slave_work(model, job)
        put!(results, result)
    end
end

function async_master_loop(
    model::PolicyGraph{T};
    slave_ids::Vector{Int},
    new_job_callback::Function,
    new_result_callback::Function
) where {T}
    jobs = Distributed.RemoteChannel(() -> Channel{AsyncJob}(Inf))
    updates = Dict(
        pid => Distributed.RemoteChannel(() -> Channel{AsyncResult{T}}(Inf))
        for pid in slave_ids
    )
    results = Distributed.RemoteChannel(() -> Channel{AsyncResult{T}}(Inf))
    for pid in slave_ids
        Distributed.remotecall_fetch(
            Core.eval, pid, SDDP, Expr(:(=), :__async_model__, model)
        )
        new_job_callback(model, jobs)
        Distributed.remote_do(init_slave_loop, pid, jobs, updates[pid], results)
    end
    while isopen(jobs)
        wait(results)
        while isready(results)
            result = take!(results)
            new_result_callback(model, result)
            for pid in slave_ids
                pid == result.pid && continue
                put!(updates[pid], result)
            end
            new_job_callback(model, jobs)
        end
    end
end

function master_loop(
    async::Asynchronous, model::PolicyGraph{T}, options::Options
) where {T}
    has_converged = false
    status = :not_solved
    async_master_loop(
        model;
        slave_ids = async.slave_ids,
        new_job_callback = (model, jobs) -> begin
            if has_converged
                close(jobs)
            else
                put!(jobs, AsyncJob(options))
            end
        end,
        new_result_callback = (model, result) -> begin
            if result.result.converged
                has_converged = true
                status = result.result.status
            end
            options.dashboard_callback(options.log[end], false)
            if options.print_level > 0
                print_helper(
                    print_iteration, options.log_file_handle, options.log[end]
                )
            end
        end
    )
    return status
end
