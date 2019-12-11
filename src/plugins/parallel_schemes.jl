#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    Serial()

Run SDDP in serial mode.
"""
struct Serial <: AbstractParallelScheme end
Base.show(io::IO, ::Serial) = print(io, "serial mode")

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
function Base.show(io::IO, a::Asynchronous)
    print(io, "Asynchronous mode with $(length(a.slave_ids)) procs.")
end

"""
    slave_update(model::PolicyGraph, result::IterationResult)

A callback called on a slave whenever a new result is available.
"""
function slave_update(model::PolicyGraph, result::IterationResult)
    for (node_index, cuts) in result.cuts
        for cut in cuts
            if cut === nothing
                error(
                    "This model uses features that are not suppored in async " *
                    "mode. Use `parallel_scheme = Serial()` instead.",
                )
            end
            _add_cut(
                model[node_index].bellman_function.global_theta,
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
    model = SDDP.__async_model__::PolicyGraph
    for (key, node) in model.nodes
        JuMP.set_optimizer(node.subproblem, node.optimizer)
    end
    return slave_loop(model, args...)
end

function slave_loop(
    model::PolicyGraph{T},
    jobs::Distributed.RemoteChannel{Channel{Options}},
    updates::Distributed.RemoteChannel{Channel{IterationResult{T}}},
    results::Distributed.RemoteChannel{Channel{IterationResult{T}}},
) where {T}
    try
        while true
            job = take!(jobs)
            while isready(updates)
                update = take!(updates)
                slave_update(model, update)
            end
            result = iteration(model, job)
            put!(results, result)
        end
    catch ex
        if isa(ex, Distributed.RemoteException) &&
           isa(ex.captured.ex, InvalidStateException)
            # The master process must have closed on us. Bail out without
            # consequence.
            return
        end
        @show dump(ex)
        rethrow(ex)
    end
end

function send_to(pid, key, val)
    Distributed.remotecall_fetch(Core.eval, pid, SDDP, Expr(:(=), key, val))
    return
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
        for pid in async.slave_ids
    )
    results = Distributed.RemoteChannel(() -> Channel{IterationResult{T}}(Inf))
    for pid in async.slave_ids
        send_to(pid, :__async_model__, model)
        put!(jobs, options)
        Distributed.remote_do(init_slave_loop, pid, jobs, updates[pid], results)
    end

    while true
        result = take!(results)
        # Add the result to the current model!
        slave_update(model, result)
        push!(
            options.log,
            Log(
                length(options.log) + 1,
                result.bound,
                result.cumulative_value,
                time() - options.start_time,
                result.pid,
            ),
        )
        options.dashboard_callback(options.log[end], false)
        if options.print_level > 0
            print_helper(print_iteration, options.log_file_handle, options.log[end])
        end
        if result.has_converged
            close(jobs)
            close(results)
            for (_, ch) in updates
                close(ch)
            end
            return result.status
        end
        for pid in async.slave_ids
            pid == result.pid && continue
            put!(updates[pid], result)
        end
        put!(jobs, options)
    end
end
