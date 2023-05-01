#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function _should_log(options, log_frequency::Int)
    return options.print_level > 0 &&
           mod(length(options.log), log_frequency) == 0
end

function _should_log(options, log_frequency::Function)
    return options.print_level > 0 && log_frequency(options.log)
end

function log_iteration(options)
    options.dashboard_callback(options.log[end], false)
    if _should_log(options, options.log_frequency)
        print_helper(print_iteration, options.log_file_handle, options.log[end])
        flush(options.log_file_handle)
    end
    return
end

"""
    Serial()

Run SDDP in serial mode.
"""
struct Serial <: AbstractParallelScheme end

Base.show(io::IO, ::Serial) = print(io, "serial mode")

interrupt(::Serial) = nothing

function master_loop(
    ::Serial,
    model::PolicyGraph{T},
    options::Options,
) where {T}
    _initialize_solver(model; throw_error = false)
    while true
        result = iteration(model, options)
        log_iteration(options)
        if result.has_converged
            return result.status
        end
    end
    return
end

function _simulate(
    model::PolicyGraph,
    ::Serial,
    number_replications::Int,
    variables::Vector{Symbol};
    kwargs...,
)
    _initialize_solver(model; throw_error = false)
    return map(
        i -> _simulate(model, variables; kwargs...),
        1:number_replications,
    )
end

struct Asynchronous <: AbstractParallelScheme
    init_callback::Function
    slave_ids::Vector{Int}
    use_master::Bool
end

"""
    Asynchronous(
        [init_callback::Function,]
        slave_pids::Vector{Int} = workers();
        use_master::Bool = true,
    )

Run SDDP in asynchronous mode workers with pid's `slave_pids`.

After initializing the models on each worker, call `init_callback(model)`. Note
that `init_callback` is run _locally on the worker_ and _not_ on the master
thread.

If `use_master` is `true`, iterations are also conducted on the master process.
"""
function Asynchronous(
    init_callback::Function,
    slave_ids::Vector{Int} = Distributed.workers();
    use_master::Bool = true,
)
    return Asynchronous(init_callback, slave_ids, use_master)
end

function Asynchronous(
    slave_ids::Vector{Int} = Distributed.workers();
    use_master::Bool = true,
)
    return Asynchronous(slave_ids; use_master = use_master) do model
        return _initialize_solver(model; throw_error = true)
    end
end

"""
    Asynchronous(
        solver::Any,
        slave_pids::Vector{Int} = workers();
        use_master::Bool = true,
    )

Run SDDP in asynchronous mode workers with pid's `slave_pids`.

Set the optimizer on each worker by calling `JuMP.set_optimizer(model, solver)`.
"""
function Asynchronous(
    solver::Any,
    slave_ids::Vector{Int} = Distributed.workers();
    use_master::Bool = true,
)
    return Asynchronous(slave_ids; use_master = use_master) do model
        return JuMP.set_optimizer(model, solver)
    end
end

interrupt(a::Asynchronous) = Distributed.interrupt(a.slave_ids)

function Base.show(io::IO, a::Asynchronous)
    return print(io, "Asynchronous mode with $(length(a.slave_ids)) workers.")
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
                cut.obj_y,
                cut.belief_y;
                cut_selection = true,
            )
        end
    end
    return
end

function slave_loop(
    async::Asynchronous,
    model::PolicyGraph{T},
    options::Options,
    updates::Distributed.RemoteChannel{Channel{IterationResult{T}}},
    results::Distributed.RemoteChannel{Channel{IterationResult{T}}},
) where {T}
    try
        async.init_callback(model)
        results_to_add = IterationResult{T}[]
        while true
            result = iteration(model, options)
            # The next four lines are subject to a race condition: if the master closes
            # `results` _after_ the call to `isopen` and _before_` the call to `put!` has
            # executed, we get an `InvalidStateException`. This gets trapped in the outer
            # try-catch.
            if !isopen(results)
                break
            end
            put!(results, result)
            # Instead of pulling a result from `updates` and adding it immediately, we want
            # to pull as many as possible in a short amount of time, the add them all and
            # start the loop again. Otherwise, by the time we've finished updating the
            # slave, there might be a new update :(
            while isready(updates)
                push!(results_to_add, take!(updates))
            end
            for result in results_to_add
                slave_update(model, result)
            end
            empty!(results_to_add)
        end
    catch ex
        trap_error(ex)
    end
    return
end

trap_error(ex::Exception) = throw(ex)
trap_error(::InterruptException) = nothing
trap_error(::InvalidStateException) = nothing
trap_error(ex::CapturedException) = trap_error(ex.ex)
trap_error(ex::Distributed.RemoteException) = trap_error(ex.captured)

function master_loop(
    async::Asynchronous,
    model::PolicyGraph{T},
    options::Options,
) where {T}
    # Initialize the remote channels. There are two types:
    # 1) updates: master -> slaves[i]: a unique channel for each slave, which
    #       is used to distribute results found by other slaves.
    # 2) results: slaves -> master: a channel which slaves collectively push to
    #       to feed the master new results.
    updates = Dict(
        pid => Distributed.RemoteChannel(
            () -> Channel{IterationResult{T}}(Inf),
        ) for pid in async.slave_ids
    )
    results = Distributed.RemoteChannel(() -> Channel{IterationResult{T}}(Inf))
    futures = Distributed.Future[]
    _uninitialize_solver(model; throw_error = true)
    for pid in async.slave_ids
        let model_pid = model, options_pid = options
            f = Distributed.remotecall(
                slave_loop,
                pid,
                async,
                model_pid,
                options_pid,
                updates[pid],
                results,
            )
            push!(futures, f)
        end
    end
    _initialize_solver(model; throw_error = true)
    while true
        # Starting workers has a high overhead. We have to copy the models across, and then
        # precompile all the methods on every process :(. While that's happening, let's
        # start running iterations on master. It has the added benefit that if the master
        # is ever idle waiting for a result from a slave, it will do some useful work :).
        #
        # It also means that Asynchronous() can be our default setting, since if there are
        # no workers, ther should be no overhead, _and_ this inner loop is just the serial
        # implementation anyway.
        while async.use_master && !isready(results)
            result = iteration(model, options)
            for (_, ch) in updates
                put!(ch, result)
            end
            log_iteration(options)
            if result.has_converged
                close(results)
                wait.(futures)
                return result.status
            end
        end
        while !isready(results)
            sleep(1.0)
        end
        # We'll only reach here is isready(results) == true, so we won't hang waiting for a
        # new result on take!. After we receive a new result from a slave, there are a few
        # things to do:
        # 1) send the result to the other slaves
        # 2) update the master problem with the new cuts
        # 3) compute the revised bound, update the log, and print to screen
        # 4) test for convergence (e.g., bound stalling, time limit, iteration limit)
        # 5) Exit, killing the running task on the workers.
        result = take!(results)
        for pid in async.slave_ids
            if pid != result.pid
                put!(updates[pid], result)
            end
        end
        slave_update(model, result)
        bound = calculate_bound(model)
        push!(
            options.log,
            Log(
                length(options.log) + 1,
                bound,
                result.cumulative_value,
                time() - options.start_time,
                result.pid,
                model.ext[:total_solves],
                duality_log_key(options.duality_handler),
                result.numerical_issue,
            ),
        )
        log_iteration(options)
        has_converged, status =
            convergence_test(model, options.log, options.stopping_rules)
        if has_converged
            close(results)
            wait.(futures)
            return status
        end
    end
    return
end

function _simulate(
    model::PolicyGraph,
    async::Asynchronous,
    number_replications::Int,
    variables::Vector{Symbol};
    kwargs...,
)
    _uninitialize_solver(model; throw_error = true)
    wp = Distributed.CachingPool(async.slave_ids)
    let model = model,
        init = false,
        async = async,
        variables = variables,
        kwargs = kwargs

        return Distributed.pmap(wp, 1:number_replications) do _
            if !init
                async.init_callback(model)
                init = true
            end
            return _simulate(model, variables; kwargs...)
        end
    end
    return
end
