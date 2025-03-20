# Use multithreading

SDDP.jl can take advantage of the parallel nature of modern computers to solve
problems across multiple threads.

## Enabling threads

By default, starting Julia with `julia` enables only a single thread. To enable
more threads, you must start Julia from a command line using the `--threads`
flag, where `N` is the number of threads to use:
```julia
julia --threads N
```
Alteratively, you may set the `JULIA_NUM_THREADS` environment variable to `N`
**before** starting Julia. (There is no way to change the number of threads once
Julia has started.)

Once Julia has started, check how many threads you have available using
`Threads.nthreads()`. For this online documentation, we set the
`JULIA_NUM_THREADS` environment variable to `4` **before** starting Julia:

```jldoctest
julia> ENV["JULIA_NUM_THREADS"]
"4"

julia> Threads.nthreads()
4
```

## Multithreading with SDDP

To enable the multithreading algorithm in SDDP.jl, pass an instance of
[`SDDP.Threaded`](@ref) to the `parallel_scheme` argument of [`SDDP.train`](@ref)
and [`SDDP.simulate`](@ref).

```@repl
using SDDP, HiGHS
model = SDDP.LinearPolicyGraph(
    stages = 12,
    lower_bound = 0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(sp, x >= 0, SDDP.State, initial_value = 1)
    @stageobjective(sp, x.out)
end
SDDP.train(
    model;
    iteration_limit = 10,
    log_every_iteration = true,
    parallel_scheme = SDDP.Threaded(),
)
simulations = SDDP.simulate(
    model,
    100;
    parallel_scheme = SDDP.Threaded(),
    custom_recorders =
        Dict{Symbol,Function}(:thread_id => sp -> Threads.threadid()),
);
simulations[1][1][:thread_id]
simulations[26][1][:thread_id]
simulations[51][1][:thread_id]
simulations[76][1][:thread_id]
```

## Choosing the number of threads

As a rule of thumb, choose the number of threads to be less than or equal to the
number of physical cores you have on your machine or the number of nodes in the
graph, whichever is smaller.

The performance of the multithreaded algorithm in SDDP.jl is strongly limited by
the number of nodes in your policy graph; more nodes is better. The number of
threads that can be used is upper bounded by the number of nodes in the graph.
For example, if you have a graph with three nodes, SDDP.jl will never use more
than three threads, even if `--threads 4` is passed.

Even if you have a large number of nodes, using more threads than you have
physical cores is likely to lead to bad performance. For example, on a quad-core
laptop, using `--threads 52`, is likely to lead to bad performance, even if you
have a graph with 52 nodes.
