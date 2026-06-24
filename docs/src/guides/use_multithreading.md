# How to use multithreading

SDDP.jl can take advantage of the parallel nature of modern computers to solve
problems across multiple threads.

## Enabling threads

By default, starting Julia with `julia` enables only a single thread. To enable
more threads, you must start Julia from a command line using the `--threads`
flag, where `N` is the number of threads to use:
```julia
julia --threads N
```
Alternatively, you may set the `JULIA_NUM_THREADS` environment variable to `N`
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

## Threading with Gurobi

Gurobi's `Gurobi.Env` objects are not thread-safe. This means that you cannot
use the `optimizer = () -> Gurobi.Optimizer(env)` syntax if you are also using
multithreading. (You can use `optimizer = Gurobi.Optimizer` because this creates
a new environment for each node, but some license types have a limit on the
number of concurrent licenses that can be created.)

As a work-around, use the following function to partition the nodes into `N`
disjoint sets and use a shared lock across the nodes within a group so that at
most one subproblem within each group can be solved at any one time.
```julia
function shard_gurobi_licenses!(model::SDDP.PolicyGraph, groups...)
    for group in groups
        env, lock = Gurobi.Env(), ReentrantLock()
        for node in group
            set_optimizer(model[node].subproblem, () -> Gurobi.Optimizer(env))
            set_silent(model[node].subproblem)
            model[node].lock = lock
        end
    end
    return
end
```
For example, if you have two licenses you may want to group the nodes into "odd"
and "even":
```julia
julia> shard_gurobi_licenses!(model, 1:2:52, 2:2:52)
```
or, if you have three licenses, you may want to split the nodes into "start",
"middle", and "end" sets:
```julia
julia> shard_gurobi_licenses!(model, 1:17, 18:35, 36:52)
```

Here's a full example:
```julia
julia> using SDDP, Gurobi

julia> function shard_gurobi_licenses!(model::SDDP.PolicyGraph, groups...)
           for group in groups
               env, lock = Gurobi.Env(), ReentrantLock()
               for node in group
                   set_optimizer(model[node].subproblem, () -> Gurobi.Optimizer(env))
                   set_silent(model[node].subproblem)
                   model[node].lock = lock
               end
           end
           return
       end
shard_gurobi_licenses! (generic function with 1 method)

julia> model = SDDP.LinearPolicyGraph(;
           stages = 52,
           lower_bound = 0.0,
       ) do sp, t
           @variable(sp, 0 <= x <= 100, SDDP.State, initial_value = 0)
           @variable(sp, 0 <= u_production <= 200)
           @variable(sp, u_overtime >= 0)
           @constraint(sp, demand, x.in - x.out + u_production + u_overtime == 0)
           Ω = [100.0, 300.0]
           SDDP.parameterize(ω -> set_normalized_rhs(demand, ω), sp, Ω)
           @stageobjective(sp, 100 * u_production + 300 * u_overtime + 50 * x.out)
       end;

julia> shard_gurobi_licenses!(model, 1:2:52, 2:2:52)

julia> SDDP.train(model; parallel_scheme = SDDP.Threaded())
-------------------------------------------------------------------
         SDDP.jl (c) Oscar Dowson and contributors, 2017-26
-------------------------------------------------------------------
problem
  nodes           : 52
  state variables : 1
  scenarios       : 4.50360e+15
  existing cuts   : false
options
  solver          : Threaded()
  risk measure    : SDDP.Expectation()
  sampling scheme : SDDP.InSampleMonteCarlo
subproblem structure
  VariableRef                             : [5, 5]
  AffExpr in MOI.EqualTo{Float64}         : [1, 1]
  VariableRef in MOI.GreaterThan{Float64} : [4, 4]
  VariableRef in MOI.LessThan{Float64}    : [2, 3]
numerical stability report
  matrix range     [1e+00, 1e+00]
  objective range  [1e+00, 3e+02]
  bounds range     [1e+02, 2e+02]
  rhs range        [1e+02, 3e+02]
-------------------------------------------------------------------
 iteration    simulation      bound        time (s)     solves  pid
-------------------------------------------------------------------
         1   1.440000e+06  1.383984e+06  2.383709e-02       309   3
        22   1.607500e+06  1.432500e+06  1.177620e+00     13832   2
        62   1.405000e+06  1.432500e+06  2.469422e+00     30472   3
-------------------------------------------------------------------
status         : simulation_stopping
total time (s) : 2.469422e+00
total solves   : 30472
best bound     :  1.432500e+06
numeric issues : 0
-------------------------------------------------------------------
```

!!! tip
    If you need help with Gurobi licenses and multithreading, please open a
    GitHub issue.
