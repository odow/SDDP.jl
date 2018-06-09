# Tutorial Ten: parallelism

The SDDP algorithm is highly parallelizable. In SDDP.jl, we chose to implement
an approach that minimizes inter-process communication. We call this
asynchronous SDDP.

In our implementation, one process is designated the *master* process and
the remainder are designated as *slaves*. Each slave receives a full copy of the
SDDP model and is set to work, performing iterations. At the end of each
iteration, the slave passes the master the cuts it discovered during the
iteration and receives any new cuts discovered by other slaves. The slave also
queries the master as to whether it should terminate, perform another iteration,
or perform a simulation. If the master requests a simulation (for example, to
calculate a confidence interval in order to test for convergence), the slave
returns the objective value of the simulation rather than a new set of cuts.

In this tutorial, we explain how to use the asynchronous solve feature of
SDDP.jl.

First, we need to add some extra processes to Julia. This can be done two ways.
We can start Julia using `julia -p N`, where `N` is the number of worker
processes to add, or we can use the `addprocs(N)` function while Julia is
running. The main process that co-ordinates everything is called the *master*
process, and the remote  processes are called *workers*. For example, to add two
workers, we run:
```julia
julia> addprocs(2)
2-element Array{Int64, 1}:
 2
 3
```
One way of running a command on all of the processes is to use the `@everywhere`
macro. We can also use the `myid()` function to query the index of the process.
For example:
```julia
julia> @everywhere println("Called from: $(myid())")
Called from: 1
        From worker 3:  Called from: 3
        From worker 2:  Called from: 2
```
Note that the order we receive things from the worker processes is not
deterministic.

Now what we need to do is to initialize the random number generator on each
process. However, we need to be careful not to use the same seed on each process
or each process will perform identical SDDP iterations!
```julia
julia> @everywhere srand(123 * myid())
```
This will set `srand(123)` on process `1`, `srand(246)` on process `2`, and so
on.

!!!note
    If you are using any data or function in the subproblem definition, these
    need to be copied to every process. The easiest way to do this is to place
    everything in a file and then run `@everywhere include("path/to/my/file")`.

Recall from [Tutorial Four: Markovian policy graphs](@ref) that our model is:
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0,
      markov_transition = Array{Float64, 2}[
          [ 1.0 ]',
          [ 0.75 0.25 ],
          [ 0.75 0.25 ; 0.25 0.75 ]
      ]
                                        ) do sp, t, i
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
        mupliplier * fuel_cost[t] * thermal_generation
    )
    if i == 1  # wet climate state
        setnoiseprobability!(sp, [1/6, 1/3, 0.5])
    else       # dry climate state
        setnoiseprobability!(sp, [0.5, 1/3, 1/6])
    end
end
```
We can solve this using SDDP.jl's asynchronous feature by passing an instance of
[`Asynchronous`](@ref) to the `solve_type` keyword in [`solve`](@ref):
```julia
status = solve(m,
    max_iterations = 10,
    solve_type     = Asynchronous()
)
```
If you have multiple processes, SDDP.jl will detect this and choose asynchronous
by default. You can force the serial solution by passing an instance of
[`Serial`](@ref) to `solve_type`.

The log is:
```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Asynchronous solver with 2 slave processors
    Model:
        Stages:         3
        States:         1
        Subproblems:    5
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------
       27.000K         5.818K        |     2    0.0      0    0.0    0.0
        0.000          6.198K        |     1    0.0      0    0.0    0.0
        9.000K         6.952K        |     3    0.0      0    0.0    0.0
        2.000K         7.135K        |     4    0.0      0    0.0    0.0
        2.000K         7.135K        |     5    0.0      0    0.0    0.0
        5.000K         7.135K        |     6    0.0      0    0.0    0.0
        5.000K         7.135K        |     7    0.0      0    0.0    0.0
        2.000K         7.135K        |     8    0.0      0    0.0    0.1
       24.000K         7.135K        |     9    0.0      0    0.0    0.1
       20.000K         7.135K        |    10    0.1      0    0.0    0.1
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         10
        Termination Status: max_iterations
===============================================================================
```
Note that the order of the *Cut #* column is not sequential because they are
numbered in order of when they were created.

That concludes our tenth tutorial for SDDP.jl.
