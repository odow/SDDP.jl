# Tutorial Two: RHS noise

In [Tutorial One: first steps](@ref), we formulated a simple hydrothermal
scheduling problem. In this tutorial, we extend the model to include
stagewise-independent noise in the right-hand side of the constraints.

!!! note
    Notably, SDDP.jl does not allow stagewise-independent noise terms in the
    constraint matrix. However, this can be modelled using a Markovian policy
    graph like the one in [Tutorial Four: Markovian policy graphs](@ref).

Recall that our model for the hydrothermal scheduling problem  from
[Tutorial One: first steps](@ref) is:
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0
                                        ) do sp, t
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
     end)
    inflow = [50.0, 50.0, 50.0]
    @constraints(sp, begin
        incoming_volume + inflow[t] - hydro_generation - hydro_spill == outgoing_volume
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

## Formulating the problem

In this tutorial, we are going to model inflows that are stagewise-independent.
Specifically, we assume that in each stage, there is an even probability of
sampling an inflow of `0.0`, `50.0`, or `100.0`. To add this noise term to the
model, we need to use the [`@rhsnoise`](@ref) macro provided by SDDP.jl.

[`@rhsnoise`](@ref) is similar to the JuMP `@constraint` macro. It takes three
arguments. The first is the subproblem `sp`. The second argument is of the form
`name=[realizations]`, where `name` is a descriptive name, and `realizations`
is a vector of elements in the sample space. The third argument is any valid
JuMP constraint that utilizes `name` in the right-hand side. For our example, we
have:
```julia
@rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
    outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
)
```
However, the realizations do not have to be the full right-hand side term. The
following is also valid:
```julia
inflows = [0.0, 50.0, 100.0]
@rhsnoise(sp, i = [1,2,3],
    outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflows[i]
)
```
We can set the probability of sampling each element in the sample space using
the [`setnoiseprobability!`](@ref) function. If `setnoiseprobability!` isn't
called, the distribution is assumed to be uniform. Despite this, for the sake of
completeness, we set the probability for our example as:
```julia
setnoiseprobability!(sp, [1/3, 1/3, 1/3])
```

Our model is now:
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0
                                        ) do sp, t
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

## Solving the problem

Now we need to solve the problem. As in [Tutorial One: first steps](@ref), we
use the [`solve`](@ref) function. However, this time we utilize some additional
arguments.

Since our problem is stochastic, we often want to simulate the policy in order
to estimate the upper (lower if maximizing) bound. This can be controlled via
the `simulation` keyword to [`solve`](@ref). The syntax has a lot going on so
we're going to give an example of how it is used, and then walk through the
different components.
```julia
status = solve(m,
    simulation = MonteCarloSimulation(
        frequency   = 2,
        confidence  = 0.95,
        termination = true
        min         = 50,
        step        = 50,
        max         = 100,
    )
)
```

First, the `frequency` argument specifies how often the  Monte Carlo simulation
is conducted (iterations/simulation). For this example, we conduct a Monte Carlo
simulation every two iterations. Second, the `confidence` specifies the level at
which to conduct the confidence interval. In this example, we construct a 95%
confidence interval. Third, the `termination` argument is a Boolean defining if
we should terminate the method if the lower limit of the confidence interval is
less than the lower bound (upper limit and bound for maximization problems). The
final three arguments implement the method of *sequential sampling*: `min` gives
the minimum number of replications to conduct before the construction of a
confidence interval. If there is evidence of convergence, another `step`
replications are conducted. This continues until either: (1) `max` number of
replications have been conducted; or (2) there is no evidence of convergence.
For our example, we conduct 50 replications, and if there is no evidence of
convergence, we conduct another 50 replications.

The output from the log is now:
```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Serial solver
    Model:
        Stages:         3
        States:         1
        Subproblems:    3
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------
       17.500K         3.438K        |     1    0.0      0    0.0    0.0
   7.606K   10.894K    7.500K   1.4  |     2    0.0     50    0.0    0.0
        7.500K         8.333K        |     3    0.0     50    0.0    0.0
   7.399K    9.651K    8.333K -11.2  |     4    0.0    150    0.1    0.1
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         4
        Termination Status: converged
===============================================================================
```

Compared with the log of a solve without using the `simulation` keyword, a few
things have changed. First, in the second and fourth rows (i.e. the iterations
in which a Monte Carlo simulation was conducted) the *Simulation* column now
gives two values. This pair is the confidence interval for the estimate of the
upper bound.

Second, in iterations in which a Monte Carlo simulation is
conducted, there is an entry in the *% Gap* column. This gaps measures the
difference between the lower limit of the simulated confidence interval and the
lower bound (given in the *Bound*) column. If the gap is positive, there is
evidence that the model has not converged. Once the gap is negative, the lower
bound lies above the lower limit of the confidence interval and we can terminate
the algorithm.

The third difference is that the *Simulations* column now records
the number of Monte Replications conducted to estimate the upper bound (in *#*)
and time performing those Monte Carlo replications (in *Time*). You can use this
information to tune the frequency at which the policy is tested for convergence.

Also observe that the first time we performed the Monte Carlo simulation, we
only conducted 50 replications; however, the second time we conducted 100. This
demonstrates the *sequential sampling* method at work.

Finally, the termination status is now `:converged` instead of
`:max_iterations`.

## Understanding the solution

The first thing we want to do is to query the lower (upper if maximizing) bound
of the solution. This can be done via the [`getbound`](@ref) function:
```julia
getbound(m)
```
This returns the value of the *Bound* column in the last row in the output table
above. In this example, the bound is `8333.0`.

Then, we can perform a Monte Carlo simulation of the policy using the
[`simulate`](@ref) function. We perform 500 replications and record the same
variables as we did in [Tutorial One: first steps](@ref).
```julia
simulation_result = simulate(m,
    500,
    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill]
)
```
This time, `length(simulation_result) = 500`. In addition to the variables, we
also record some additional fields. This includes `:stageobjective`, the value
of the stage-objective in each stage. We can calculate the cumulative objective
of each replication by summing the stage-objectives as follows:
```julia
julia> sum(simulation_result[100][:stageobjective])
2500.0
```
We can calculate the objective of all of each replication using Julia's
generator syntax:
```julia
julia> objectives = [sum(replication[:stageobjective]) for replication in simulation_result]
500-element Array{Float64, 1}:
  5000.0
 20000.0
 15000.0
 ⋮
```
Then, we can calculate the mean and standard deviation of these objectives:
```julia
julia> mean(objectives), std(objectives)
(8025.0, 5567.66)
```

We can query the noise that was sampled in each stage using the `:noise` key.
This returns the index of the noise from the vector of `realizations`. For
example:
```julia
julia> simulation_result[100][:noise]
3-element Array{Int, 1}:
 1
 3
 2
```

This concludes our second tutorial for SDDP.jl. In the next tutorial,
[Tutorial Three: objective noise](@ref), we introduce stagewise-independent
noise into the objective function.
