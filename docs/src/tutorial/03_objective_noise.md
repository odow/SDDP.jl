# Tutorial Three: objective noise

In [Tutorial One: first steps](@ref), we formulated a simple, deterministic
hydrothermal scheduling problem. Then, [Tutorial Two: RHS noise](@ref), we
extended this model to include stagewise-independent noise to the right-hand
side of a constraint. Now, in this tutorial, we extend the model to include
stagewise-independent noise in the objective function.

!!! note
    Notably, SDDP.jl does not allow stagewise-independent noise terms in the
    constraint matrix. However, this can be modelled using a Markovian policy
    graph like the one in [Tutorial Four: Markovian policy graphs](@ref).

Recall that our model for the hydrothermal scheduling problem  from
[Tutorial Two: RHS noise](@ref) is:
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

## Formulating the problem

In this tutorial, we are going to model the fuel cost of the thermal generator
by a stagewise-independent process. Specifically, we assume that in each stage,
there is an even probability of sampling a fuel cost of 80%, 100%, or 120% of
the usual fuel costs of \\\$50/MWh in the first stage, \\\$100/MWh in the second
stage, and \\\$150/MWh in the third stage. To add this noise term to the model,
we need to use a modified version of the [`@stageobjective`](@ref) macro
provided by SDDP.jl.

This version of [`@stageobjective`](@ref) is similar to the [`@rhsnoise`](@ref)
macro that we discussed in [Tutorial Two: RHS noise](@ref). It takes three
arguments. The first is the subproblem `sp`. The second argument is of the form
`name=[realizations]`, where `name` is a descriptive name, and `realizations` is
a vector of elements in the sample space. The third argument is any valid input
to the normal `@stageobjective` macro.

It is important to note that there must be the same number of realizations in
the objective as there are realizations in the right-hand-side random variable
(created using `@rhsnoise`). The two noise terms will be sampled in unison, so
that when the first element of the right-hand side noise is sampled, so to will
the first element of the objective noise. If the two noise terms should be
sampled independently, the user should form the Cartesian product.

For our example, the price multiplier and the inflows are negatively correlated.
Therefore, when `inflow=0.0`, the multiplier is `1.2`, when `inflow=50.0`, the
multiplier is `1.0`, and when `inflow=100.0`, the multiplier is `0.8`. Thus, we
have:
```julia
fuel_cost = [50.0, 100.0, 150.0]
@stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
    mupliplier * fuel_cost[t] * thermal_generation
)
```
As in [Tutorial Two: RHS noise](@ref), the noise terms are sampled using the
probability distribution set by the [`setnoiseprobability!`](@ref) function.

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
    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
        mupliplier * fuel_cost[t] * thermal_generation
    )
end
```

## Solving the problem

Now we need to solve the problem. As in the previous two tutorials, we use the
[`solve`](@ref) function. However, this time we use the bound stalling stopping
rule. This can be controlled via the `bound_stalling` keyword to
[`solve`](@ref). The syntax has a lot going on so we're going to give an example
of how it is used, and then walk through the different components.
```julia
status = solve(m,
    bound_stalling = BoundStalling(
        iterations = 5,
        rtol       = 0.0,
        atol       = 1e-6
    )
)
```

First, the `iterations` argument specifies how many iterations that the bound
must change by less than `atol` or `rtol` before terminating. For this example,
we choose to terminate the SDDP algorithm after the bound has failed to improve
for `5` iterations. Second, the `rtol` and `atol` keywords determine the
absolute and relative tolerances by which we compare the equality of the bound
in consecutive iterations. In this case, since the model is simple we choose an
absolute convergence tolerance of `1e-6`.

The termination status is `:bound_stalling`, and the output from the log is
now:
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
        7.500K         5.733K        |     1    0.0      0    0.0    0.0
       11.800K         8.889K        |     2    0.0      0    0.0    0.0
       14.000K         9.167K        |     3    0.0      0    0.0    0.0
       11.000K         9.167K        |     4    0.0      0    0.0    0.0
        9.000K         9.167K        |     5    0.0      0    0.0    0.0
        2.000K         9.167K        |     6    0.0      0    0.0    0.0
       14.000K         9.167K        |     7    0.0      0    0.0    0.0
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         7
        Termination Status: bound_stalling
===============================================================================
```

## Understanding the solution

Instead of performing a Monte Carlo simulation, you may want to simulate one
particular sequence of noise realizations. This *historical* simulation can
also be conducted using the [`simulate`](@ref) function.
```julia
simulation_result = simulate(m,
    [:outgoing_volume, :thermal_generation, :hydro_generation, :hydro_spill],
    noises = [1, 1, 3]
)
```
This time, `simulation_result` is a single dictionary. We can query the
objective of the simulation as follows:
```julia
julia> simulation_result[:objective]
9000.0
```
Interestingly, despite sampling the low-inflow, high-price realization in the
first stage, the model generates 150 MWh at a price of \\\$60/MWh:
```julia
julia> simulation_result[:thermal_generation]
3-element Array{Any, 1}:
 150.0
   0.0
   0.0
```

This concludes our third tutorial for SDDP.jl. In the next tutorial,
[Tutorial Four: Markovian policy graphs](@ref), we introduce stagewise-dependent
noise via a Markov chain.
