# Tutorial Six: distributionally robust SDDP

In [Tutorial Five: risk](@ref), we saw how risk measures can be used within our model.
In this tutorial we will learn how to incorporate a distributionally robust
optimization approach to SDDP in `SDDP.jl`.

Distributionally robust optimization (DRO) is a paradigm for optimizing under uncertainty.
In our setup, DRO is equivalent to a coherent risk measure, and we can apply DRO
by using the `risk_measure` keyword we saw previously.

## A little motivation on the concept
When we build a policy using SDDP, we impose a model on uncertain parameters.
When we come to use or evaluate our policy, our uncertainty may not resemble the
model we assumed.
For example, the hydrothermal scheduling model from the
previous tutorials assumed that inflows were independent between stages.
However, a real sequence of inflows is likely to exhibit correlation between stages.
Furthermore, we may wish to hold out a set of historical inflow sequences other than those
included while generating the policy, and evaluate the performance of the
policy based on its performance in the held out set. The held out set may also
correspond to inflows that are not stagewise independent.
A key motivation for a distributionally robust approach, is to avoid assuming an explicit model
on the probabilities of the scenarios we consider.
Instead, each time we come to add a cut, we assume that the probabilities
associated with each scenario are the worst case probabilities possible
(with respect to our objective), within some ambiguity set.

The implementation of distributionally robust SDDP here comes from the paper:
*A.B. Philpott, V.L. de Matos, L. Kapelevich* (2018): Distributionally Robust SDDP,
Computational Management Science,
[(link)](http://link.springer.com/article/10.1007/s10287-018-0314-0 "")
where the details of the approach are described.

## Formulating the problem
Recall that our model for the hydrothermal scheduling problem with RHS
uncertainty from [Tutorial Two: RHS noise](@ref) is:
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = GurobiSolver(OutputFlag=0),
        objective_bound = 0.0
                                        ) do sp, t
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    if t == 1
        @rhsnoise(sp, inflow = [50.0, 50.0, 50.0],
            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
      )
    else
        @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
      )
    end
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```
This model assumes that the probability of each inflow is equally likely, with
probability 0.33. The lower bound we converged to was 8.33k.

To describe a distributionally robust version of this problem, we need to choose a
*radius of uncertainty*. This is the maximum distance around the default probability vector ([0.33, 0.33, 0.33]) we will consider in our ambiguity set.
For a problem with S noises, this radius should be less than sqrt((S-1)/S)
`\\sqrt{\\frac{S-1}{S}}`.
(which would be the same as `TheWorstCase()` from [Tutorial Five: risk](@ref)).

!!! note
    We currently assume our uncertainty set is a ball centered around the
    probability vector assigning equal probabilities to all noises. We could
    generalize the algorithm to alter the center of the ball, but this is not
    currently implemented (see [this issue](https://github.com/odow/SDDP.jl/issues/117)).

Suppose, for example, we choose the radius of uncertainty to be `1/6`. We can
implement this by inserting a `DRO(1/6)` for the keyword argument `risk_measure`
we saw earlier.

This gives the new model
```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0,
           risk_measure = DRO(1/6)
                                        ) do sp, t
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    if t == 1
        @rhsnoise(sp, inflow = [50.0, 50.0, 50.0],
            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
        )
    else
        @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
        )
    end
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```

## Solving the problem
We can solve the above problem terminating with a maximum number of iterations
criterion. For example,
```julia
solve(m, iteration_limit = 10)
```
gives the following output log:
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
       17.500K        10.023K        |     1    0.0      0    0.0    0.0
        5.000K        10.023K        |     2    0.0      0    0.0    0.0
       17.500K        10.023K        |     3    0.0      0    0.0    0.0
       12.500K        10.023K        |     4    0.0      0    0.0    0.0
        5.000K        10.023K        |     5    0.0      0    0.0    0.0
        5.000K        10.023K        |     6    0.0      0    0.0    0.0
        5.000K        10.023K        |     7    0.0      0    0.0    0.0
        5.000K        10.023K        |     8    0.0      0    0.0    0.0
        5.000K        10.023K        |     9    0.0      0    0.0    0.0
       10.000K        10.023K        |    10    0.0      0    0.0    0.0
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         10
        Termination Status: iteration_limit
===============================================================================
:iteration_limit
```
We have converged to a lower bound of roughly 10.023k. This is a little lower
than the worst case of $15k.

To run a sanity check, let us set the radius to be sufficiently large in order
to match `TheWorstCase()`:

```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0,
           risk_measure = DRO(sqrt(2/3)-1e-4)
                                        ) do sp, t
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    if t == 1
        @rhsnoise(sp, inflow = [50.0, 50.0, 50.0],
            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
        )
    else
        @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
            outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
        )
    end
    setnoiseprobability!(sp, [1/3, 1/3, 1/3])
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end
```
