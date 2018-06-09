# Tutorial Six: distributionally robust SDDP

In the previous tutorial, we saw how risk measures can be used within our model.
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
Recall that our model for the hydrothermal scheduling problem  with risk
aversion from [Tutorial Five: risk](@ref) is:
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
      ],
           risk_measure = EAVaR(lambda=0.5, beta=0.1)
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

We will now replace the `EAVaR` with a description for a distributionally
robust set. All we need to provide is a *radius of uncertainty* for our
distributionally
