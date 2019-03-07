# Basic III: objective uncertainty

In the previous tutorial, [Basic II: adding uncertainty](@ref), we created a
stochastic hydro-thermal scheduling model. In this tutorial, we extend the
problem by adding uncertainty to the fuel costs.

Previously, we assumed that the fuel cost was deterministic: \\\$50/MWh in the
first stage, \\\$100/MWh in the second stage, and \\\$150/MWh in the third
stage. For this tutorial, we assume that in addition to these base costs, the
actual fuel cost is correlated with the inflows.

Our new model for the uncertinty is given by the following table:

| ω               |   1 |   2 |    3 |
| ----            | --- | --- | ---- |
| P(ω)            | 1/3 | 1/3 |  1/3 |
| inflow          |   0 |  50 |  100 |
| fuel_multiplier | 1.5 | 1.0 | 0.75 |


In stage `t`, the objective is not to minimize

`fuel_multiplier * fuel_cost[t] * thermal_generation`

## Creating a model

To add an uncertain objective, we can simply call [`@stageobjective`](@ref) from
inside the [`Kokako.parameterize`](@ref) function.

```jldoctest tutorial_two
using Kokako, GLPK

model = Kokako.LinearPolicyGraph(
            stages = 3,
            sense = :Min,
            lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    # Define the state variable.
    @variable(subproblem, 0 <= volume <= 200, Kokako.State, initial_value = 200)
    # Define the control variables.
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    # Define the constraints
    @constraints(subproblem, begin
        volume.out == volume.in + inflow - hydro_generation - hydro_spill
        thermal_generation + hydro_generation == 150.0
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    # Parameterize the subproblem.
    Ω = [
        (inflow = 0.0, fuel_multiplier = 1.5),
        (inflow = 50.0, fuel_multiplier = 1.0),
        (inflow = 100.0, fuel_multiplier = 0.75)
    ]
    Kokako.parameterize(subproblem, Ω, [1/3, 1/3, 1/3]) do ω
        JuMP.fix(inflow, ω.inflow)
        @stageobjective(subproblem,
            ω.fuel_multiplier * fuel_cost[t] * thermal_generation)
    end
end

# output

A policy graph with 3 nodes.
 Node indices: 1, 2, 3
```

## Training and simulating the policy

As in the previous two tutorials, we train the policy:
```jldoctest tutorial_two; filter=[r"\|.+?\n", r"Confidence interval.+?\n"]
Kokako.train(model; iteration_limit = 10)

simulations = Kokako.simulate(model, 500)

objective_values = [
    sum(stage[:stage_objective] for stage in sim) for sim in simulations
]

using Statistics

μ = round(mean(objective_values), digits = 2)
ci = round(1.96 * std(objective_values) / sqrt(500), digits = 2)

println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", round(Kokako.calculate_bound(model), digits = 2))

# output

———————————————————————————————————————————————————————————————————————————————
                        SDDP.jl - © Oscar Dowson, 2017-19.
———————————————————————————————————————————————————————————————————————————————
 Iteration | Simulation |      Bound |   Time (s)
———————————————————————————————————————————————————————————————————————————————
         1 |     7.500K |     8.173K |     0.046
         2 |    13.654K |    10.506K |     0.047
         3 |    31.607K |    10.625K |     0.049
         4 |    22.500K |    10.625K |     0.051
         5 |     1.875K |    10.625K |     0.053
         6 |     1.875K |    10.625K |     0.054
         7 |    24.375K |    10.625K |     0.055
         8 |    27.500K |    10.625K |     0.058
         9 |    11.250K |    10.625K |     0.060
        10 |    11.250K |    10.625K |     0.061
———————————————————————————————————————————————————————————————————————————————
 Terminating training with status: iteration_limit
———————————————————————————————————————————————————————————————————————————————
Confidence interval: 11342.5 ± 753.02
Lower bound: 10625.0
```

This concludes our third tutorial for `SDDP.jl`. In the next tutorial,
[Basic IV: Markov uncertainty](@ref), we add stagewise-dependence to the
inflows using a Markov chain.
