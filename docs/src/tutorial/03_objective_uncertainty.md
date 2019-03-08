```@meta
CurrentModule = Kokako
```

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
```jldoctest tutorial_two; filter=r"Confidence interval.+?\n"
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

----------------------------------------------------
         SDDP.jl (c) Oscar Dowson, 2017-19

Numerical stability report
  Non-zero Matrix range     [1e+00, 1e+00]
  Non-zero Objective range  [1e+00, 2e+02]
  Non-zero Bounds range     [2e+02, 2e+02]
  Non-zero RHS range        [2e+02, 2e+02]

 Iteration   Simulation      Bound        Time (s)
         1   3.37500e+04   8.17308e+03   5.00000e-02
         2   1.42308e+04   1.05060e+04   9.06000e-01
         3   1.12500e+04   1.06250e+04   9.07000e-01
         4   1.12500e+04   1.06250e+04   9.08000e-01
         5   9.37500e+03   1.06250e+04   9.09000e-01
         6   2.43750e+04   1.06250e+04   9.10000e-01
         7   3.37500e+04   1.06250e+04   9.11000e-01
         8   2.75000e+04   1.06250e+04   9.12000e-01
         9   1.87500e+03   1.06250e+04   9.40000e-01
        10   1.87500e+03   1.06250e+04   9.41000e-01

Terminating training with status: iteration_limit
----------------------------------------------------
Confidence interval: 11477.5 ± 780.7
Lower bound: 10625.0
```

This concludes our third tutorial for `SDDP.jl`. In the next tutorial,
[Basic IV: Markov uncertainty](@ref), we add stagewise-dependence to the
inflows using a Markov chain.
