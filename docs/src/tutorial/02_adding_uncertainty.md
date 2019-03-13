```@meta
CurrentModule = Kokako
```

# Basic II: adding uncertainty

In the previous tutorial, [Basic I: first steps](@ref), we created a
deterministic  hydro-thermal scheduling model. In this tutorial, we extend the
problem by adding uncertainty.

Notably missing from our previous model were inflows. Inflows are the water that
flows into the reservoir through rainfall or rivers. These inflows are
uncertain, and are the cause of the main trade-off in hydro-thermal scheduling:
the desire to use water now to generate cheap electricity, against the risk that
future inflows will be low, leading to blackouts or expensive thermal
generation.

For our simple model, we assume that the inflows can be modelled by a discrete
distribution with the three outcomes given in the following table:

| ω    |   0 |  50 | 100 |
| ---- | --- | --- | --- |
| P(ω) | 1/3 | 1/3 | 1/3 |

The value of the noise (the random variable) is observed by the agent at the
start of each stage. This makes the problem a _wait-and-see_ or
_hazard-decision_ formulation.

To represent this, we can draw the following picture. The wavy lines denote the
uncertainty arriving into the start of each stage (node).

![Linear policy graph](../assets/stochastic_linear_policy_graph.png)

In addition to adding this uncertainty to the model, we also need to modify the
dynamics to include `inflow`:

`volume.out = volume.in + inflow - hydro_generation - hydro_spill`


## Creating a model

To add an uncertain variable to the model, we create a new JuMP variable
`inflow`, and then call the function [`Kokako.parameterize`](@ref). The
[`Kokako.parameterize`](@ref) function takes three arguments: the subproblem,
a vector of realizations, and a corresponding vector of probabilities.

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
        demand_constraint, thermal_generation + hydro_generation == 150.0
    end)
    # Define the objective for each stage `t`. Note that we can use `t` as an
    # index for t = 1, 2, 3.
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(subproblem, fuel_cost[t] * thermal_generation)
    # Parameterize the subproblem.
    Kokako.parameterize(subproblem, [0.0, 50.0, 100.0], [1/3, 1/3, 1/3]) do ω
        JuMP.fix(inflow, ω)
    end
end

# output

A policy graph with 3 nodes.
 Node indices: 1, 2, 3
```

Note how we use the JuMP function [`JuMP.fix`](http://www.juliaopt.org/JuMP.jl/v0.19/variables/#JuMP.fix)
to set the value of the `inflow` variable to `ω`.

!!! note
    [`Kokako.parameterize`](@ref) can only be called once in each subproblem
    definition!

## Training and simulating the policy

As in [Basic I: first steps](@ref), we train the policy:

```jldoctest tutorial_two
julia> Kokako.train(model; iteration_limit = 10)
-------------------------------------------------------
         SDDP.jl (c) Oscar Dowson, 2017-19

Numerical stability report
  Non-zero Matrix range     [1e+00, 1e+00]
  Non-zero Objective range  [1e+00, 2e+02]
  Non-zero Bounds range     [2e+02, 2e+02]
  Non-zero RHS range        [2e+02, 2e+02]
No problems detected

 Iteration    Simulation       Bound         Time (s)
        1    1.750000e+04   3.437500e+03   5.721000e+00
        2    1.093750e+04   7.500000e+03   6.552000e+00
        3    1.000000e+04   8.333333e+03   6.553000e+00
        4    1.250000e+04   8.333333e+03   6.554000e+00
        5    2.500000e+03   8.333333e+03   6.554000e+00
        6    5.000000e+03   8.333333e+03   6.555000e+00
        7    1.500000e+04   8.333333e+03   6.556000e+00
        8    1.250000e+04   8.333333e+03   6.557000e+00
        9    5.000000e+03   8.333333e+03   6.558000e+00
       10    1.250000e+04   8.333333e+03   6.558000e+00

Terminating training with status: iteration_limit
-------------------------------------------------------
```

!!! note
    Since SDDP is a stochastic algorithm, you might get slightly different
    numerical results.

We can also simulate the policy. Note that this time, the simulation is
stochastic. One common approach to quantify the quality of the policy is to
perform  a Monte Carlo simulation and then form a confidence interval for the
expected cost. This confidence interval is an estimate for the upper bound.

In addition to the confidence interval, we can calculate the lower bound using
[`Kokako.calculate_bound`](@ref).

```jldoctest tutorial_two; filter=r"Confidence interval.+"
julia> simulations = Kokako.simulate(model, 500);

julia> objective_values = [
           sum(stage[:stage_objective] for stage in sim) for sim in simulations
       ];

julia> using Statistics

julia> μ = round(mean(objective_values), digits = 2);

julia> ci = round(1.96 * std(objective_values) / sqrt(500), digits = 2);

julia> println("Confidence interval: ", μ, " ± ", ci)
Confidence interval: 8400.00 ± 409.34

julia> println("Lower bound: ", round(Kokako.calculate_bound(model), digits = 2))
Lower bound: 8333.33
```

In addition to simulating the primal values of variables, we can also pass
`SDDP.jl` custom recorder functions. Each of these functions takes one
argument, the JuMP subproblem, and returns anything you can compute. For example,
the dual of the demand constraint (which we named `demand_constraint`)
corresponds to the price we should charge for electricity, since it represents
the cost of each additional unit of demand. To calculate this, we can go

```jldoctest tutorial_two; filter = r"\s+?\-?\d+\.0"
julia> simulations = Kokako.simulate(
           model,
           1,
           custom_recorders = Dict{Symbol, Function}(
               :price => (sp) -> JuMP.dual(sp[:demand_constraint])
           )
       );

julia> [stage[:price] for stage in simulations[1]]
3-element Array{Float64,1}:
  50.0
 100.0
  -0.0
```

This concludes our second tutorial for `SDDP.jl`. In the next tutorial,
[Basic III: objective uncertainty](@ref), we extend the uncertainty to the
fuel cost.
