# # Basic II: adding uncertainty

# In the previous tutorial, [Basic I: first steps](@ref), we created a
# deterministic  hydro-thermal scheduling model. In this tutorial, we extend the
# problem by adding uncertainty.

# Notably missing from our previous model were inflows. Inflows are the water
# that flows into the reservoir through rainfall or rivers. These inflows are
# uncertain, and are the cause of the main trade-off in hydro-thermal
# scheduling: the desire to use water now to generate cheap electricity, against
# the risk that future inflows will be low, leading to blackouts or expensive
# thermal generation.

# For our simple model, we assume that the inflows can be modelled by a discrete
# distribution with the three outcomes given in the following table:
#
# | ω    |   0 |  50 | 100 |
# | ---- | --- | --- | --- |
# | P(ω) | 1/3 | 1/3 | 1/3 |

# The value of the noise (the random variable) is observed by the agent at the
# start of each stage. This makes the problem a _wait-and-see_ or
# _hazard-decision_ formulation.

# To represent this, we can draw the following picture. The wavy lines denote
# the uncertainty arriving into the start of each stage (node).
#
# ![Linear policy graph](../assets/stochastic_linear_policy_graph.png)

# In addition to adding this uncertainty to the model, we also need to modify
# the dynamics to include `inflow`:
#
# `volume.out = volume.in + inflow - hydro_generation - hydro_spill`

# ## Creating a model

# To add an uncertain variable to the model, we create a new JuMP variable
# `inflow`, and then call the function [`SDDP.parameterize`](@ref). The
# [`SDDP.parameterize`](@ref) function takes three arguments: the subproblem, a
# vector of realizations, and a corresponding vector of probabilities.

using SDDP, GLPK

model = SDDP.LinearPolicyGraph(
    stages = 3,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer
) do subproblem, t
    ## Define the state variable.
    @variable(subproblem, 0 <= volume <= 200, SDDP.State, initial_value = 200)
    ## Define the control variables.
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
        inflow
    end)
    ## Define the constraints
    @constraints(subproblem, begin
        volume.out == volume.in + inflow - hydro_generation - hydro_spill
        demand_constraint, thermal_generation + hydro_generation == 150.0
    end)
    ## Define the objective for each stage `t`. Note that we can use `t` as an
    ## index for t = 1, 2, 3.
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(subproblem, fuel_cost[t] * thermal_generation)
    ## Parameterize the subproblem.
    SDDP.parameterize(subproblem, [0.0, 50.0, 100.0], [1/3, 1/3, 1/3]) do ω
        JuMP.fix(inflow, ω)
    end
end

# Note how we use the JuMP function
# [`JuMP.fix`](http://www.juliaopt.org/JuMP.jl/v0.19/variables/#JuMP.fix) to set
# the value of the `inflow` variable to `ω`.

# !!! note
#     [`SDDP.parameterize`](@ref) can only be called once in each subproblem
#     definition!

# ## Training and simulating the policy

# As in [Basic I: first steps](@ref), we train the policy:

SDDP.train(model; iteration_limit = 10)

# We can also simulate the policy. Note that this time, the simulation is
# stochastic. One common approach to quantify the quality of the policy is to
# perform  a Monte Carlo simulation and then form a confidence interval for the
# expected cost. This confidence interval is an estimate for the upper bound.

# In addition to the confidence interval, we can calculate the lower bound using
# [`SDDP.calculate_bound`](@ref).

using Statistics

simulations = SDDP.simulate(model, 500)

objective_values = [
    sum(stage[:stage_objective] for stage in sim) for sim in simulations
]

μ = round(mean(objective_values), digits = 2)

ci = round(1.96 * std(objective_values) / sqrt(500), digits = 2)

println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", round(SDDP.calculate_bound(model), digits = 2))

# In addition to simulating the primal values of variables, we can also pass
# `SDDP.jl` custom recorder functions. Each of these functions takes one
# argument, the JuMP subproblem, and returns anything you can compute. For
# example, the dual of the demand constraint (which we named
# `demand_constraint`) corresponds to the price we should charge for
# electricity, since it represents the cost of each additional unit of demand.
# To calculate this, we can go

simulations = SDDP.simulate(
    model,
    1,
    custom_recorders = Dict{Symbol, Function}(
        :price => (sp) -> JuMP.dual(sp[:demand_constraint])
    )
)

prices = [stage[:price] for stage in simulations[1]]

# This concludes our second tutorial for `SDDP.jl`. In the next tutorial, [Basic
# III: objective uncertainty](@ref), we extend the uncertainty to the fuel cost.
