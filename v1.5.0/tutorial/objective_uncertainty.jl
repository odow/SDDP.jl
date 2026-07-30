#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Uncertainty in the objective function

# In the previous tutorial, [An introduction to SDDP.jl](@ref), we created a
# stochastic hydro-thermal scheduling model. In this tutorial, we extend the
# problem by adding uncertainty to the fuel costs.

# Previously, we assumed that the fuel cost was deterministic: \$50/MWh in the
# first stage, \$100/MWh in the second stage, and \$150/MWh in the third
# stage. For this tutorial, we assume that in addition to these base costs, the
# actual fuel cost is correlated with the inflows.

# Our new model for the uncertainty is given by the following table:
#
# | ω               |   1 |   2 |    3 |
# | ----            | --- | --- | ---- |
# | P(ω)            | 1/3 | 1/3 |  1/3 |
# | inflow          |   0 |  50 |  100 |
# | fuel multiplier | 1.5 | 1.0 | 0.75 |

# In stage `t`, the objective is now to minimize:
#
# `fuel_multiplier * fuel_cost[t] * thermal_generation`

# ## Creating a model

# To add an uncertain objective, we can simply call [`@stageobjective`](@ref)
# from inside the [`SDDP.parameterize`](@ref) function.

using SDDP, HiGHS

model = SDDP.LinearPolicyGraph(
    stages = 3,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do subproblem, t
    ## Define the state variable.
    @variable(subproblem, 0 <= volume <= 200, SDDP.State, initial_value = 200)
    ## Define the control variables.
    @variables(subproblem, begin
        thermal_generation >= 0
        hydro_generation >= 0
        hydro_spill >= 0
        inflow
    end)
    ## Define the constraints
    @constraints(
        subproblem,
        begin
            volume.out == volume.in + inflow - hydro_generation - hydro_spill
            thermal_generation + hydro_generation == 150.0
        end
    )
    fuel_cost = [50.0, 100.0, 150.0]
    ## Parameterize the subproblem.
    Ω = [
        (inflow = 0.0, fuel_multiplier = 1.5),
        (inflow = 50.0, fuel_multiplier = 1.0),
        (inflow = 100.0, fuel_multiplier = 0.75),
    ]
    SDDP.parameterize(subproblem, Ω, [1 / 3, 1 / 3, 1 / 3]) do ω
        JuMP.fix(inflow, ω.inflow)
        @stageobjective(
            subproblem,
            ω.fuel_multiplier * fuel_cost[t] * thermal_generation
        )
    end
end

# ## Training and simulating the policy

# As in the previous two tutorials, we train and simulate the policy:

SDDP.train(model; stopping_rules = [SDDP.SimulationStoppingRule()])

simulations = SDDP.simulate(model, 500)

objective_values =
    [sum(stage[:stage_objective] for stage in sim) for sim in simulations]

using Statistics

μ = round(mean(objective_values), digits = 2)
ci = round(1.96 * std(objective_values) / sqrt(500), digits = 2)

println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", round(SDDP.calculate_bound(model), digits = 2))
