#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: batteries

# The purpose of this tutorial is to demonstrate how to model the operation of a
# battery in power system with uncertain load.

# This example was originally developed by Andy Philpott as part of the
# Winter School Kvitfjell 2025, _Planning Under Uncertainty in Energy Markets_,
# organized by NORDAB - The Norwegian Doctoral Academy in Business.

# ## Packages

# This tutorial requires the following packages:

using SDDP
import HiGHS
import Plots

# ## The model

# Our model is a day in the life of a simple power system with one thermal
# generator and one battery. There are 24 hourly stages that form a linear
# policy graph.

# There are two state variables: `x_soc`, the state-of-charge of the battery;
# and `x_thermal`, the output of the thermal generator.

# There are four control variables: `u_charge` and `u_discharge` for charging
# and discharging the battery in each stage, and `u_slack` and `u_surplus` for
# measuring the slack and surplus energy generation.

# There is one random variable: `w_load`, which is the energy load of the system
# in each stage.

# The objective is to minimize cost, which is comprised of the cost of thermal
# generation and the value of lost-load (`u_slack`).

# There are three constraints: a balance on the state-of-charge of the battery,
# which accounts for charging inefficiency; ramping limits on the thermal
# generator; and an energy balance constraint on the energy produced.

# Here's the model in SDDP.jl:

model = SDDP.LinearPolicyGraph(;
    stages = 24,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    ## State variables
    @variable(sp, 0 <= x_soc <= 30, SDDP.State, initial_value = 4)
    @variable(sp, 0 <= x_thermal <= 70, SDDP.State, initial_value = 35)
    ## Control variables
    @variable(sp, 0 <= u_charge <= 15)
    @variable(sp, 0 <= u_discharge <= 15)
    @variable(sp, u_slack >= 0)
    @variable(sp, u_surplus >= 0)
    ## Random variables
    @variable(sp, w_load)
    Ω = [-4.0, -2.0, 0.0, 2.0, 4.0]
    d = vcat(
        [40, 41, 42, 43, 35, 40, 40, 25, 10, 8, 6, 5],  # Hours 01-12
        [5, 6, 8, 10, 20, 30, 55, 72, 75, 70, 64, 60],  # Hours 13-24
    )
    SDDP.parameterize(ω -> JuMP.fix(w_load, d[t] + ω), sp, Ω)
    ## Objective function
    @stageobjective(sp, 70 * x_thermal.out + 500 * u_slack)
    ## Constraints
    @constraint(sp, x_soc.out == x_soc.in + 0.8 * u_charge - u_discharge)
    @constraint(sp, x_thermal.out - x_thermal.in <= 10)
    @constraint(
        sp,
        λ,
        x_thermal.out + u_discharge - u_charge + u_slack == w_load + u_surplus,
    )
    return
end

# ## Training

# We train the model for 500 iterations using the [`SDDP.Threaded`](@ref)
# parallel scheme.

SDDP.train(model; iteration_limit = 500, parallel_scheme = SDDP.Threaded())

results = SDDP.simulate(
    model,
    100,
    [:x_soc, :x_thermal, :u_slack, :w_load];
    custom_recorders = Dict{Symbol,Function}(:λ => sp -> JuMP.dual(sp[:λ])),
);

# ## Analyzing the solution

# The load follows the classic "duck curve". There are two defining features:
# large amounts of solar during the middle of the day reduce the net load to
# form the belly of the duck, and the drop in solar combined with high gross
# load in the evening forms a steep ramp along the neck of the duck.

SDDP.publication_plot(results; xlabel = "Hour", ylabel = "Net load") do d
    return d[:w_load]
end

# Thermal generation mostly follows net load, except the thermal generation is
# limited by the ramp up constraints during the afternoon as the generator
# increases to reach maximum production during the 21 hour:

SDDP.publication_plot(
    results;
    xlabel = "Hour",
    ylabel = "Thermal generation",
) do d
    return d[:x_thermal].out
end

# The battery starts with 4 units of charge, empties it during the 10th hour,
# charges during the afternoon, and then discharges during the evening peak.
# Because we have modeled a finite horizon problem with no terminal value
# function, the battery ends the 24th hour with an empty state of charge.

SDDP.publication_plot(results; xlabel = "Hour", ylabel = "State of charge") do d
    return d[:x_soc].out
end

# The batteries charging decisions are driven by the dual value of the
# `demand_constraint` which can be interpreted as the price of energy. During
# the first part of the day, the price of energy is \$70/unit, which is the
# fuel price of the thermal generator. During the afternoon, the price
# decreases below the fuel price, and the battery arbitrages this difference by
# charging. The peak price is incurred during the evening, when the price
# sometimes spikes to \$500/unit, which is the value of lost load (VOLL).

SDDP.publication_plot(results; xlabel = "Hour", ylabel = "Dual λ") do d
    return d[:λ]
end

# The `u_slack` variable is the lost load:

SDDP.publication_plot(results; xlabel = "Hour", ylabel = "Lost load") do d
    return d[:u_slack]
end

# ## Next steps

# Can you use the information in the previous plots to explain why the energy
# price falls below the fuel cost?
