
#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: deterministic to stochastic

# The purpose of this tutorial is to explain how we can go from a deterministic
# time-staged optimal control model in JuMP to a multistage stochastic
# optimization model in `SDDP.jl`. As a motivating problem, we consider the
# hydro-thermal problem with a single reservoir.

# _This tutorial was written by Oscar Dowson and Andy Philpott for the 2023
# Stochastic Programming Winter School, held March XX to XX in XXX, Norway._

# ## Packages

# This tutorial requires the following packages:

using JuMP
using SDDP

import CSV
import DataFrames
import HiGHS
import Plots

# ## Data

# The data for this tutorial is contained in the `example_reservoir.csv` file in
# the SDDP.jl repository. To run locally, save the CSV file, then change
# `filename` to point to the location where you downloaded it to.

filename = joinpath(@__DIR__, "example_reservoir.csv")
data = CSV.read("data52.csv", DataFrames.DataFrame)

# It's easier to visualize if we plot it:

Plots.plot(
    Plots.plot(data[!, :inflow], ylabel = "Inflow"),
    Plots.plot(data[!, :demand], ylabel = "Demand"),
    Plots.plot(data[!, :cost], ylabel = "Cost", xlabel = "Week");
    layout = (3, 1),
    legend = false,
)

# The number of weeks will be useful later:

T = size(data, 1)

# ## Deterministic JuMP model

# To start, we construct a deterministic model in pure JuMP.

# Create a JuMP model, using `HiGHS` as the optimizer:

model = Model(HiGHS.Optimizer)
set_silent(model)

# `x[t]`: the amount of water in the reservoir at the start of stage `t`:

reservoir_max = 320.0

@variable(model, 0 <= x[1:T+1] <= reservoir_max)

# We need an initial condition for `x[1]`. Fix it to 300 units:

reservoir_initial = 300

fix(x[1], reservoir_initial; force = true)

# `u[t]`: the amount of water to flow through the turbine in stage `t`:

flow_max = 12

@variable(model, 0 <= u[i = 1:T] <= flow_max)

# `s[t]`: the amount of water to spill from the reservoir in stage `t`,
# bypassing the turbine:

@variable(model, 0 <= s[1:T])

# `r[t]`: the amount of thermal generation in stage `t`:
@variable(model, 0 <= r[1:T])

# `i[t]`: the amount of inflow to the reservoir in stage `t`:

@variable(model, i[1:T])

# For this model, our inflow is fixed, so we fix it to the data we have:

for t in 1:T
    fix(i[t], data[t, :inflow])
end

# The water balance constraint says that the water in the reservoir at the start
# of stage `t+1` is the water in the reservoir at the start of stage `t`, less
# the amount flowed through the turbine, `u[t]`, less the amount spilled,
# `s[t]`, plus the amount of inflow, `i[t]`, into the reservoir:

@constraint(model, [t = 1:T], x[t+1] == x[t] - u[t] - s[t] + i[t])

# We also need a `supply = demand` constraint. In practice, the units of this
# would be in MWh, and there would be a conversion factor between the amount of
# water flowing through the turbine and the power output. To simplify, we assume
# that power and water have the same units, so that one "unit" of demand is
# equal to one "unit" of the reservoir `x[t]`:

@constraint(model, [t = 1:T], u[t] + r[t] == data[t, :demand])

# Our objective is to minimize the cost of thermal generation:

@objective(model, Min, sum(data[t, :cost] * r[t] for t in 1:T))

# Let's optimize and check the solution

optimize!(model)
solution_summary(model)

# The total cost is:

objective_value(model)

# Here's a plot of demand and generation:

Plots.plot(data[!, :demand]; label = "Demand", xlabel = "Week")
Plots.plot!(value.(r), label = "Thermal")
Plots.plot!(value.(u), label = "Hydro")

# And here's the storage over time:

Plots.plot(value.(x); label = "Storage", xlabel = "Week")

# ## Deterministic SDDP model

# For the next step, we show how to decompose our JuMP model into SDDP.jl. It
# should obtain the same solution!

model = SDDP.LinearPolicyGraph(
    stages = T,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(
        sp,
        0 <= x <= reservoir_max,
        SDDP.State,
        initial_value = reservoir_initial,
    )
    @variable(sp, 0 <= u <= flow_max)
    @variable(sp, 0 <= r)
    @variable(sp, 0 <= s)
    @variable(sp, i)
    fix(i, data[t, :inflow])
    @constraint(sp, x.out == x.in - u - s + i)
    @constraint(sp, u + r == data[t, :demand])
    @stageobjective(sp, data[t, :cost] * r)
end

# Can you see the differences?

# Instead of calling `JuMP.optimize!`, SDDP.jl uses a `train` method. With our
# machine learning hat on, you can think of SDDP.jl as training a function for
# each stage that accepts the current reservoir state as input and returns the
# optimal actions as output. It is also an iterative algorithm, so we need to
# specify when it should terminate:

SDDP.train(model; iteration_limit = 10)

# As a quick sanity check, did we get the same cost as our JuMP model?

SDDP.calculate_bound(model)

# That's good. Next, to check the value of the decision variables. This isn't as
# straight forward as our JuMP model. Instead, we need to _simulate_ the policy,
# and then extract the values of the decision variables from the results of the
# simulation.

# Since our model is determministic, we need only 1 replication of the
# simulation, and we want to record the values of the `x`, `u`, and `r`
# variables:

n_replications = 1
simulations = SDDP.simulate(model, n_replications, [:x, :u, :r]);

# The `simulations` vector is too big to show. But it contains one element for
# each replication, and each replication contains one dictionary for each stage.

# For example, the data corresponding to the tenth stage in the first
# replication is:

simulations[1][10]

# Let's grab the trace of the `r` and `u` variables in the first replication,
# and then plot them:

r_sim = [sim[:r] for sim in simulations[1]]
u_sim = [sim[:u] for sim in simulations[1]]

Plots.plot(data[!, :demand]; label = "Demand", xlabel = "Week")
Plots.plot!(r_sim, label = "Thermal")
Plots.plot!(u_sim, label = "Hydro")

# Perfect. That's the same as we got before.

# Now let's look at `x`. This is a little more complicated, because we need to
# grab the outgoing value of the state variable in each stage:

x_sim = [sim[:x].out for sim in simulations[1]]

Plots.plot(x_sim; label = "Storage", xlabel = "Week")

# ## Stochastic SDDP model

model = SDDP.LinearPolicyGraph(
    stages = T,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(
        sp,
        0 <= x <= reservoir_max,
        SDDP.State,
        initial_value = reservoir_initial,
    )
    @variable(sp, 0 <= u <= flow_max)
    @variable(sp, 0 <= r)
    @variable(sp, 0 <= s)
    @variable(sp, i)
    ## <--- This bit is new
    Ω = [data[t, :inflow] - 2, data[t, :inflow], data[t, :inflow] + 5]
	P = [0.3, 0.4, 0.3]
    SDDP.parameterize(sp, Ω, P) do ω
        fix(i,ω)
    end
    ## --->
    @constraint(sp, x.out == x.in - u - s + i)
    @constraint(sp, u + r == data[t, :demand])
    @stageobjective(sp, data[t, :cost] * r)
end

# Can you see the differences?

# Let's train our new model. We need more iterations because of the
# stochasticity:

SDDP.train(model; iteration_limit = 100)

# Now simulate the policy. This time we do 100 replications because the policy
# is now stochastic instead of deterministic:

n_replications = 100
simulations = SDDP.simulate(model, n_replications, [:x, :u, :r, :i]);

# And let's plot the use of thermal generation in each replication:

plt = Plots.plot(data[!, :demand]; label = "Demand", xlabel = "Week")
for simulation in simulations
    Plots.plot!(plt, [sim[:r] for sim in simulation], label = "")
end
plt

# Viewing an interpreting static plots like this is difficult, particularly as
# the number of simulations grows. SDDP.jl includes an interactive
# `SpaghettiPlot` that makes things easier:

plt = SDDP.SpaghettiPlot(simulations)
SDDP.add_spaghetti(plt; title = "Storage") do sim
   return sim[:x].out
end
SDDP.add_spaghetti(plt; title = "Hydro") do sim
    return sim[:u]
end
SDDP.add_spaghetti(plt; title = "Inflow") do sim
    return sim[:i]
end
SDDP.plot(plt, "spaghetti_plot.html")

# ## Next steps

# Take this model and modify it following the suggestions below. The SDDP.jl
# documentation has a range of similar examples and hints for how to achieve
# them.

# * Can you add a second reservoir to make a river chain?
# * Can you add random demand or cost data as well as inflows?
# * Can you add a risk measure to make the policy risk-averse?
# * Was our stopping rule correct? What happens if we use fewer or more
#   iterations? What other stopping rules could you try?
# * The model ends with an empty reservoir. That isn't ideal for the following
#   year. Can you solve an infinite horizon model instead of finite horizon?
