
#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: deterministic to stochastic

# The purpose of this tutorial is to explain how we can go from a deterministic
# time-staged optimal control model in JuMP to a multistage stochastic
# optimization model in `SDDP.jl`. As a motivating problem, we consider the
# hydro-thermal problem with a single reservoir.

# ## Packages

# This tutorial requires the following packages:

using JuMP
using SDDP
import CSV
import DataFrames
import HiGHS
import Plots

# ## Data

# First, we need some data for the problem. For this tutorial, we'll write CSV
# files to a temporary directory from Julia. If you have an existing file, you
# could change the filename to point to that instead.

dir = mktempdir()
filename = joinpath(dir, "example_reservoir.csv")

# Here is the data

csv_data = """
week,inflow,demand,cost
1,3,7,10.2\n2,2,7.1,10.4\n3,3,7.2,10.6\n4,2,7.3,10.9\n5,3,7.4,11.2\n
6,2,7.6,11.5\n7,3,7.8,11.9\n8,2,8.1,12.3\n9,3,8.3,12.7\n10,2,8.6,13.1\n
11,3,8.9,13.6\n12,2,9.2,14\n13,3,9.5,14.5\n14,2,9.8,14.9\n15,3,10.1,15.3\n
16,2,10.4,15.8\n17,3,10.7,16.2\n18,2,10.9,16.6\n19,3,11.2,17\n20,3,11.4,17.4\n
21,3,11.6,17.7\n22,2,11.7,18\n23,3,11.8,18.3\n24,2,11.9,18.5\n25,3,12,18.7\n
26,2,12,18.9\n27,3,12,19\n28,2,11.9,19.1\n29,3,11.8,19.2\n30,2,11.7,19.2\n
31,3,11.6,19.2\n32,2,11.4,19.2\n33,3,11.2,19.1\n34,2,10.9,19\n35,3,10.7,18.9\n
36,2,10.4,18.8\n37,3,10.1,18.6\n38,2,9.8,18.5\n39,3,9.5,18.4\n40,3,9.2,18.2\n
41,2,8.9,18.1\n42,3,8.6,17.9\n43,2,8.3,17.8\n44,3,8.1,17.7\n45,2,7.8,17.6\n
46,3,7.6,17.5\n47,2,7.4,17.5\n48,3,7.3,17.5\n49,2,7.2,17.5\n50,3,7.1,17.6\n
51,3,7,17.7\n52,3,7,17.8\n
"""
write(filename, csv_data);

# And here we read it into a DataFrame:

data = CSV.read(filename, DataFrames.DataFrame)

# It's easier to visualize the data if we plot it:

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

# `x_storage[t]`: the amount of water in the reservoir at the start of stage `t`:

reservoir_max = 320.0
@variable(model, 0 <= x_storage[1:T+1] <= reservoir_max)

# We need an initial condition for `x_storage[1]`. Fix it to 300 units:

reservoir_initial = 300
fix(x_storage[1], reservoir_initial; force = true)

# `u_flow[t]`: the amount of water to flow through the turbine in stage `t`:

flow_max = 12
@variable(model, 0 <= u_flow[1:T] <= flow_max)

# `u_spill[t]`: the amount of water to spill from the reservoir in stage `t`,
# bypassing the turbine:

@variable(model, 0 <= u_spill[1:T])

# `u_thermal[t]`: the amount of thermal generation in stage `t`:

@variable(model, 0 <= u_thermal[1:T])

# `ω_inflow[t]`: the amount of inflow to the reservoir in stage `t`:

@variable(model, ω_inflow[1:T])

# For this model, our inflow is fixed, so we fix it to the data we have:

for t in 1:T
    fix(ω_inflow[t], data[t, :inflow])
end

# The water balance constraint says that the water in the reservoir at the start
# of stage `t+1` is the water in the reservoir at the start of stage `t`, less
# the amount flowed through the turbine, `u_flow[t]`, less the amount spilled,
# `u_spill[t]`, plus the amount of inflow, `ω_inflow[t]`, into the reservoir:

@constraint(
    model,
    [t in 1:T],
    x_storage[t+1] == x_storage[t] - u_flow[t] - u_spill[t] + ω_inflow[t],
)

# We also need a `supply = demand` constraint. In practice, the units of this
# would be in MWh, and there would be a conversion factor between the amount of
# water flowing through the turbine and the power output. To simplify, we assume
# that power and water have the same units, so that one "unit" of demand is
# equal to one "unit" of the reservoir `x_storage[t]`:

@constraint(model, [t in 1:T], u_flow[t] + u_thermal[t] == data[t, :demand])

# Our objective is to minimize the cost of thermal generation:

@objective(model, Min, sum(data[t, :cost] * u_thermal[t] for t in 1:T))

# Let's optimize and check the solution

optimize!(model)
solution_summary(model)

# The total cost is:

objective_value(model)

# Here's a plot of demand and generation:

Plots.plot(data[!, :demand]; label = "Demand", xlabel = "Week")
Plots.plot!(value.(u_thermal), label = "Thermal")
Plots.plot!(value.(u_flow), label = "Hydro")

# And here's the storage over time:

Plots.plot(value.(x_storage); label = "Storage", xlabel = "Week")

# ## Deterministic SDDP model

# For the next step, we show how to decompose our JuMP model into SDDP.jl. It
# should obtain the same solution.

model = SDDP.LinearPolicyGraph(
    stages = T,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(
        sp,
        0 <= x_storage <= reservoir_max,
        SDDP.State,
        initial_value = reservoir_initial,
    )
    @variable(sp, 0 <= u_flow <= flow_max)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    fix(ω_inflow, data[t, :inflow])
    @constraint(sp, x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow)
    @constraint(sp, u_flow + u_thermal == data[t, :demand])
    @stageobjective(sp, data[t, :cost] * u_thermal)
    return
end

# Can you see how the JuMP model maps to this syntax? We have created a
# [`SDDP.LinearPolicyGraph`](@ref) with `T` stages, we're minimizing, and we're
# using `HiGHS.Optimizer` as the optimizer.

# A few bits might be non-obvious:
#
# * We need to provide a lower bound for the objective function. Since our costs
#   are always positive, a valid lower bound for the total cost is `0.0`.
# * We define `x_storage` as a state variable using `SDDP.State`. A state
#   variable is any variable that flows through time, and for which we need to
#   know the value of it in stage `t-1` to compute the best action in stage `t`.
#   The state variable `x_storage` is actually two decision variables,
#   `x_storage.in` and `x_storage.out`, which represent `x_storage[t]` and
#   `x_storage[t+1]` respectively.
# * We need to use `@stageobjective` instead of `@objective`.

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

# Since our model is deterministic, we need only 1 replication of the
# simulation, and we want to record the values of the `x_storage`, `u_flow`, and
# `u_thermal` variables:

simulations = SDDP.simulate(
    model,
    1,  # Number of replications
    [:x_storage, :u_flow, :u_thermal],
);

# The `simulations` vector is too big to show. But it contains one element for
# each replication, and each replication contains one dictionary for each stage.

# For example, the data corresponding to the tenth stage in the first
# replication is:

simulations[1][10]

# Let's grab the trace of the `u_thermal` and `u_flow` variables in the first
# replication, and then plot them:

r_sim = [sim[:u_thermal] for sim in simulations[1]]
u_sim = [sim[:u_flow] for sim in simulations[1]]

Plots.plot(data[!, :demand]; label = "Demand", xlabel = "Week")
Plots.plot!(r_sim, label = "Thermal")
Plots.plot!(u_sim, label = "Hydro")

# Perfect. That's the same as we got before.

# Now let's look at `x_storage`. This is a little more complicated, because we
# need to grab the outgoing value of the state variable in each stage:

x_sim = [sim[:x_storage].out for sim in simulations[1]]

Plots.plot(x_sim; label = "Storage", xlabel = "Week")

# ## Stochastic SDDP model

# Now we add some randomness to our model. In each stage, we assume that the
# inflow could be: 2 units lower, with 30% probability; the same as before, with
# 40% probability; or 5 units higher, with 30% probability.

model = SDDP.LinearPolicyGraph(
    stages = T,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(
        sp,
        0 <= x_storage <= reservoir_max,
        SDDP.State,
        initial_value = reservoir_initial,
    )
    @variable(sp, 0 <= u_flow <= flow_max)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    ## <--- This bit is new
    Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
    SDDP.parameterize(sp, Ω, P) do ω
        fix(ω_inflow, data[t, :inflow] + ω)
        return
    end
    ## --->
    @constraint(sp, x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow)
    @constraint(sp, u_flow + u_thermal == data[t, :demand])
    @stageobjective(sp, data[t, :cost] * u_thermal)
    return
end

# Can you see the differences?

# Let's train our new model. We need more iterations because of the
# stochasticity:

SDDP.train(model; iteration_limit = 100)

# Now simulate the policy. This time we do 100 replications because the policy
# is now stochastic instead of deterministic:

simulations =
    SDDP.simulate(model, 100, [:x_storage, :u_flow, :u_thermal, :ω_inflow]);

# And let's plot the use of thermal generation in each replication:

plot = Plots.plot(data[!, :demand]; label = "Demand", xlabel = "Week")
for simulation in simulations
    Plots.plot!(plot, [sim[:u_thermal] for sim in simulation], label = "")
end
plot

# Viewing an interpreting static plots like this is difficult, particularly as
# the number of simulations grows. SDDP.jl includes an interactive
# `SpaghettiPlot` that makes things easier:

plot = SDDP.SpaghettiPlot(simulations)
SDDP.add_spaghetti(plot; title = "Storage") do sim
    return sim[:x_storage].out
end
SDDP.add_spaghetti(plot; title = "Hydro") do sim
    return sim[:u_flow]
end
SDDP.add_spaghetti(plot; title = "Inflow") do sim
    return sim[:ω_inflow]
end
SDDP.plot(
    plot,
    "spaghetti_plot.html";
    ## We need this to build the documentation. Set to true if running locally.
    open = false,
)

# ```@raw html
# <iframe src="../spaghetti_plot.html" style="width:100%;height:500px;"></iframe>
# ```

# !!! info
#     If you have trouble viewing the plot, you can
#     [open it in a new window](spaghetti_plot.html).

# ## Cyclic graphs

# One major problem with our model is that the reservoir is empty at the end of
# the time horizon. This is because our model does not consider the cost of
# future years after the `T` weeks.

# We can fix this using a cyclic policy graph. One way to construct a graph is
# with the [`SDDP.UnicyclicGraph`](@ref) constructor:

SDDP.UnicyclicGraph(0.7; num_nodes = 2)

# This graph has two nodes, and a loop from node 2 back to node 1 with
# probability 0.7.

# We can construct a cyclic policy graph as follows:

graph = SDDP.UnicyclicGraph(0.95; num_nodes = T)
model = SDDP.PolicyGraph(
    graph;
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(
        sp,
        0 <= x_storage <= reservoir_max,
        SDDP.State,
        initial_value = reservoir_initial,
    )
    @variable(sp, 0 <= u_flow <= flow_max)
    @variable(sp, 0 <= u_thermal)
    @variable(sp, 0 <= u_spill)
    @variable(sp, ω_inflow)
    Ω, P = [-2, 0, 5], [0.3, 0.4, 0.3]
    SDDP.parameterize(sp, Ω, P) do ω
        fix(ω_inflow, data[t, :inflow] + ω)
        return
    end
    @constraint(sp, x_storage.out == x_storage.in - u_flow - u_spill + ω_inflow)
    @constraint(sp, u_flow + u_thermal == data[t, :demand])
    @stageobjective(sp, data[t, :cost] * u_thermal)
    return
end

# Notice how the only thing that has changed is our graph; the subproblems
# remain the same.

# Let's train a policy:

SDDP.train(model; iteration_limit = 100)

# When we simulate now, each trajectory will be a different length, because
# each cycle has a 95% probability of continuing and a 5% probability of
# stopping.

simulations = SDDP.simulate(model, 3);
length.(simulations)

# We can simulate a fixed number of cycles by passing a `sampling_scheme`:

simulations = SDDP.simulate(
    model,
    100,
    [:x_storage, :u_flow];
    sampling_scheme = SDDP.InSampleMonteCarlo(;
        max_depth = 5 * T,
        terminate_on_dummy_leaf = false,
    ),
);
length.(simulations)

# Let's visualize the policy:

Plots.plot(
    SDDP.publication_plot(simulations; ylabel = "Storage") do sim
        return sim[:x_storage].out
    end,
    SDDP.publication_plot(simulations; ylabel = "Hydro") do sim
        return sim[:u_flow]
    end;
    layout = (2, 1),
)

# ## Next steps

# Our model is very basic. There are many aspects that we could improve:
#
# * Can you add a second reservoir to make a river chain?
#
# * Can you modify the problem and data to use proper units, including a
#   conversion between the volume of water flowing through the turbine and the
#   electrical power output?
