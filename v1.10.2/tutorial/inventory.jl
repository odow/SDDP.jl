#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: inventory management

# The purpose of this tutorial is to demonstrate a well-known inventory
# management problem with a finite- and infinite-horizon policy.

# ## Required packages

# This tutorial requires the following packages:

using SDDP
import Distributions
import HiGHS
import Plots
import Statistics

# ## Background

# Consider a periodic review inventory problem involving a single product. The
# initial inventory is denoted by $x_0 \geq 0$, and a decision-maker can place
# an order at the start of each stage. The objective is to minimize expected
# costs over the planning horizon. The following parameters define the cost
# structure:
#
#  * `c` is the unit cost for purchasing each unit
#  * `h` is the holding cost per unit remaining at the end of each stage
#  * `p` is the shortage cost per unit of unsatisfied demand at the end of each
#    stage
#
# There are no fixed ordering costs, and the demand at each stage is assumed to
# follow an independent and identically distributed random variable with
# cumulative distribution function (CDF) $\Phi(\cdot)$. Any unsatisfied demand
# is backlogged and carried forward to the next stage.

# At each stage, an agent must decide how many items to order. The per-stage
# costs are the sum of the order costs, shortage and holding costs incurred at
# the end of the stage, after demand is realized.

# Following Chapter 19 of Introduction to Operations Research by Hillier and
# Lieberman (7th edition), we use the following parameters: $c=15, h=1, p=15$.

x_0 = 10        # initial inventory
c = 35          # unit inventory cost
h = 1           # unit inventory holding cost
p = 15          # unit order cost

# Demand follows a continuous uniform distribution between 0 and 800. We
# construct a sample average approximation with 20 scenarios:

Ω = range(0, 800; length = 20);

# This is a well-known inventory problem with a closed-form solution. The
# optimal policy is a simple order-up-to policy: if the inventory level is
# below a certain number of units, the decision-maker orders up to that number
# of units. Otherwise, no order is placed. For a detailed analysis, refer
# to Foundations of Stochastic Inventory Theory by Evan Porteus (Stanford
# Business Books, 2002).

# ## Finite horizon

# For a finite horizon of length $T$,  the problem is to minimize the total
# expected cost over all stages.

# In the last stage, the decision-maker can recover the unit cost `c` for each
# leftover item, or buy out any remaining backlog, also at the unit cost `c`.

T = 10 # number of stages
model = SDDP.LinearPolicyGraph(;
    stages = T + 1,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(sp, x_inventory >= 0, SDDP.State, initial_value = x_0)
    @variable(sp, x_demand >= 0, SDDP.State, initial_value = 0)
    ## u_buy is a Decision-Hazard control variable. We decide u.out for use in
    ## the next stage
    @variable(sp, u_buy >= 0, SDDP.State, initial_value = 0)
    @variable(sp, u_sell >= 0)
    @variable(sp, w_demand == 0)
    @constraint(sp, x_inventory.out == x_inventory.in + u_buy.in - u_sell)
    @constraint(sp, x_demand.out == x_demand.in + w_demand - u_sell)
    if t == 1
        fix(u_sell, 0; force = true)
        @stageobjective(sp, c * u_buy.out)
    elseif t == T + 1
        fix(u_buy.out, 0; force = true)
        @stageobjective(sp, -c * x_inventory.out + c * x_demand.out)
        SDDP.parameterize(ω -> JuMP.fix(w_demand, ω), sp, Ω)
    else
        @stageobjective(sp, c * u_buy.out + h * x_inventory.out + p * x_demand.out)
        SDDP.parameterize(ω -> JuMP.fix(w_demand, ω), sp, Ω)
    end
    return
end

# Train and simulate the policy:

SDDP.train(model)
simulations = SDDP.simulate(model, 200, [:x_inventory, :u_buy])
objective_values = [sum(t[:stage_objective] for t in s) for s in simulations]
μ, ci = round.(SDDP.confidence_interval(objective_values, 1.96); digits = 2)
lower_bound = round(SDDP.calculate_bound(model); digits = 2)
println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", lower_bound)

# Plot the optimal inventory levels:

plt = SDDP.publication_plot(
    simulations;
    title = "x_inventory.out + u_buy.out",
    xlabel = "Stage",
    ylabel = "Quantity",
    ylims = (0, 1_000),
) do data
    return data[:x_inventory].out + data[:u_buy].out
end

# In the early stages, we indeed recover an order-up-to policy. However,
# there are end-of-horizon effects as the agent tries to optimize their
# decision making knowing that they have 10 realizations of demand.

# ## Infinite horizon

# We can remove the end-of-horizonn effects by considering an infinite
# horizon model. We assume a discount factor $\alpha=0.95$:

α = 0.95
graph = SDDP.LinearGraph(2)
SDDP.add_edge(graph, 2 => 2, α)
graph

# The objective in this case is to minimize the discounted expected costs over
# an infinite planning horizon.

model = SDDP.PolicyGraph(
    graph;
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(sp, x_inventory >= 0, SDDP.State, initial_value = x_0)
    @variable(sp, x_demand >= 0, SDDP.State, initial_value = 0)
    ## u_buy is a Decision-Hazard control variable. We decide u.out for use in
    ## the next stage
    @variable(sp, u_buy >= 0, SDDP.State, initial_value = 0)
    @variable(sp, u_sell >= 0)
    @variable(sp, w_demand == 0)
    @constraint(sp, x_inventory.out == x_inventory.in + u_buy.in - u_sell)
    @constraint(sp, x_demand.out == x_demand.in + w_demand - u_sell)
    if t == 1
        fix(u_sell, 0; force = true)
        @stageobjective(sp, c * u_buy.out)
    else
        @stageobjective(sp, c * u_buy.out + h * x_inventory.out + p * x_demand.out)
        SDDP.parameterize(ω -> JuMP.fix(w_demand, ω), sp, Ω)
    end
    return
end

SDDP.train(model; iteration_limit = 400)
simulations = SDDP.simulate(
    model,
    200,
    [:x_inventory, :u_buy];
    sampling_scheme = SDDP.InSampleMonteCarlo(;
        max_depth = 50,
        terminate_on_dummy_leaf = false,
    ),
);

# Plot the optimal inventory levels:

plt = SDDP.publication_plot(
    simulations;
    title = "x_inventory.out + u_buy.out",
    xlabel = "Stage",
    ylabel = "Quantity",
    ylims = (0, 1_000),
) do data
    return data[:x_inventory].out + data[:u_buy].out
end
Plots.hline!(plt, [662]; label = "Analytic solution")

# We again recover an order-up-to policy. The analytic solution is to
# order-up-to 662 units. We do not precisely recover this solution because
# we used a sample average approximation of 20 elements. If we increased the
# number of samples, our solution would approach the analytic solution.
