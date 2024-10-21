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
#  * $c$ is the unit cost for purchasing each unit
#  * $h$ is the holding cost per unit remaining at the end of each stage
#  * $p$ is the shortage cost per unit of unsatisfied demand at the end of each
#    stage
#
# There are no fixed ordering costs, and the demand at each stage is assumed to
# follow an independent and identically distributed random variable with
# cumulative distribution function (CDF) $\Phi(\cdot)$. Any unsatisfied demand
# is backlogged and carried forward to the next stage.

# At each stage, an agent must decide how many items to order or, equivalently,
# what the inventory level should be at the beginning of the next stage. The
# per-stage costs are the sum of the order costs, shortage and holding costs
# incurred at the end of the stage, after demand is realized. If $x$ represents
# the inventory level at the beginning of a stage, and $y$ is the level of
# inventory after ordering, the stage costs are given by
#
# ```math
# c(y-x) + \int_{0}^y h(y-\xi) d\Phi(\xi) + \int_{y}^{\infty} p(\xi - y) d\Phi(\xi).
# ```

# Following Chapter 19 of Introduction to Operations Research by Hillier and
# Lieberman (7th edition), we use the following parameters: $c=15, h=1, p=15$,
# and demand follows a continuous uniform distribution between 0 and 800.

x0 = 100        # Initial inventory
α = 0.995       # discount factor
c = 35          # unit inventory cost
h = 1           # unit inventory holding cost
p = 15          # unit order cost
UD = 800        # upper bound for demand
N = 10          # size of sample space
Ω = rand(Distributions.Uniform(0, UD), N);

# This is a well-known inventory problem with a closed-form solution. The
# optimal policy is a simple order-up-to policy: if the inventory level is
# below 741 units, the decision-maker orders up to 741 units. Otherwise, no
# order is placed. For a detailed analysis, refer to Foundations of Stochastic
# Inventory Theory by Evan Porteus (Stanford Business Books, 2002).

# ## Finite horizon

# For a finite horizon of length $T$,  the problem is to minimize the total
# expected cost over all stages, with a terminal value function given by:
# ```math
# v_T(x) = cx.
# ```
# This reflects that at the end of the last stage, the decision-maker can
# recover the unit cost for each leftover item while incurring any backlog costs
# and shortage penalties from the previous period.

T = 10 # number of stages

model = SDDP.LinearPolicyGraph(
    stages = T + 1,
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(sp, x, SDDP.State, initial_value = x0)
    @variable(sp, y, SDDP.State, initial_value = x0)
    @constraint(sp, y.out >= x.out)
    @constraint(sp, con_balance, x.out == y.in)
    if t == 1
        @stageobjective(sp, 0)
    elseif t == T + 1
        @stageobjective(sp, α^(t - 1) * (-c * x.in))
    else
        @stageobjective(
            sp,
            α^(t - 1) * (
                p * UD / 2 + (c - p) * y.in - c * x.in + (h + p) / (2 * UD) * y.in^2
            )
        )
        SDDP.parameterize(ω -> JuMP.set_normalized_rhs(con_balance, -ω), sp, Ω)
    end
    return
end

# Train and simulate the policy:

SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(model, 200, [:y])
objective_values = [sum(t[:stage_objective] for t in s) for s in simulations]
μ, ci = round.(SDDP.confidence_interval(objective_values, 1.96); digits = 2)
lower_bound = round(SDDP.calculate_bound(model); digits = 2)
println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", lower_bound)

# Plot the optimal inventory levels:

SDDP.publication_plot(
    simulations;
    title = "Inventory level",
    xlabel = "Stage",
    ylabel = "Inventory level",
) do data
    return data[:y].in
end

# ## Infinite horizon

# For the infinite horizon case, assume a discount factor $\alpha=0.995$. The
# objective in this case is to minimize the discounted expected costs over an
# infinite planning horizon.

model = SDDP.PolicyGraph(
    SDDP.UnicyclicGraph(α; num_nodes = 1);
    sense = :Min,
    lower_bound = 0.0,
    optimizer = HiGHS.Optimizer,
) do sp, t
    @variable(sp, x, SDDP.State, initial_value = x0)
    @variable(sp, y, SDDP.State, initial_value = x0)
    @constraint(sp, y.out >= x.out)
    @constraint(sp, con_balance, x.out == y.in)
    @stageobjective(
        sp,
        p * UD / 2 + (c - p) * y.in - c * x.in + (h + p) / (2 * UD) * y.in^2
    )
    SDDP.parameterize(ω -> JuMP.set_normalized_rhs(con_balance, -ω), sp, Ω)
    return
end

# Train and simulate the policy:

SDDP.train(model; iteration_limit = 100)
simulations = SDDP.simulate(model, 200, [:y])
objective_values = [sum(t[:stage_objective] for t in s) for s in simulations]
μ, ci = round.(SDDP.confidence_interval(objective_values, 1.96); digits = 2)
lower_bound = round(SDDP.calculate_bound(model); digits = 2)
println("Confidence interval: ", μ, " ± ", ci)
println("Lower bound: ", lower_bound)

# Plot the optimal inventory levels:

SDDP.publication_plot(
    simulations;
    title = "Inventory level",
    xlabel = "Stage",
    ylabel = "Inventory level",
) do data
    return data[:y].in
end
