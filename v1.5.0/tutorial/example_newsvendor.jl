#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: Two-stage Newsvendor

# This example is based on the classical newsvendor problem.

using SDDP
import HiGHS
import Plots

# First, we need some discretized distribution of demand. For simplicity, we're
# going to sample 10 points from the uniform `[0, 1)` distribution three times
# and add them together. This is a rough approximation of a normal distribution.

Ω = rand(10) .+ rand(10) .+ rand(10)

# We also need some price data. We assume the agent can buy a newspaper for \$1
# and sell it for \$1.20.

buy_price, sales_price = 1.0, 1.2

# Now we can formulate and train a policy for the two-stage newsvendor problem:

model = SDDP.LinearPolicyGraph(
    stages = 2,
    sense = :Max,
    upper_bound = maximum(Ω) * sales_price,
    optimizer = HiGHS.Optimizer,
) do subproblem, stage
    @variable(subproblem, inventory >= 0, SDDP.State, initial_value = 0)
    if stage == 1
        @variable(subproblem, buy >= 0)
        @constraint(subproblem, inventory.out == inventory.in + buy)
        @stageobjective(subproblem, -buy_price * buy)
    else
        @variable(subproblem, sell >= 0)
        @constraint(subproblem, sell <= inventory.in)
        SDDP.parameterize(subproblem, Ω) do ω
            return JuMP.set_upper_bound(sell, ω)
        end
        @stageobjective(subproblem, sales_price * sell)
    end
end

SDDP.train(model; stopping_rules = [SDDP.SimulationStoppingRule()])

# To check the first-stage buy decision, we need to obtain a decision rule for
# the first-stage node `1`:

first_stage_rule = SDDP.DecisionRule(model, node = 1)

# Then we can evaluate it, passing in a starting point for the incoming state:

solution = SDDP.evaluate(
    first_stage_rule;
    incoming_state = Dict(:inventory => 0.0),
    controls_to_record = [:buy],
)

# The optimal value of the `buy` variable is stored here:

solution.controls[:buy]

# ## Introducing risk

# The solution of a single newsvendor problem offers little insight about how
# a decision-maker should act. In particular, they may be averse to bad
# outcomes, such as when they purchase a larger number of newspapers only for
# there to be little demand.

# We can explore how the optimal decision changes with risk by creating a
# function:

function solve_risk_averse_newsvendor(Ω, risk_measure)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = maximum(Ω) * sales_price,
        optimizer = HiGHS.Optimizer,
    ) do subproblem, stage
        @variable(subproblem, inventory >= 0, SDDP.State, initial_value = 0)
        if stage == 1
            @variable(subproblem, buy >= 0)
            @constraint(subproblem, inventory.out == inventory.in + buy)
            @stageobjective(subproblem, -buy_price * buy)
        else
            @variable(subproblem, sell >= 0)
            @constraint(subproblem, sell <= inventory.in)
            SDDP.parameterize(subproblem, Ω) do ω
                return JuMP.set_upper_bound(sell, ω)
            end
            @stageobjective(subproblem, sales_price * sell)
        end
    end
    SDDP.train(
        model;
        risk_measure = risk_measure,
        stopping_rules = [SDDP.SimulationStoppingRule()],
        print_level = 0,
    )
    first_stage_rule = SDDP.DecisionRule(model, node = 1)
    solution = SDDP.evaluate(
        first_stage_rule;
        incoming_state = Dict(:inventory => 0.0),
        controls_to_record = [:buy],
    )
    return solution.controls[:buy]
end

# Now we can see how many units a risk-neutral decision maker would order:

solve_risk_averse_newsvendor(Ω, SDDP.Expectation())

# as well as a decision-maker who cares only about the worst-case outcome:

solve_risk_averse_newsvendor(Ω, SDDP.WorstCase())

# In general, the decision-maker will be somewhere between the two extremes.
# The [`SDDP.Entropic`](@ref) risk measure is a risk measure that has a single
# parameter that lets us explore the space of policies between the two extremes.
# When the parameter is small, the measure acts like [`SDDP.Expectation`](@ref),
# and when it is large, it acts like [`SDDP.WorstCase`](@ref).

# Here is what we get if we solve our problem multiple times for different
# values of the risk aversion parameter ``\gamma``:

Γ = [10^i for i in -3:0.2:3]
buy = [solve_risk_averse_newsvendor(Ω, SDDP.Entropic(γ)) for γ in Γ]
Plots.plot(
    Γ,
    buy;
    xaxis = :log,
    xlabel = "Risk aversion parameter γ",
    ylabel = "First-stage buy decision",
    legend = false,
)

# ## Things to try

# There are a number of things you can try next:

#  * Experiment with different buy and sales prices
#  * Experiment with different distributions of demand
#  * Explore how the optimal policy changes if you use a different risk measure
#  * What happens if you can only buy and sell integer numbers of newspapers?
#    Try this by adding `Int` to the variable definitions:
#    `@variable(subproblem, buy >= 0, Int)`
