#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Two-stage stochastic programs

# The purpose of this tutorial is to demonstrate how to model and solve a
# two-stage stochastic program.

# This tutorial uses the following packages

using JuMP
using SDDP
import Distributions
import HiGHS
import Plots
import StatsPlots
import Statistics

# ## Background

# During the week, you are a busy practitioner of Operations Research. To escape
# the drudgery of mathematics, you decide to invest in a food truck that you
# park at a popular park on a Saturday afternoon. Your delicacy of choice is
# a creamy mushroom pie with puff pastry. After a few weeks, it quickly becomes
# apparent that operating a food business is not so easy.

# The pies must be prepared on a Saturday morning, _before_ you arrive at the
# location and can gauge the level of demand. If you bake too many, the unsold
# pies at the end of the day must be discarded and you have wasted time and
# money on their production. But if you bake too few, then there may be
# un-served customers and you could have made more money by baking more pies.

# After a few weeks of poor decision making, you decide to put your knowledge of
# Operations Research to good use, starting with some data collection.

# Each pie costs you $2 to make, and you sell them at $5 each. Disposal of an
# unsold pie costs $0.10. Based on three weeks of data collected, in which you
# made 200 pies each week, you sold 150, 190, and 200 pies. Thus, as a guess,
# you assume a triangular distribution of demand with a minimum of 150, a median
# of 200, and a maximum of 250.

# We can model this problem by a two-stage stochastic program. In the first
# stage, we decide a quantity of pies to make ``x``. We make this decision
# before we observe the demand ``d_\omega``. In the second stage, we sell
# ``y_\omega`` pies, and incur any costs for unsold pies.

# We can formulate this problem as follows:
# ```math
# \begin{aligned}
# \min 2x + \mathbb{E}_\omega[-5y_\omega + 0.1(x - y_\omega)] \\
# y_\omega \le x \quad \forall \omega \in \Omega \\
# 0 \le y_\omega \le d_\omega \quad \forall \omega \in \Omega \\
# x \ge 0.
# \end{aligned}
# ```
# Note that we are minimizingn cost, so the profit from selling a pie is a
# negative cost.

# ## Sample Average approximation

# If the distribution of demand is continuous, then our problem has an infinite
# number of variables and constraints. To form a computationally tractable
# problem, we instead use a finite set of samples drawn from the distribution.
# This is called sample average approximation (SAA).

D = Distributions.TriangularDist(150.0, 250.0, 200.0)
N = 100
d = sort!(rand(D, N));
Ω = 1:N
P = fill(1 / N, N);
StatsPlots.histogram(d; bins = 20)

# ## JuMP model

# The implementation of our two-stage stochastic program in JuMP is:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@variable(model, 0 <= y[ω in Ω] <= d[ω])
@constraint(model, [ω in Ω], y[ω] <= x)
@expression(model, z[ω in Ω], -5y[ω] + 0.1 * (x - y[ω]))
@objective(model, Min, 2x + sum(P[ω] * z[ω] for ω in Ω))
optimize!(model)
solution_summary(model)

# The optimal number of pies to make is:

value(x)

# The distribution of total profit is:

total_cost = [2 * value(x) + value(z[ω]) for ω in Ω]

bins = extrema(total_cost)
plot = StatsPlots.histogram(
    -total_cost;
    bins = 20,
    label = "",
    xlabel = "Profit [\$]",
    ylabel = "Number of outcomes",
)
μ = Statistics.mean(-total_cost)
Plots.vline!(
    plot,
    [μ];
    label = "Expected profit (\$$(round(Int, μ)))",
    linewidth = 3,
)

# ## Exercises

#  * Try solving this problem for different numbers of samples and different
#    distributions.
#  * Refactor the example to avoid hard-coding the costs. What happens to the
#    solution if the cost of disposing unsold pies increases?

# ## Risk measures

# A risk measure is a function which maps a random variable to a real number.
# Common risk measures include the mean (expectation), median, mode, and
# maximum. We need a risk measure to convert the distribution of second stage
# costs into a single number that can be optimized.

# Our model currently uses the expectation risk measure, but others are possible
# too. One popular risk measure is the conditional value at risk (CVaR).

# CVaR has a parameter $\gamma$. Thhe risk measure computes the expectation of
# the worst $\gamma$ fraction of outcomes.

# The definition of CVaR is:
# ```math
# CVaR_{\gamma}[Z] = \min\limits_{\xi} \xi + \frac{1}{\gamma}\mathbb{E}_\omega[(Z - \xi)_+]
# ```
# which can be formulated as the linear program:
# ```math
# \begin{aligned}
# CVaR_{\gamma}[Z] = \min \xi + \frac{1}{\gamma}\sum P_\omega z_\omega\\
#  & z_\omega \ge Z_\omega - \xi \\
#  & z_\omega \ge 0 \forall \omega
# \end{aligned}
# ```

function CVaR(Z::Vector{Float64}, P::Vector{Float64}; γ::Float64)
    @assert 0 < γ <= 1
    N = length(Z)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, ξ)
    @variable(model, z[1:N] >= 0)
    @constraint(model, [i in 1:N], z[i] >= Z[i] - ξ)
    @objective(model, Min, ξ + 1 / γ * sum(P[i] * z[i] for i in 1:N))
    optimize!(model)
    return objective_value(model)
end

# When `γ` is `1.0`, we compute the mean of the cost:

cvar_10 = CVaR(total_cost, P; γ = 1.0)

#-

Statistics.mean(total_cost)

# As `γ` approaches `0.0`, we compute the worst-case (maximum) cost:

cvar_00 = CVaR(total_cost, P; γ = 0.0001)

#-

maximum(total_cost)

# By varying `γ` between `0` and `1` we can compute some trade-off of these two
# extremes:

cvar_05 = CVaR(total_cost, P; γ = 0.5)

# Let's plot these outcomes on our distribution:

plot = StatsPlots.histogram(
    -total_cost;
    bins = 20,
    label = "",
    xlabel = "Profit [\$]",
    ylabel = "Number of outcomes",
)
Plots.vline!(
    plot,
    -[cvar_10 cvar_05 cvar_00];
    label = ["γ = 1.0" "γ = 0.5" "γ = 0.0"],
    linewidth = 3,
)

# ## Risk averse sample average approximation

# Because CVaR can be formulated as a minimization linear program, we can form a
# risk averse sample average approximation model by combining the two
# formulations:

γ = 0.4
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x >= 0)
@variable(model, 0 <= y[ω in Ω] <= d[ω])
@constraint(model, [ω in Ω], y[ω] <= x)
@expression(model, Z[ω in Ω], -5 * y[ω] + 0.1(x - y[ω]))
@variable(model, ξ)
@variable(model, z[ω in Ω] >= 0)
@constraint(model, [ω in Ω], z[ω] >= Z[ω] - ξ)
@objective(model, Min, 2x + ξ + 1 / γ * sum(P[ω] * z[ω] for ω in Ω))
optimize!(model)

# When ``\gamma = 0.4``, the optimal number of pies to bake is:

value(x)

# The distribution of total profit is:

risk_averse_total_profit = [-value(2x + Z[ω]) for ω in Ω]
plot = StatsPlots.histogram(
    -total_cost;
    label = "Expectation",
    bins = 300:25:650,
)
StatsPlots.histogram!(
    risk_averse_total_profit;
    label = "CV@R",
    bins = 300:25:650,
    alpha = 0.5,
)

# ## Policy Graph

# Now we can formulate and train a policy for the two-stage newsvendor problem.

# First, we need to construct the graph:

graph = SDDP.LinearGraph(2)

# Then, we need to write a function which builds a JuMP model for eachh node in
# the graph:

function build_subproblem(subproblem::JuMP.Model, stage::Int)
    @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
    if stage == 1
        @stageobjective(subproblem, 2 * x.out)
    else
        @variable(subproblem, y >= 0)
        @constraint(subproblem, y <= x.in)
        SDDP.parameterize(subproblem, d, P) do ω
            set_upper_bound(y, ω)
            return
        end
        @stageobjective(subproblem, -5 * y + 0.1 * (x.in - y))
    end
    return
end

# Then, we can combine the graph and the subproblem builder into a policy graph:

model = SDDP.PolicyGraph(
    build_subproblem,
    graph;
    sense = :Min,
    lower_bound = -5 * maximum(d),
    optimizer = HiGHS.Optimizer,
)

# Use [`SDDP.train`](@ref) to construct the policy:

SDDP.train(model)

# To check the first-stage buy decision, we need to obtain a decision rule for
# the first-stage node `1`:

first_stage_rule = SDDP.DecisionRule(model, node = 1)

# Then we can evaluate it, passing in a starting point for the incoming state:

solution = SDDP.evaluate(
    first_stage_rule;
    incoming_state = Dict(:x => 0.0),
)

# The optimal value of the `buy` variable is stored here:

solution.outgoing_state[:x]

# We can simplify the model construction by using [`SDDP.LinearPolicyGraph`](@ref):

model = SDDP.LinearPolicyGraph(
    build_subproblem;
    stages = 2,
    sense = :Min,
    lower_bound = -5 * maximum(d),
    optimizer = HiGHS.Optimizer,
)

# and we can use Julia's `do` syntax to avoid writing a separate functionn:

model = SDDP.LinearPolicyGraph(;
    stages = 2,
    sense = :Min,
    lower_bound = -5 * maximum(d),
    optimizer = HiGHS.Optimizer,
) do subproblem::JuMP.Model, stage::Int
    @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
    if stage == 1
        @stageobjective(subproblem, 2 * x.out)
    else
        @variable(subproblem, y >= 0)
        @constraint(subproblem, y <= x.in)
        SDDP.parameterize(subproblem, d, P) do ω
            set_upper_bound(y, ω)
            return
        end
        @stageobjective(subproblem, -5 * y + 0.1 * (x.in - y))
    end
    return
end

# ## Introducing risk

# The solution of a single newsvendor problem offers little insight about how
# a decision-maker should act. In particular, they may be averse to bad
# outcomes, such as when they purchase a larger number of newspapers only for
# there to be little demand.

# We can explore how the optimal decision changes with risk by creating a
# function:

function solve_risk_averse_newsvendor(risk_measure)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Min,
        lower_bound = -5 * maximum(d),
        optimizer = HiGHS.Optimizer,
    ) do subproblem, stage
        @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
        if stage == 1
            @stageobjective(subproblem, 2 * x.out)
        else
            @variable(subproblem, y >= 0)
            @constraint(subproblem, y <= x.in)
            SDDP.parameterize(subproblem, d, P) do ω
                set_upper_bound(y, ω)
                return
            end
            @stageobjective(subproblem, -5 * y + 0.1 * (x.in - y))
        end
        return
    end
    SDDP.train(model; risk_measure = risk_measure, print_level = 0)
    first_stage_rule = SDDP.DecisionRule(model; node = 1)
    solution = SDDP.evaluate(
        first_stage_rule;
        incoming_state = Dict(:x => 0.0),
    )
    return solution.outgoing_state[:x]
end

# Now we can see how many units a risk-neutral decision maker would order:

solve_risk_averse_newsvendor(SDDP.AVaR(0.4))

# as well as a decision-maker who cares only about the worst-case outcome:

solve_risk_averse_newsvendor(SDDP.WorstCase())

# In general, the decision-maker will be somewhere between the two extremes.
# The [`SDDP.Entropic`](@ref) risk measure is a risk measure that has a single
# parameter that lets us explore the space of policies between the two extremes.
# When the parameter is small, the measure acts like [`SDDP.Expectation`](@ref),
# and when it is large, it acts like [`SDDP.WorstCase`](@ref).

# Here is what we get if we solve our problem multiple times for different
# values of the risk aversion parameter ``\gamma``:

Γ = [10^i for i in -5:0.5:2]
buy = [solve_risk_averse_newsvendor(SDDP.Entropic(γ)) for γ in Γ]
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
