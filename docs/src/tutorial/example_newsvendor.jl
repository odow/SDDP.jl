#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: two-stage newsvendor

# The purpose of this tutorial is to demonstrate how to model and solve a
# two-stage stochastic program.

# It is based on the [Two stage stochastic programs](https://jump.dev/JuMP.jl/dev/tutorials/applications/two_stage_stochastic/)
# tutorial in JuMP.

# This tutorial uses the following packages

using JuMP
using SDDP
import Distributions
import HiGHS
import Plots
import StatsPlots
import Statistics

# ## Background

# The data for this problem is:

D = Distributions.TriangularDist(150.0, 250.0, 200.0)
N = 100
d = sort!(rand(D, N));
Ω = 1:N
P = fill(1 / N, N);
StatsPlots.histogram(d; bins = 20, label = "", xlabel = "Demand")

# ## The L-Shaped method

# The L-Shaped method is a way of solving two-stage stochastic programs by
# Benders' decomposition. It takes the problem:

# ```math
# \begin{aligned}
# \max\limits_{x,y_\omega} \;\; & -2x + \mathbb{E}_\omega[5y_\omega - 0.1(x - y_\omega)] \\
#   & y_\omega \le x              & \quad \forall \omega \in \Omega \\
#   & 0 \le y_\omega \le d_\omega & \quad \forall \omega \in \Omega \\
#   & x \ge 0.
# \end{aligned}
# ```

# and decomposes it into a second-stage problem:

# ```math
# \begin{aligned}
# V_2(\bar{x}, d_\omega) = \max\limits_{x,x^\prime,y_\omega} \;\; & 5y_\omega - x^\prime \\
#   & y_\omega \le x \\
#   & x^\prime = x - y_\omega \\
#   & 0 \le y_\omega \le d_\omega \\
#   & x = \bar{x} & [\lambda]
# \end{aligned}
# ```

# and a first-stage problem:

# ```math
# \begin{aligned}
# V^K = \max\limits_{x,\theta} \;\; & -2x + \theta \\
#   & \theta \le \mathbb{E}_\omega[V_2(x^k, \omega) + \lambda^k(x - x^k)] & \quad k = 1,\ldots,K\\
#   & x \ge 0,
# \end{aligned}
# ```

# Here's a function to compute the second-stage problem;

function solve_second_stage(x̅, d_ω)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x_in)
    @variable(model, x_out >= 0)
    fix(x_in, x̅)
    @variable(model, 0 <= u_sell <= d_ω)
    @constraint(model, x_out == x_in - u_sell)
    @constraint(model, u_sell <= x_in)
    @objective(model, Max, 5 * u_sell - 0.1 * x_out)
    optimize!(model)
    return (
        V = objective_value(model),
        λ = reduced_cost(x_in),
        x = value(x_out),
        u = value(u_sell),
    )
end

solve_second_stage(200, 170)

# Here's the first-stage subproblem:

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x_in == 0)
@variable(model, x_out >= 0)
@variable(model, u_make >= 0)
@constraint(model, x_out == x_in + u_make)
@variable(model, θ <= 10_000)
@objective(model, Max, -2 * u_make + θ)

# Importantly, to ensure we have a bounded solution, we need to add an upper
# bound to the variable `θ`.

kIterationLimit = 100
for k in 1:kIterationLimit
    println("Solving iteration k = $k")
    optimize!(model)
    xᵏ = value(x_out)
    println("  xᵏ = $xᵏ")
    ub = objective_value(model)
    println("  V̅ = $ub")
    ret = [solve_second_stage(xᵏ, d[ω]) for ω in Ω]
    lb = value(-2 * u_make) + sum(p * r.V for (p, r) in zip(P, ret))
    println("  V̲ = $lb")
    if ub - lb < 1e-6
        println("Terminating with near-optimal solution")
        break
    end
    c = @constraint(
        model,
        θ <= sum(p * (r.V + r.λ * (x_out - xᵏ)) for (p, r) in zip(P, ret)),
    )
    println("  Added cut: $c")
end

# ## Policy Graph

# Now we can formulate and train a policy for the two-stage newsvendor problem.

model = SDDP.LinearPolicyGraph(;
    stages = 2,
    sense = :Max,
    upper_bound = 5 * maximum(d),  # The `M` in θ <= M
    optimizer = HiGHS.Optimizer,
) do subproblem::JuMP.Model, stage::Int
    @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
    if stage == 1
        @variable(subproblem, u_make >= 0)
        @constraint(subproblem, x.out == x.in + u_make)
        @stageobjective(subproblem, -2 * u_make)
    else
        @variable(subproblem, u_sell >= 0)
        @constraint(subproblem, u_sell <= x.in)
        @constraint(subproblem, x.out == x.in - u_sell)
        SDDP.parameterize(subproblem, d, P) do ω
            set_upper_bound(u_sell, ω)
            return
        end
        @stageobjective(subproblem, 5 * u_sell - 0.1 * x.out)
    end
    return
end

SDDP.train(model)

# One way to query the optimal policy is with [`SDDP.DecisionRule`](@ref):

first_stage_rule = SDDP.DecisionRule(model; node = 1)

#-

solution_1 = SDDP.evaluate(first_stage_rule; incoming_state = Dict(:x => 0.0))

# Here's the second stage:

second_stage_rule = SDDP.DecisionRule(model; node = 2)
solution = SDDP.evaluate(
    second_stage_rule;
    incoming_state = Dict(:x => solution_1.outgoing_state[:x]),
    noise = 170.0,  # A value of d[ω], can be out-of-sample.
    controls_to_record = [:u_sell],
)

# ## Simulation

# Querying the decision rules is tedious. It's often more useful to simulate the
# policy:

simulations = SDDP.simulate(
    model,
    10,  #= number of replications =#
    [:x, :u_sell, :u_make];  #= variables to record =#
    skip_undefined_variables = true,
);

# `simulations` is a vector with 10 elements

length(simulations)

# and each element is a vector with two elements (one for each stage)

length(simulations[1])

# The first stage contains:

simulations[1][1]

# The second stage contains:

simulations[1][2]

# We can compute aggregated statistics across the simulations:

objectives = map(simulations) do simulation
    return sum(data[:stage_objective] for data in simulation)
end
μ, t = SDDP.confidence_interval(objectives)
println("Simulation ci : $μ ± $t")

# ## Risk aversion revisited

# SDDP.jl contains a number of risk measures. One example is:

0.5 * SDDP.Expectation() + 0.5 * SDDP.WorstCase()

# You can construct a risk-averse policy by passing a risk measure to the
# `risk_measure` keyword argument of [`SDDP.train`](@ref).

# We can explore how the optimal decision changes with risk by creating a
# function:

function solve_newsvendor(risk_measure::SDDP.AbstractRiskMeasure)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        sense = :Max,
        upper_bound = 5 * maximum(d),
        optimizer = HiGHS.Optimizer,
    ) do subproblem, node
        @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
        if node == 1
            @stageobjective(subproblem, -2 * x.out)
        else
            @variable(subproblem, u_sell >= 0)
            @constraint(subproblem, u_sell <= x.in)
            @constraint(subproblem, x.out == x.in - u_sell)
            SDDP.parameterize(subproblem, d, P) do ω
                set_upper_bound(u_sell, ω)
                return
            end
            @stageobjective(subproblem, 5 * u_sell - 0.1 * x.out)
        end
        return
    end
    SDDP.train(model; risk_measure = risk_measure, print_level = 0)
    first_stage_rule = SDDP.DecisionRule(model; node = 1)
    solution = SDDP.evaluate(first_stage_rule; incoming_state = Dict(:x => 0.0))
    return solution.outgoing_state[:x]
end

# Now we can see how many units a decision maker would order using `CVaR`:

solve_newsvendor(SDDP.CVaR(0.4))

# as well as a decision-maker who cares only about the worst-case outcome:

solve_newsvendor(SDDP.WorstCase())

# In general, the decision-maker will be somewhere between the two extremes.
# The [`SDDP.Entropic`](@ref) risk measure is a risk measure that has a single
# parameter that lets us explore the space of policies between the two extremes.
# When the parameter is small, the measure acts like [`SDDP.Expectation`](@ref),
# and when it is large, it acts like [`SDDP.WorstCase`](@ref).

# Here is what we get if we solve our problem multiple times for different
# values of the risk aversion parameter ``\gamma``:

Γ = [10^i for i in -4:0.5:1]
buy = [solve_newsvendor(SDDP.Entropic(γ)) for γ in Γ]
Plots.plot(
    Γ,
    buy;
    xaxis = :log,
    xlabel = "Risk aversion parameter γ",
    ylabel = "Number of pies to make",
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
#  * What happens if you use a different upper bound? Try an invalid one like
#    `-100`, and a very large one like `1e12`.
