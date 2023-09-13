#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: Two-stage Newsvendor

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

# ## Policy Graph

# Now we can formulate and train a policy for the two-stage newsvendor problem.

# First, we need to construct the graph:

graph = SDDP.LinearGraph(2)

# Then, we need to write a function which builds a JuMP model for each node in
# the graph:

function build_subproblem(subproblem::JuMP.Model, stage::Int)
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

# Then, we can combine the graph and the subproblem builder into a policy graph:

model = SDDP.PolicyGraph(
    build_subproblem,
    graph;
    sense = :Max,
    upper_bound = 5 * maximum(d), #= some large upper bound =#
    optimizer = HiGHS.Optimizer,
)

# Use [`SDDP.train`](@ref) to construct the policy:

SDDP.train(model)

# To check the first-stage buy decision, we need to obtain a decision rule for
# the first-stage node `1`:

first_stage_rule = SDDP.DecisionRule(model, node = 1)

# Then we can evaluate it, passing in a starting point for the incoming state:

solution = SDDP.evaluate(first_stage_rule; incoming_state = Dict(:x => 0.0))

# The optimal value of the state variable is stored here:

solution.outgoing_state[:x]

# We can simplify the model construction by using [`SDDP.LinearPolicyGraph`](@ref):

model = SDDP.LinearPolicyGraph(
    build_subproblem;
    stages = 2,
    sense = :Max,
    upper_bound = 5 * maximum(d),
    optimizer = HiGHS.Optimizer,
)

# and we can use Julia's `do` syntax to avoid writing a separate function:

model = SDDP.LinearPolicyGraph(;
    stages = 2,
    sense = :Max,
    upper_bound = 5 * maximum(d),
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

# ## Simulation

# Querying the decision rules is tedious. It's often more useful to simulate the
# policy:

SDDP.train(model)

#-

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

# ## Other graphs

# The two-stage newsvendor problem can be extended to other graphs. For example,
# here is a three-stage graph:

graph = SDDP.Graph(:root)
SDDP.add_node(graph, :make_only)
SDDP.add_node(graph, :sell_only)
SDDP.add_node(graph, :make_and_sell)
SDDP.add_edge(graph, :root => :make_only, 1.0)
SDDP.add_edge(graph, :make_only => :make_and_sell, 1.0)
SDDP.add_edge(graph, :make_and_sell => :sell_only, 1.0)
graph

# with the model formulation:

model = SDDP.PolicyGraph(
    graph;
    sense = :Max,
    upper_bound = 2 * 5 * maximum(d),
    optimizer = HiGHS.Optimizer,
) do subproblem::JuMP.Model, node::Symbol
    @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
    if node == :make_only
        @variable(subproblem, u_make >= 0)
        @constraint(subproblem, x.out == x.in + u_make)
        @stageobjective(subproblem, -2 * u_make)
    elseif node == :sell_only
        @variable(subproblem, u_sell >= 0)
        @constraint(subproblem, u_sell <= x.in)
        @constraint(subproblem, x.out == x.in - u_sell)
        SDDP.parameterize(subproblem, d, P) do ω
            set_upper_bound(u_sell, ω)
            return
        end
        @stageobjective(subproblem, 5 * u_sell - 0.1 * x.out)
    else
        @assert node == :make_and_sell
        @variable(subproblem, u_sell >= 0)
        @variable(subproblem, u_make >= 0)
        @constraint(subproblem, u_sell <= x.in)
        @constraint(subproblem, x.out == x.in - u_sell + u_make)
        SDDP.parameterize(subproblem, d, P) do ω
            set_upper_bound(u_sell, ω)
            return
        end
        @stageobjective(subproblem, 5 * u_sell - 0.1 * x.out - 2 * u_make)
    end
    return
end

SDDP.train(model)
simulations = SDDP.simulate(
    model,
    10,  #= number of replications =#
    [:x, :u_sell, :u_make];  #= variables to record =#
    skip_undefined_variables = true,
);

#-

simulations[1][1]

#-

simulations[1][2]

#-

simulations[1][3]

# Try changing the graph. What happens if you add a cycle?
