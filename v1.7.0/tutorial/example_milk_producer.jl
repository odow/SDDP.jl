#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: the milk producer

# The purpose of this tutorial is to demonstrate how to fit a Markovian policy
# graph to a univariate stochastic process.

# This tutorial uses the following packages:

using SDDP
import HiGHS
import Plots

# ## Background

# A company produces milk for sale on a spot market each month. The quantity of
# milk they produce is uncertain, and so too is the price on the spot market.
# The company can store unsold milk in a stockpile of dried milk powder.

# The spot price is determined by an auction system, and so varies from month to
# month, but demonstrates serial correlation. In each auction, there is
# sufficient demand that the milk producer finds a buyer for all their
# milk, regardless of the quantity they supply. Furthermore, the spot price
# is independent of the milk producer (they are a small player in the market).

# The spot price is highly volatile, and is the result of a process that is out
# of the control of the company. To counteract their price risk, the company
# engages in a forward contracting programme.

# The forward contracting programme is a deal for physical milk four months in
# the future.

# The futures price is the current spot price, plus some forward contango (the
# buyers gain certainty that they will receive the milk in the future).

# In general, the milk company should forward contract (since they reduce
# their price risk), however they also have production risk. Therefore, it may
# be the case that they forward contract a fixed amount, but find that they do
# not produce enough milk to meet the fixed demand. They are then forced to
# buy additional milk on the spot market.

# The goal of the milk company is to choose the extent to which they forward
# contract in order to maximise (risk-adjusted) revenues, whilst managing their
# production risk.

# ## A stochastic process for price

# It is outside the scope of this tutorial, but assume that we have gone away
# and analysed historical data to fit a stochastic process to the sequence of
# monthly auction spot prices.

# One plausible model is a multiplicative auto-regressive model of order one,
# where the white noise term is modeled by a finite distribution of empirical
# residuals. We can simulate this stochastic process as follows:

function simulator()
    residuals = [0.0987, 0.199, 0.303, 0.412, 0.530, 0.661, 0.814, 1.010, 1.290]
    residuals = 0.1 * vcat(-residuals, 0.0, residuals)
    scenario = zeros(12)
    y, μ, α = 4.5, 6.0, 0.05
    for t in 1:12
        y = exp((1 - α) * log(y) + α * log(μ) + rand(residuals))
        scenario[t] = clamp(y, 3.0, 9.0)
    end
    return scenario
end

simulator()

# It may be helpful to visualize a number of simulations of the price process:

plot = Plots.plot(
    [simulator() for _ in 1:500];
    color = "gray",
    opacity = 0.2,
    legend = false,
    xlabel = "Month",
    ylabel = "Price [\$/kg]",
    xlims = (1, 12),
    ylims = (3, 9),
)

# The prices gradually revert to the mean of \$6/kg, and there is high
# volatility.

# We can't incorporate this price process directly into SDDP.jl, but we can fit
# a [`SDDP.MarkovianGraph`](@ref) directly from the simulator:

graph = SDDP.MarkovianGraph(simulator; budget = 30, scenarios = 10_000);
nothing  # hide

# Here `budget` is the number of nodes in the policy graph, and `scenarios` is
# the number of simulations to use when estimating the transition probabilities.

# The graph contains too many nodes to be show, but we can plot it:

for ((t, price), edges) in graph.nodes
    for ((t′, price′), probability) in edges
        Plots.plot!(
            plot,
            [t, t′],
            [price, price′];
            color = "red",
            width = 3 * probability,
        )
    end
end

plot

# That looks okay. Try changing `budget` and `scenarios` to see how different
# Markovian policy graphs can be created.

# ## Model

# Now that we have a Markovian graph, we can build the model. See if you can
# work out how we arrived at this formulation by reading the background
# description. Do all the variables and constraints make sense?

model = SDDP.PolicyGraph(
    graph;
    sense = :Max,
    upper_bound = 1e2,
    optimizer = HiGHS.Optimizer,
) do sp, node
    ## Decompose the node into the month (::Int) and spot price (::Float64)
    t, price = node::Tuple{Int,Float64}
    ## Transactions on the futures market cost 0.01
    c_transaction = 0.01
    ## It costs the company +50% to buy milk on the spot market and deliver to
    ## their customers
    c_buy_premium = 1.5
    ## Buyer is willing to pay +5% for certainty
    c_contango = 1.05
    ## Distribution of production
    Ω_production = range(0.1, 0.2; length = 5)
    c_max_production = 12 * maximum(Ω_production)
    ## x_stock: quantity of milk in stock pile
    @variable(sp, 0 <= x_stock, SDDP.State, initial_value = 0)
    ## x_forward[i]: quantity of milk for delivery in i months
    @variable(sp, 0 <= x_forward[1:4], SDDP.State, initial_value = 0)
    ## u_spot_sell: quantity of milk to sell on spot market
    @variable(sp, 0 <= u_spot_sell <= c_max_production)
    ## u_spot_buy: quantity of milk to buy on spot market
    @variable(sp, 0 <= u_spot_buy <= c_max_production)
    ## u_spot_buy: quantity of milk to sell on futures market
    c_max_futures = t <= 8 ? c_max_production : 0.0
    @variable(sp, 0 <= u_forward_sell <= c_max_futures)
    ## ω_production: production random variable
    @variable(sp, ω_production)
    ## Forward contracting constraints:
    @constraint(sp, [i in 1:3], x_forward[i].out == x_forward[i+1].in)
    @constraint(sp, x_forward[4].out == u_forward_sell)
    ## Stockpile balance constraint
    @constraint(
        sp,
        x_stock.out ==
        x_stock.in + ω_production + u_spot_buy - x_forward[1].in - u_spot_sell
    )
    ## The random variables. `price` comes from the Markov node
    ##
    ## !!! warning
    ##     The elements in Ω MUST be a tuple with 1 or 2 values, where the first
    ##     value is `price` and the second value is the random variable for the
    ##     current node. If the node is deterministic, use Ω = [(price,)].
    Ω = [(price, p) for p in Ω_production]
    SDDP.parameterize(sp, Ω) do ω
        ## Fix the ω_production variable
        fix(ω_production, ω[2])
        @stageobjective(
            sp,
            ## Sales on spot market
            ω[1] * (u_spot_sell - c_buy_premium * u_spot_buy) +
            ## Sales on futures smarket
            (ω[1] * c_contango - c_transaction) * u_forward_sell
        )
        return
    end
    return
end

# ## Training a policy

# Now we have a model, we train a policy. The [`SDDP.SimulatorSamplingScheme`](@ref)
# is used in the forward pass. It generates an out-of-sample sequence of prices
# using `simulator` and traverses the closest sequence of nodes in the policy
# graph. When calling [`SDDP.parameterize`](@ref) for each subproblem, it uses
# the new out-of-sample price instead of the price associated with the Markov
# node.

SDDP.train(
    model;
    time_limit = 20,
    risk_measure = 0.5 * SDDP.Expectation() + 0.5 * SDDP.AVaR(0.25),
    sampling_scheme = SDDP.SimulatorSamplingScheme(simulator),
)

# !!! warning
#     We're intentionally terminating the training early so that the
#     documentation doesn't take too long to build. If you run this example
#     locally, increase the time limit.

# ## Simulating the policy

# When simulating the policy, we can also use the
# [`SDDP.SimulatorSamplingScheme`](@ref).

simulations = SDDP.simulate(
    model,
    200,
    Symbol[:x_stock, :u_forward_sell, :u_spot_sell, :u_spot_buy];
    sampling_scheme = SDDP.SimulatorSamplingScheme(simulator),
);
nothing  # hide

# To show how the sampling scheme uses the new out-of-sample price instead of
# the price associated with the Markov node, compare the index of the Markov
# state visited in stage 12 of the first simulation:

simulations[1][12][:node_index]

# to the realization of the noise `(price, ω)` passed to [`SDDP.parameterize`](@ref):

simulations[1][12][:noise_term]

# ## Visualizing the policy

# Finally, we can plot the policy to gain insight (although note that we
# terminated the training early, so we should run the re-train the policy for
# more iterations before making too many judgements).

plot = Plots.plot(
    SDDP.publication_plot(simulations; title = "x_stock.out") do data
        return data[:x_stock].out
    end,
    SDDP.publication_plot(simulations; title = "u_forward_sell") do data
        return data[:u_forward_sell]
    end,
    SDDP.publication_plot(simulations; title = "u_spot_buy") do data
        return data[:u_spot_buy]
    end,
    SDDP.publication_plot(simulations; title = "u_spot_sell") do data
        return data[:u_spot_sell]
    end;
    layout = (2, 2),
)

# ## Next steps

# * Train the policy for longer. What do you observe?
# * Try creating different Markovian graphs. What happens if you add more nodes?
# * Try different risk measures
