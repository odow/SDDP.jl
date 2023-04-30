#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Example: the perishable widget producer

# A company produces perishable, divisible widgets (i.e. a continuous product,
# rather than discrete units) for sale on a spot market each month. The quantity
# of widgets they produce is uncertain and comes from some finite distribution.

# The company can store the widgets, but, over time, the value of the widgets
# decreases. Eventually the widgets perish and are worthless.

# The spot price is determined by an auction system, and so varies from month to
# month, but demonstrates serial correlation. In each auction, there is
# sufficient demand that the widget producer finds a buyer for all their
# widgets, regardless of the quantity they supply. Furthermore, the spot price
# is independent of the widget producer (they are a small player in the market).

# The spot price is highly volatile, and is the result of a process that is out
# of the control of the company. To counteract their price risk, the company
# engages in a forward contracting programme.

# The forward contracting programme is a deal for physical widgets at a future
# date in time. The company can engage in `x_forward` for sales up to four
# months in the future.

# The futures price is the current spot price, plus some forward contango (the
# buyers gain certainty that they will receive the widgets in the future).

# In general, the widget company should forward contract (since they reduce
# their price risk), however they also have production risk. Therefore, it may
# be the case that they forward contract a fixed amount, but find that they do
# not produce enough widgets to meet the fixed demand. They are then forced to
# buy additional widgets on the spot market.

# The goal of the widget company is to choose the extent to which they forward
# contract in order to maximise (risk-adjusted) revenues, whilst managing their
# production risk.

using SDDP
import HiGHS

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

model = SDDP.PolicyGraph(
    SDDP.MarkovianGraph(simulator; budget = 60, scenarios = 100_000);
    sense = :Max,
    upper_bound = 1e4,
    optimizer = HiGHS.Optimizer,
) do sp, node
    t, markov_state = node
    c_contango = [1.0, 1.025, 1.05, 1.075]
    c_transaction, c_perish_factor, c_buy_premium = 0.01, 0.95, 1.5
    F, P = length(c_contango), 5
    @variable(sp, 0 <= x_forward[1:F], SDDP.State, (initial_value = 0))
    @variable(sp, 0 <= x_stock[p = 1:P], SDDP.State, (initial_value = 0))
    @variable(sp, 0 <= u_spot_sell[1:P])
    @variable(sp, 0 <= u_spot_buy)
    @variable(sp, 0 <= u_forward_deliver[1:P])
    @variable(sp, 0 <= u_forward_sell[1:F] <= 10)
    @variable(sp, 0 <= u_production)
    for f in 1:F
        if t + f >= 12
            fix(u_forward_sell[f], 0.0; force = true)
        end
    end
    ## Contract balance
    @constraint(
        sp,
        [i = 1:F-1],
        x_forward[i].out == x_forward[i+1].in + u_forward_sell[i],
    )
    @constraint(sp, x_forward[F].out == u_forward_sell[F])
    @constraint(sp, x_forward[1].in == sum(u_forward_deliver))
    ## Inventory balance
    @constraint(
        sp,
        x_stock[1].out ==
        u_production - u_forward_deliver[1] - u_spot_sell[1] + u_spot_buy
    )
    @constraint(
        sp,
        [i = 2:P],
        x_stock[i].out ==
        c_perish_factor * x_stock[i-1].in - u_forward_deliver[i] -
        u_spot_sell[i]
    )
    Ω = [(markov_state, (production = p,)) for p in range(0.1, 0.2; length = 5)]
    SDDP.parameterize(sp, Ω) do (markov_state, ω)
        set_upper_bound(u_production, ω.production)
        @stageobjective(
            sp,
            markov_state * sum(u_spot_sell) +
            markov_state * sum(c_contango[f] * u_forward_sell[f] for f in 1:F) -
            markov_state * c_buy_premium * u_spot_buy -
            c_transaction * sum(u_forward_sell)
        )
    end
end

status = SDDP.train(
    model;
    iteration_limit = 300,
    time_limit = 20,
    risk_measure = SDDP.EAVaR(; lambda = 0.5, beta = 0.25),
    sampling_scheme = SDDP.SimulatorSamplingScheme(simulator),
)

simulations = SDDP.simulate(
    model,
    200,
    Symbol[
        :x_forward,
        :x_stock,
        :u_spot_sell,
        :u_forward_sell,
        :u_forward_deliver,
        :u_spot_buy,
        :u_production,
    ];
    sampling_scheme = SDDP.SimulatorSamplingScheme(simulator),
);

# plt = SDDP.SpaghettiPlot(simulations)
# SDDP.add_spaghetti(plt; title = "u_spot_buy") do data
#     return data[:u_spot_buy]
# end
# SDDP.add_spaghetti(plt; title = "u_production") do data
#     return data[:u_production]
# end
# SDDP.add_spaghetti(plt; title = "sum(x_stock.out)") do data
#     return sum(y.out for y in data[:x_stock])
# end
# SDDP.add_spaghetti(plt; title = "sum(x_forward.out)") do data
#     return sum(y.out for y in data[:x_forward])
# end

# SDDP.add_spaghetti(plt; title = "sum(u_forward_sell)") do data
#     return sum(data[:u_forward_sell])
# end
# SDDP.add_spaghetti(plt; title = "sum(u_spot_sell)") do data
#     return sum(data[:u_spot_sell])
# end
# SDDP.plot(plt)

# for t in 1:5
#     SDDP.add_spaghetti(plt; title = "x_foward_$t") do data
#         return data[:x_forward][t].out
#     end
#     SDDP.add_spaghetti(plt; title = "x_stock_$t") do data
#         return data[:x_stock][t].out
#     end
#     SDDP.add_spaghetti(plt; title = "u_forward_sell_$t") do data
#         return data[:u_forward_sell][t]
#     end
#     SDDP.add_spaghetti(plt; title = "u_spot_sell_$t") do data
#         return data[:u_spot_sell][t]
#     end
# end
# SDDP.plot(plt)
