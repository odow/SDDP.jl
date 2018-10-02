#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    Example: simple contracting.

    In this example, we consider the problem faced by a producer of widgets over
    a single year, which we break into 24 periods, each representing half a
    month. To complicate matters, the producer is unaware of the price they will
    receive for their widgets until the final time period. Instead, the producer
    can observe a forecast for the final price that becomes more accurate as the
    final time period approaches. Moreover, the observed forecast is unbiased,
    and there exists a forward-contract for widgets that trades at the forecast
    price in each time period, and settles against the final price.

    We assume that we can model the forecast of the final price over time by the
    log-normal, auto-regressive process with lag one:
        log(price[t]) = log(1.01) + log(price[t-1]) + σ²(t) × ϕ,
    where ϕ is a stagewise independent random variable and σ²(t) is a constant
    in t that is 1 in stage 1, and linearly declines to 0 in the final stage T.

    We introduce a state variable `production[t]` which counts the number of
    widgets produced at the end of stage `t`.
        production[t] = production[t-1] + ε,
    where ε is a stagewise independent random variable that takes the value 0,
    1, or 2 with equal probability.

    In addition, the producer can sell contracts in the futures market (this
    incurs a transaction cost for each contract). We introduce a state variable
    `contracts[t]` which counts the number of contracts sold at the end of
    stage `t`.
        contracts[t] = contracts[t-1] + sell[t].

    In each stage t={1,2,...,T-1}, the producer earns
        sell[t] × (price[t] - δ),
    where δ is the transaction cost. In the final time period `T`, the producer
    sells the total units produced at the price `price[T]`, and settles their
    contracting obligations. Therefore, they earn
        price[T] × (production[T] - contracts[T]).
    Note that no contracts are sold in the final stage.

    If the producer is risk-neutral, they would not engage in contracting due to
    the transaction cost. Therefore, we model their risk-aversion with a nested
    convex-combination of Expectation and AV@R:
        F[X] = 0.8×E[X] + 0.2×AV@R(0.5)[X].
=#

using SDDP, JuMP, Clp, Base.Test

function contracting_example(DISCRETIZATION = 1)
    srand(10)
    T             = 24
    MIN_PRICE     = 3.0
    INITIAL_PRICE = 6.0
    MAX_PRICE     = 9.0
    NOISES        = DiscreteDistribution([-0.1290, -0.1010, -0.0814, -0.0661, -0.0530,
        -0.0412, -0.0303, -0.0199, -0.00987, 0.0, 0.00987, 0.0199, 0.0303, 0.0412,
        0.0530, 0.0661, 0.0814, 0.1010, 0.1290])

    function pricedynamics(price, noise, stage, markov)
         # Decreasing variance in changes in price over time
        σ² = linspace(1, 0, 24)
        next_price = 1.01 * exp(log(price) + σ²[stage]*noise)
        min(MAX_PRICE,max(MIN_PRICE, next_price))
    end

    value_function = if DISCRETIZATION == 1
        (t,i) -> DynamicPriceInterpolation(
            dynamics       = (p,w) -> pricedynamics(p,w,t,i),
            initial_price  = INITIAL_PRICE,
            min_price      = MIN_PRICE,
            max_price      = MAX_PRICE,
            noise          = NOISES,
            cut_oracle = SDDP.NanniciniOracle(typeof(INITIAL_PRICE), 20),
            lipschitz_constant = 75.0
        )
    else
        (t,i) -> StaticPriceInterpolation(
            dynamics       = (p,w) -> pricedynamics(p,w,t,i),
            initial_price  = INITIAL_PRICE,
            rib_locations  =  collect(linspace(MIN_PRICE, MAX_PRICE, DISCRETIZATION)),
            noise          = NOISES,
                cut_oracle = LevelOneCutOracle(),
        )
    end

    m = SDDPModel(
        sense             = :Max,
        stages            = T,
        objective_bound   = 200.0,
        solver            = ClpSolver(),
        risk_measure      = EAVaR(lambda=0.8, beta=0.5),
        value_function    = value_function
                                            ) do sp, t
        @states(sp, begin
             contracts >= 0, contracts0 == 0
            production >= 0, production0 == 0
        end)
        @variable(sp, sell >= 0)
        @constraints(sp, begin
            contracts  == contracts0 + sell
            # place an upper limit to prevent attempted arbitrage when
            # value function approximation is poor
            contracts <= 2T
        end)
        @rhsnoise(sp, ε = [0, 1, 2], production == production0 + ε)
        if t == T
            @stageobjective(sp,
                price -> (production - contracts) * price
            )
            @constraint(sp, sell == 0)
        else
            δ = 0.01 # transaction cost
            @stageobjective(sp,
                price -> sell * (price - δ)
            )
        end
    end
    return m
end

# dynamic interpolation
m = contracting_example()
srand(123)
SDDP.solve(m, iteration_limit = 50, cut_selection_frequency=10, print_level=0)
@test SDDP.getbound(m) <= 175.0

# historical simulation
results = simulate(m, [:production, :contracts],
    noises=fill(2, 24), pricenoises=fill(10, 24)
)
@test isapprox(results[:objective], 181, atol=1)
results = simulate(m, 500)
@test isapprox(mean(r[:objective] for r in results), 171.0, atol=1.0)

# 3 fixed ribs
m3 = contracting_example(3)
srand(123)
SDDP.solve(m3, iteration_limit = 10, cut_selection_frequency=5, print_level=0)
@test SDDP.getbound(m3) <= 150.0

# historical simulation
results = simulate(m3, [:production, :contracts],
    noises=fill(2, 24), pricenoises=fill(10, 24)
)
@test isapprox(results[:objective], 181.0, atol=1)
results = simulate(m3, 500)
@test isapprox(mean(r[:objective] for r in results), 177.0, atol=1.0)

# 5 fixed ribs
m5 = contracting_example(5)
srand(123)
SDDP.solve(m5, iteration_limit = 10, print_level=0)
@test SDDP.getbound(m5) <= 150.0

# historical simulation
results = simulate(m5, [:production, :contracts],
    noises=fill(2, 24), pricenoises=fill(10, 24)
)
@test isapprox(results[:objective], 182.0, atol=1)
results = simulate(m5, 500)
@test isapprox(mean(r[:objective] for r in results), 179.0, atol=1.0)

@test SDDP.getbound(m5) <= SDDP.getbound(m3)
