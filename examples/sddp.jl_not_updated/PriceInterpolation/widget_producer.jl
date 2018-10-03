#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

"""
    widget_producer_example(DISCRETIZATION = 1)

Create an instance of the widget producer example. `DISCRETIZATION` sets the
number of static ribs in the price dimension. The default, `DISCRETIZATION`=1`
uses the dynamic method instead.
"""
function widget_producer_example(DISCRETIZATION = 1)

    srand(10)
    # number of stages
    T = 12
    # Initial Price
    P₀ = 1.0
    # Long-run mean
    μ = 1.0
    # Mean-reversion factor
    ρ = 0.05
    # distribution of price noises
    Φ = DiscreteDistribution([-0.050, -0.025, 0.0, 0.025, 0.050])
    # Minimum price
    Pₗ  = 0.0
    # Maximum price
    Pᵤ = 2.0
    # Price dynamics
    function pricedynamics(price, noise, stage, markov)
        if stage == 1 # Stage 1 is deterministic
            return price
        end
        return price - ρ * (price - μ) + noise
    end

    value_function = if DISCRETIZATION == 1
        (t,i) -> DynamicPriceInterpolation(
            dynamics       = (p,w) -> pricedynamics(p,w,t,i),
            initial_price  = P₀,
            min_price      = Pₗ,
            max_price      = Pᵤ,
            noise          = Φ
        )
    else
        (t,i) -> StaticPriceInterpolation(
            dynamics       = (p,w) -> pricedynamics(p,w,t,i),
            initial_price  = P₀,
            # we need to ensure the ribs are not binding on the edge of the
            # price domain because sometimes the duals can be weird / numeric
            # error happens in the interpolation. I still haven't isolated the
            # exact cause, this this seems to fix the problem.
            rib_locations  =  collect(linspace(Pₗ, Pᵤ, DISCRETIZATION)),
            noise          = Φ
        )
    end

    m = SDDPModel(
        sense             = :Max,
        stages            = T,
        objective_bound   = 50.0,
        solver            = ClpSolver(),
        # risk_measure      = EAVaR(lambda=0.5, beta=0.25),
        value_function    = value_function
                                            ) do sp, t
        #=
        Data that defines problem size
        =#
        # Perishability
        λ = [0.98, 0.98]
        P = length(λ)
        # Contango
        γ = [1.0, 1.1]
        M = length(γ)
        # production shut-down probability
        Ω = DiscreteDistribution([1.0, 0.0], [0.98, 0.02])
        # transaction costs
        ψ = 0.01
        # spot premium to buy
        κ = 5.0
        # max contracts to sell
        cₘ = 10

        @states(sp, begin
            C[1:M] >= 0, C0 == 0 # number of contracts in future months 1,2,M
            S[1:P] >= 0, S0 == 0 # stock levels in perishability levels
            w      >= 0, w0 == 1 # factory is operating 1=true, 0=false
        end)
        @variables(sp, begin
            d[1:P] >= 0 # deliver to contracts from perishability level p
            s[1:P] >= 0 # sell on spot market from perishability level p
            c[1:M] >= 0 # sell contracts m months in the future
            b      >= 0 # buy widgets from spot
        end)
        @constraints(sp, begin
            # age contracts + new sales
            [m=1:(M-1)], C[m] == C0[m+1] + c[m]
            # in final month, only new sales
            C[M] == c[M]
            # cannot sell more than cₘ
            c .<= cₘ

            # factory cannot restart once shut
            w <= w0

            # initial stock is production plus buys,
            # less spot sales and deliveries.
            # if factory is operating (w=1), one unit is produced, otherwise,
            # 0 units are produced.
            S[1] == w - s[1] - d[1] + b
            # age the stock less spot sales and deliveries
            [p=2:P], S[p] == λ[p] * S0[p-1] - s[p] - d[p]

            # ensure deliveries meet contracted amount
            C0[1] == sum(d[p] for p in 1:P)

            # can't offer a contract past final stage
            [m=1:M; m>T-t], c[m] == 0
        end)
        if t != 1
            # first stage is deterministic
            # otherwise there is a smal probabilty of shut-down
            @rhsnoise(sp, ω=Ω, w <= SDDP.observation(ω))
            setnoiseprobability!(sp, [SDDP.probability(ω) for ω in Ω])
        end

        @stageobjective(sp, price -> price * (
                # spot sales
                sum(s[p] for p in 1:P) +
                # contracts with contango
                sum(γ[m] * c[m] for m in 1:M) -
                # spot buys
                κ * b
                # transaction cost
            ) - ψ * sum(c[m] for m in 1:M)
        )
    end
    return m
end

m = widget_producer_example()
srand(123)
status = SDDP.solve(m,
    time_limit     = 5.0,
    simulation = MonteCarloSimulation(
        frequency = 100, min=100, step=100,max=1000
    ),
    print_level = 0
)
@test status == :time_limit
results = simulate(m, 1000, [:C, :S, :price, :s, :d, :b])
@test length(results) == 1000
@test getbound(m) <= 15.0

# SDDP.save!("results.julia", results)
# plt = SDDP.newplot()
# SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->results[i][:stageobjective][t], cumulative=true, title="Accumulation of Profit")
# SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->results[i][:price][t], title="Price")
# SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->sum(results[i][:S][t]), title="Stockpiled Widgets")
# SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->results[i][:b][t], title="Purchases")
# SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->sum(results[i][:s][t]), title="Sales")
# SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->results[i][:S][t][1] + results[i][:s][t][1] + results[i][:d][t][1] - results[i][:b][t], title="Production")
# for c in 1:2
#     SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->results[i][:C][t][c], title="C$(c)")
# end
# for p in 1:2
#     SDDP.addplot!(plt, 1:1000, 1:12, (i, t)->results[i][:S][t][p], title="S$(p)")
# end
# SDDP.show(plt)
# rm("results.julia")
