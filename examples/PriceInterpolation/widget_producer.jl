#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

function widget_producer_example(DISCRETIZATION = 1)

    srand(10)

    P₀ = 1.0
    μ = 1.0
    ρ = 0.05
    Φ = DiscreteDistribution([-0.050, -0.025, 0.0, 0.025, 0.050])
    Pₗ  = 0.0
    Pᵤ = 2.0

    function pricedynamics(price, noise, stage, markov)
        if stage == 1
            return price
        end
        return price - ρ * (price - μ) + noise
    end

    value_function = if DISCRETIZATION == 1
        DynamicPriceInterpolation(
            dynamics       = pricedynamics,
            initial_price  = P₀,
            min_price      = Pₗ,
            max_price      = Pᵤ,
            noise          = Φ
        )
    else
        StaticPriceInterpolation(
            dynamics       = pricedynamics,
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
        stages            = 12,
        objective_bound   = 50.0,
        solver            = ClpSolver(),
        # risk_measure      = NestedAVaR(lambda=0.5, beta=0.25),
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

        @states(sp, begin
            C[1:M] >= 0, C0 == 0
            S[1:P] >= 0, S0 == 0
            w      >= 0, w0 == 1
        end)
        @variables(sp, begin
            d[1:P] >= 0
            s[1:P] >= 0
            c[1:M] >= 0
            b      >= 0
        end)
        @constraints(sp, begin
            [m=1:(M-1)], C[m] == C0[m+1] + c[m]
            C[M] == c[M]
            c .<= 10

            w <= w0
            S[1] == w - s[1] - d[1] + b
            [p=2:P], S[p] == λ[p] * S0[p-1] - s[p] - d[p]

            C0[1] == sum(d[p] for p in 1:P)
            [m=1:M; m>12-t], c[m] == 0
        end)
        if t != 1
            @rhsnoise(sp, ω=Ω, w <= SDDP.observation(ω))
            setnoiseprobability!(sp, [SDDP.probability(ω) for ω in Ω])
        end
        @stageobjective(sp, price -> price * (
                sum(s[p] for p in 1:P) +
                sum(γ[m] * c[m] for m in 1:M) -
                κ * b
            ) - ψ * sum(c[m] for m in 1:M)
        )
    end
    return m
end

m = widget_producer_example()
status = SDDP.solve(m,
    time_limit     = 20.0,
    simulation = MonteCarloSimulation(
        frequency = 100, min=100, step=100,max=1000
    )
)
@test status == :time_limit
results = simulate(m, 1000, [:C, :S, :price, :s, :d, :b])
@test length(results) == 1000
@test getbound(m) <= 12.0

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
