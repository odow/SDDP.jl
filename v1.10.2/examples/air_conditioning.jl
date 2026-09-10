#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Air conditioning

# Taken from [Anthony Papavasiliou's notes on SDDP](https://web.archive.org/web/20200504214809/https://perso.uclouvain.be/anthony.papavasiliou/public_html/SDDP.pdf)
# This is a variation of the problem that first appears in the book 
# Introduction to Stochastic Programming by Birge and Louveaux, 1997, 
# Springer-Verlag, New York, on page 237, Example 1. For a rescaled problem,
# they reported an optimal value of 6.25 with a first-stage solution of x1 = 2 
# (production)and y1 = 1 (store production). On this variation, without rescaling, 
# it would be equivalent to 62500, 200 and 100, respectively.

# Consider the following problem
# * Produce air conditioners for 3 months
# * 200 units/month at 100 \$/unit
# * Overtime costs 300 \$/unit
# * Known demand of 100 units for period 1
# * Equally likely demand, 100 or 300 units, for periods 2, 3
# * Storage cost is 50 \$/unit
# * All demand must be met

# The known optimal solution is \$62,500

using SDDP, HiGHS, Test

function air_conditioning_model(duality_handler)
    model = SDDP.LinearPolicyGraph(;
        stages = 3,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, stage
        @variable(
            sp,
            0 <= stored_production <= 100,
            Int,
            SDDP.State,
            initial_value = 0
        )
        @variable(sp, 0 <= production <= 200, Int)
        @variable(sp, overtime >= 0, Int)
        @variable(sp, demand)
        DEMAND = [[100.0], [100.0, 300.0], [100.0, 300.0]]
        SDDP.parameterize(ω -> JuMP.fix(demand, ω), sp, DEMAND[stage])
        @constraint(
            sp,
            stored_production.out ==
            stored_production.in + production + overtime - demand
        )
        @stageobjective(
            sp,
            100 * production + 300 * overtime + 50 * stored_production.out
        )
    end
    SDDP.train(model; duality_handler = duality_handler)
    lb = SDDP.calculate_bound(model)
    println("Lower bound is: $lb")
    @test isapprox(lb, 62_500.0, atol = 0.1)
    sims = SDDP.simulate(model, 1, [:production, :stored_production])
    x1 = sims[1][1][:production]
    y1 = sims[1][1][:stored_production].out
    @test isapprox(x1, 200, atol = 0.1)
    @test isapprox(y1, 100, atol = 0.1)
    println(
        "With first stage solutions $(x1) (production) and $(y1) (stored_production).",
    )
    return
end

for duality_handler in [SDDP.LagrangianDuality(), SDDP.ContinuousConicDuality()]
    air_conditioning_model(duality_handler)
end
