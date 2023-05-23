#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Air conditioning

# Taken from [Anthony Papavasiliou's notes on SDDP](https://web.archive.org/web/20200504214809/https://perso.uclouvain.be/anthony.papavasiliou/public_html/SDDP.pdf)

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
    model = SDDP.LinearPolicyGraph(
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
    @test isapprox(SDDP.calculate_bound(model), 62_500.0, atol = 0.1)
    return
end

for duality_handler in [SDDP.LagrangianDuality(), SDDP.ContinuousConicDuality()]
    air_conditioning_model(duality_handler)
end
