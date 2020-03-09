#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    The air conditioning example from Anthony Papavasiliou
    https://perso.uclouvain.be/anthony.papavasiliou/public_html/SDDP.pdf

    Consider the following problem
        Produce air conditioners for 3 months
        200 units/month at 100 $/unit
        Overtime costs 300 $/unit
        Known demand of 100 units for period 1
        Equally likely demand, 100 or 300 units, for periods 2, 3
        Storage cost is 50 $/unit
        All demand must be met

    Optimal bound $62,500
=#
using SDDP, GLPK, Test

function air_conditioning_model(integrality_handler)
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
        integrality_handler = integrality_handler,
    ) do sp, stage
        @variable(sp, 0 <= stored_production <= 100, Int, SDDP.State, initial_value = 0)
        @variable(sp, 0 <= production <= 200, Int)
        @variable(sp, overtime >= 0, Int)
        @variable(sp, demand)
        DEMAND = [[100.0], [100.0, 300.0], [100.0, 300.0]]
        SDDP.parameterize(ω -> JuMP.fix(demand, ω), sp, DEMAND[stage])
        @constraint(
            sp,
            stored_production.out == stored_production.in + production + overtime - demand
        )
        @stageobjective(sp, 100 * production + 300 * overtime + 50 * stored_production.out)
    end
    SDDP.train(model, iteration_limit = 20, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 62_500.0
end

for integrality_handler in [SDDP.SDDiP(), SDDP.ContinuousRelaxation()]
    air_conditioning_model(integrality_handler)
end
