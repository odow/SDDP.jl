#  Copyright 2017-19, Oscar Dowson.
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

using GLPK
using SDDP
using Test

function air_conditioning_model(integrality_handler)
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = with_optimizer(GLPK.Optimizer),
        integrality_handler = integrality_handler
    ) do sp, stage
        @variable(sp, 0 <= stored_production <= 100, Int, SDDP.State, initial_value = 0)
        @variable(sp, 0 <= production <= 200, Int)
        @variable(sp, overtime >= 0, Int)
        @constraint(
            sp,
            demand_con,
            stored_production.out == stored_production.in + production + overtime
        )
        @stageobjective(
            sp,
            100 * production + 300 * overtime + 50 * stored_production.out
        )
        DEMAND = [[100.0], [100.0, 300.0], [100.0, 300.0]]
        SDDP.parameterize(sp, DEMAND[stage]) do w
            JuMP.set_normalized_rhs(demand_con, -w)
        end
    end
    SDDP.train(model, iteration_limit = 30, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 62_500.0
end

air_conditioning_model(SDDP.SDDiP())
air_conditioning_model(SDDP.ContinuousRelaxation())
