#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
The Asset Management problem taken from

    J. R. Birge,  F. Louveaux,  Introduction to Stochastic Programming,
    Springer Series in Operations Research and Financial Engineering,
    Springer New York, New York, NY, 2011
=#

using SDDP, GLPK, Test

function asset_management_simple()
    model = SDDP.PolicyGraph(
        SDDP.MarkovianGraph(Array{Float64,2}[
            [1.0]',
            [0.5 0.5],
            [0.5 0.5; 0.5 0.5],
            [0.5 0.5; 0.5 0.5],
        ]),
        bellman_function = SDDP.BellmanFunction(lower_bound = -1_000.0),
        optimizer = GLPK.Optimizer,
    ) do subproblem, index
        (stage, markov_state) = index
        rstock = [1.25, 1.06]
        rbonds = [1.14, 1.12]
        @variable(subproblem, stocks >= 0, SDDP.State, initial_value = 0.0)
        @variable(subproblem, bonds >= 0, SDDP.State, initial_value = 0.0)
        if stage == 1
            @constraint(subproblem, stocks.out + bonds.out == 55)
            @stageobjective(subproblem, 0)
        elseif 1 < stage < 4
            @constraint(
                subproblem,
                rstock[markov_state] * stocks.in + rbonds[markov_state] * bonds.in ==
                stocks.out + bonds.out
            )
            @stageobjective(subproblem, 0)
        else
            @variable(subproblem, over >= 0)
            @variable(subproblem, short >= 0)
            @constraint(
                subproblem,
                rstock[markov_state] * stocks.in + rbonds[markov_state] * bonds.in - over + short == 80
            )
            @stageobjective(subproblem, -over + 4 * short)
        end
    end
    SDDP.train(model, iteration_limit = 25, print_level = 0)
    @test SDDP.termination_status(model) == :iteration_limit
    @test SDDP.calculate_bound(model) ≈ 1.514 atol = 1e-4
end

asset_management_simple()
