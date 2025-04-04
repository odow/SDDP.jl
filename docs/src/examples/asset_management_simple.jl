#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Asset management

# Taken from the book
# J.R. Birge, F. Louveaux, Introduction to Stochastic Programming,
# Springer Series in Operations Research and Financial Engineering,
# Springer New York, New York, NY, 2011

using SDDP, HiGHS, Test

function asset_management_simple()
    model = SDDP.PolicyGraph(
        SDDP.MarkovianGraph(
            Array{Float64,2}[
                [1.0]',
                [0.5 0.5],
                [0.5 0.5; 0.5 0.5],
                [0.5 0.5; 0.5 0.5],
            ],
        );
        lower_bound = -1_000.0,
        optimizer = HiGHS.Optimizer,
    ) do subproblem, index
        (stage, markov_state) = index
        r_stock = [1.25, 1.06]
        r_bonds = [1.14, 1.12]
        @variable(subproblem, stocks >= 0, SDDP.State, initial_value = 0.0)
        @variable(subproblem, bonds >= 0, SDDP.State, initial_value = 0.0)
        if stage == 1
            @constraint(subproblem, stocks.out + bonds.out == 55)
            @stageobjective(subproblem, 0)
        elseif 1 < stage < 4
            @constraint(
                subproblem,
                r_stock[markov_state] * stocks.in +
                r_bonds[markov_state] * bonds.in == stocks.out + bonds.out
            )
            @stageobjective(subproblem, 0)
        else
            @variable(subproblem, over >= 0)
            @variable(subproblem, short >= 0)
            @constraint(
                subproblem,
                r_stock[markov_state] * stocks.in +
                r_bonds[markov_state] * bonds.in - over + short == 80
            )
            @stageobjective(subproblem, -over + 4 * short)
        end
    end
    SDDP.train(model; log_frequency = 5)
    @test SDDP.calculate_bound(model) ≈ 1.514 atol = 1e-4
    return
end

asset_management_simple()
