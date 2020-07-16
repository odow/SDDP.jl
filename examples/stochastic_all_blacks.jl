#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Test

function stochastic_all_blacks()
    # Number of time periods
    T = 3
    # Number of seats
    N = 2
    # R_ij = price of seat i at time j
    R = [3 3 6; 3 3 6]
    # Number of noises
    s = 3
    offers = [[[1, 1], [0, 0], [1, 1]], [[1, 0], [0, 0], [0, 0]], [[0, 1], [1, 0], [1, 1]]]

    model = SDDP.LinearPolicyGraph(
        stages = T,
        sense = :Max,
        upper_bound = 100.0,
        optimizer = GLPK.Optimizer,
        integrality_handler = SDDP.SDDiP(),
    ) do sp, stage

        # Seat remaining?
        @variable(sp, 0 <= x[1:N] <= 1, SDDP.State, Bin, initial_value = 1)
        # Action: accept offer, or don't accept offer
        # We are allowed to accept some of the seats offered but not others
        @variable(sp, accept_offer[1:N], Bin)
        @variable(sp, offers_made[1:N])
        # Balance on seats
        @constraint(sp, balance[i in 1:N], x[i].in - x[i].out == accept_offer[i])
        @stageobjective(sp, sum(R[i, stage] * accept_offer[i] for i = 1:N))
        SDDP.parameterize(sp, offers[stage]) do o
            JuMP.fix.(offers_made, o)
        end
        @constraint(sp, accept_offer .<= offers_made)
    end

    SDDP.train(model, iteration_limit = 10, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 8.0
end

stochastic_all_blacks()
