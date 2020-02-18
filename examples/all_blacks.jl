#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Test

function all_blacks()
    # Number of time periods, number of seats, R_ij = evenue from selling seat i
    # at time j, offer_ij = whether an offer for seat i will come at time j
    (T, N, R, offer) = (3, 2, [3 3 6; 3 3 6], [1 1 0; 1 0 1])

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
        @variable(sp, accept_offer, Bin)
        # Balance on seats
        @constraint(sp, [i in 1:N], x[i].out == x[i].in - offer[i, stage] * accept_offer)
        @stageobjective(sp, sum(R[i, stage] * offer[i, stage] * accept_offer for i = 1:N))
    end
    SDDP.train(model, iteration_limit = 10, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 9.0
end

all_blacks()
