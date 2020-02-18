#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Test

function no_strong_duality()
    model = SDDP.PolicyGraph(
        SDDP.Graph(:root, [:node], [(:root => :node, 1.0), (:node => :node, 0.5)]),
        optimizer = GLPK.Optimizer,
        lower_bound = 0.0,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 1.0)
        @stageobjective(sp, x.out)
        @constraint(sp, x.in == x.out)
    end

    SDDP.train(model, iteration_limit = 20, print_level = 0)

    @test SDDP.calculate_bound(model) ≈ 2.0 atol = 1e-8
end

no_strong_duality()
