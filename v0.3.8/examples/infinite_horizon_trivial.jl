#  Copyright 2017-21, Oscar Dowson.                                     #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Infinite horizon trivial

using SDDP, GLPK, Test

function infinite_trivial()
    graph = SDDP.Graph(
        :root_node, [:week], [(:root_node => :week, 1.0), (:week => :week, 0.9)]
    )
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, node
        @variable(subproblem, state, SDDP.State, initial_value = 0)
        @constraint(subproblem, state.in == state.out)
        @stageobjective(subproblem, 2.0)
    end
    SDDP.train(model; iteration_limit = 100, log_frequency = 10)
    @test SDDP.calculate_bound(model) ≈ 2.0 / (1 - 0.9) atol = 1e-3
    return
end

infinite_trivial()
