#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    An implementation of the QuickStart example from FAST
    https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/demo
=#

using SDDP, GLPK, Test

function fast_quickstart()
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        bellman_function = SDDP.BellmanFunction(lower_bound = -5),
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x >= 0, SDDP.State, initial_value = 0.0)
        if t == 1
            @stageobjective(sp, x.out)
        else
            @variable(sp, s >= 0)
            @constraint(sp, s <= x.in)
            SDDP.parameterize(sp, [2, 3]) do ω
                JuMP.set_upper_bound(s, ω)
            end
            @stageobjective(sp, -2s)
        end
    end

    det = SDDP.deterministic_equivalent(model, GLPK.Optimizer)
    JuMP.optimize!(det)
    @test JuMP.objective_value(det) == -2

    SDDP.train(model, iteration_limit = 3, print_level = 0)
    @test SDDP.calculate_bound(model) == -2
end

fast_quickstart()
