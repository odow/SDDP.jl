#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    An implementation of the Hydro-thermal example from FAST
    https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/hydro%20thermal
=#

using SDDP, GLPK, Test

function fast_hydro_thermal()
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 8, SDDP.State, initial_value = 0.0)
        @variables(sp, begin
            y >= 0
            p >= 0
            ξ
        end)
        @constraints(sp, begin
            p + y >= 6
            x.out <= x.in - y + ξ
        end)
        RAINFALL = (t == 1 ? [6] : [2, 10])
        SDDP.parameterize(sp, RAINFALL) do ω
            JuMP.fix(ξ, ω)
        end
        @stageobjective(sp, 5 * p)
    end

    det = SDDP.deterministic_equivalent(model, GLPK.Optimizer)
    JuMP.optimize!(det)
    @test JuMP.objective_value(det) == 10

    SDDP.train(model, iteration_limit = 10, print_level = 0)
    @test SDDP.calculate_bound(model) == 10
end

fast_hydro_thermal()
