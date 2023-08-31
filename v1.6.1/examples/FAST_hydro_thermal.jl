#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # FAST: the hydro-thermal problem

# An implementation of the Hydro-thermal example from [FAST](https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/hydro%20thermal)

using SDDP, HiGHS, Test

function fast_hydro_thermal()
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
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
            return JuMP.fix(ξ, ω)
        end
        @stageobjective(sp, 5 * p)
    end

    det = SDDP.deterministic_equivalent(model, HiGHS.Optimizer)
    set_silent(det)
    JuMP.optimize!(det)
    @test JuMP.objective_value(det) == 10
    SDDP.train(model)
    @test SDDP.calculate_bound(model) == 10
    return
end

fast_hydro_thermal()
