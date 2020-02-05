#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    An implementation of the Production Management example from FAST
    https://github.com/leopoldcambier/FAST/blob/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/production management multiple stages/
=#

using SDDP, GLPK, Test

function fast_production_management(; cut_type)
    DEMAND = [2, 10]
    H = 3
    N = 2
    C = [0.2, 0.7]
    S = 2 .+ [0.33, 0.54]
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(H),
        bellman_function = SDDP.BellmanFunction(lower_bound = -50.0, cut_type = cut_type),
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x[1:N] >= 0, SDDP.State, initial_value = 0.0)
        @variables(sp, begin
            s[i = 1:N] >= 0
            d
        end)
        @constraints(sp, begin
            [i = 1:N], s[i] <= x[i].in
            sum(s) <= d
        end)
        SDDP.parameterize(sp, t == 1 ? [0] : DEMAND) do ω
            JuMP.fix(d, ω)
        end
        @stageobjective(sp, sum(C[i] * x[i].out for i = 1:N) - S's)
    end
    SDDP.train(model, iteration_limit = 10, print_level = 2)
    @test SDDP.calculate_bound(model) ≈ -23.96 atol = 1e-2
end

fast_production_management(cut_type = SDDP.SINGLE_CUT)
fast_production_management(cut_type = SDDP.MULTI_CUT)
