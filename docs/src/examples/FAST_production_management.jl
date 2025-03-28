#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # FAST: the production management problem

# An implementation of the Production Management example from [FAST](https://github.com/leopoldcambier/FAST/blob/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/production management multiple stages/)

using SDDP, HiGHS, Test

function fast_production_management(; cut_type)
    DEMAND = [2, 10]
    H = 3
    N = 2
    C = [0.2, 0.7]
    S = 2 .+ [0.33, 0.54]
    model = SDDP.LinearPolicyGraph(;
        stages = H,
        lower_bound = -50.0,
        optimizer = HiGHS.Optimizer,
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
            return JuMP.fix(d, ω)
        end
        @stageobjective(sp, sum(C[i] * x[i].out for i in 1:N) - S's)
    end
    SDDP.train(model; cut_type = cut_type, print_level = 2, log_frequency = 5)
    @test SDDP.calculate_bound(model) ≈ -23.96 atol = 1e-2
end

fast_production_management(; cut_type = SDDP.SINGLE_CUT)
fast_production_management(; cut_type = SDDP.MULTI_CUT)
