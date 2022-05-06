#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # SLDP: example 1

# This example is derived from Section 4.2 of the paper:
# Ahmed, S., Cabral, F. G., & da Costa, B. F. P. (2019). Stochastic Lipschitz
# Dynamic Programming. Optimization Online. [PDF](http://www.optimization-online.org/DB_FILE/2019/05/7193.pdf)

using SDDP, GLPK, Test

function sldp_example_one()
    model = SDDP.LinearPolicyGraph(
        stages = 8,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 2.0)
        @variables(sp, begin
            x⁺ >= 0
            x⁻ >= 0
            0 <= u <= 1, Bin
            ω
        end)
        @stageobjective(sp, 0.9^(t - 1) * (x⁺ + x⁻))
        @constraints(sp, begin
            x.out == x.in + 2 * u - 1 + ω
            x⁺ >= x.out
            x⁻ >= -x.out
        end)
        points = [
            -0.3089653673606697,
            -0.2718277412744214,
            -0.09611178608243474,
            0.24645863921577763,
            0.5204224537256875,
        ]
        return SDDP.parameterize(φ -> JuMP.fix(ω, φ), sp, [points; -points])
    end
    SDDP.train(model, iteration_limit = 100, log_frequency = 10)
    @test SDDP.calculate_bound(model) <= 1.1675
    return
end

sldp_example_one()
