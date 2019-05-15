#  Copyright 2019, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This example is derived from Section 4.2 of the paper:
#   Ahmed, S., Cabral, F. G., & da Costa, B. F. P. (2019). Stochastic Lipschitz
#   Dynamic Programming. Optimization Online.
#   URL: http://www.optimization-online.org/DB_FILE/2019/05/7193.pdf

using SDDP, GLPK, Test

function sldp_example_one()
    model = SDDP.LinearPolicyGraph(
            stages = 8,
            lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
            ) do sp, t
        @variable(sp, x, SDDP.State, initial_value = 2.0)
        @variables(sp, begin
            x⁺ >= 0
            x⁻ >= 0
            0 <= u <= 1, Bin
            ω
        end)
        @stageobjective(sp, 0.9^(t-1) * (x⁺ + x⁻))
        @constraints(sp, begin
            x.out == x.in + 2 * u - 1 + ω
            x⁺ >= x.out
            x⁻ >= x.in
        end)
        points = [-0.5, 0.3, 0.4, 1.1, 1.5]
        SDDP.parameterize(φ -> JuMP.fix(ω, φ), sp, [points; -points])
    end
    SDDP.train(model, iteration_limit = 100, print_level = 0)
    # TODO(odow): include the actual set of points when known, and update this
    # bound. The paper reports a bound of [3.085, 3.313].
    @test SDDP.calculate_bound(model) ≈ 5.49 atol=0.1
end

sldp_example_one()
