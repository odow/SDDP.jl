#  Copyright 2017-20, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This example is derived from Section 4.3 of the paper:
#   Ahmed, S., Cabral, F. G., & da Costa, B. F. P. (2019). Stochastic Lipschitz
#   Dynamic Programming. Optimization Online.
#   URL: http://www.optimization-online.org/DB_FILE/2019/05/7193.pdf

using SDDP, GLPK, Test

function sldp_example_two(; first_stage_integer::Bool = true, N = 2)
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = -100.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x[1:2] <= 5, SDDP.State, initial_value = 0.0)
        if t == 1
            if first_stage_integer
                @variable(sp, 0 <= u[1:2] <= 5, Int)
                @constraint(sp, [i = 1:2], u[i] == x[i].out)
            end
            @stageobjective(sp, -1.5 * x[1].out - 4 * x[2].out)
        else
            @variable(sp, 0 <= y[1:4] <= 1, Bin)
            @variable(sp, ω[1:2])
            @stageobjective(sp, -16 * y[1] - 19 * y[2] - 23 * y[3] - 28 * y[4])
            @constraint(sp, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= ω[1] - x[1].in)
            @constraint(sp, 6 * y[1] + 1 * y[2] + 3 * y[3] + 2 * y[4] <= ω[2] - x[2].in)
            steps = range(5, stop = 15, length = N)
            SDDP.parameterize(sp, [[i, j] for i in steps for j in steps]) do φ
                JuMP.fix.(ω, φ)
            end
        end
    end
    SDDP.train(model, iteration_limit = 100, print_level = 0)
    bound = SDDP.calculate_bound(model)

    if N == 2
        @test bound <= -57.0
    elseif N == 3
        @test bound <= -59.33
    elseif N == 6
        @test bound <= -61.22
    end
end

sldp_example_two(N = 2)
sldp_example_two(N = 3)
sldp_example_two(N = 6)
