#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Training with a different forward model

using SDDP
import HiGHS
import Test

function create_air_conditioning_model(; convex::Bool)
    return SDDP.LinearPolicyGraph(;
        stages = 3,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do sp, t
        @variable(sp, 0 <= x <= 100, SDDP.State, initial_value = 0)
        @variable(sp, 0 <= u_production <= 200)
        @variable(sp, u_overtime >= 0)
        if !convex
            set_integer(x.out)
            set_integer(u_production)
            set_integer(u_overtime)
        end
        @constraint(sp, demand, x.in - x.out + u_production + u_overtime == 0)
        Ω = [[100.0], [100.0, 300.0], [100.0, 300.0]]
        SDDP.parameterize(ω -> JuMP.set_normalized_rhs(demand, ω), sp, Ω[t])
        @stageobjective(sp, 100 * u_production + 300 * u_overtime + 50 * x.out)
    end
end

convex = create_air_conditioning_model(; convex = true)
non_convex = create_air_conditioning_model(; convex = false)
SDDP.train(
    convex;
    forward_pass = SDDP.AlternativeForwardPass(non_convex),
    post_iteration_callback = SDDP.AlternativePostIterationCallback(non_convex),
    iteration_limit = 10,
)
Test.@test isapprox(SDDP.calculate_bound(non_convex), 62_500.0, atol = 0.1)
Test.@test isapprox(SDDP.calculate_bound(convex), 62_500.0, atol = 0.1)
