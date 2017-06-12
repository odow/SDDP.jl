#==
#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

    This example comes from github.com/JuliaOpt/StochDynamicProgramming.jl/tree/
        f68b9da541c2f811ce24fc76f6065803a0715c2f/examples/stock-example.jl

    Vincent Leclere, Francois Pacaud and Henri Gerard

 Compare different ways of solving a stock problem :
 Min   E [\sum_{t=1}^TF c_t u_t]
 s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
         0 <= s_t <= 1
         u_min <= u_t <= u_max
         u_t choosen knowing xi_1 .. xi_t

==##############################################################################

using SDDP, JuMP, Clp, Base.Test

m = SDDPModel(
                  sense = :Min,
                 stages = 5,
                 solver = ClpSolver(),
        objective_bound = -2
                                ) do sp, stage

    @state(sp, 0 <= state <= 1, state0 == 0.5)

    @variable(sp, 0 <= control <= 0.5)

    @noise(sp,
        noise = linspace(0, 0.3, 10),
        state == state0 + control - noise
    )

    stageobjective!(sp, (sin(3 * stage) - 1) * control)
end

@time status = SDDP.solve(m, max_iterations = 50)
@test isapprox(SDDP.getbound(m), -1.471, atol=0.001)

results = simulate(m, 1000)
@test isapprox(mean(r[:objective] for r in results), -1.471, atol=0.05)
