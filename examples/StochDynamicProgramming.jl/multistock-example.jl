#==
#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

    This example comes from github.com/JuliaOpt/StochDynamicProgramming.jl/tree/
        f68b9da541c2f811ce24fc76f6065803a0715c2f/examples/multistock-example.jl

    Vincent Leclere, Francois Pacaud and Henri Gerard

 Compare different ways of solving a stock problem :
 Min   E [\sum_{t=1}^TF c_t u_t]
 s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
         0 <= s_t <= 1
         u_min <= u_t <= u_max
         u_t choosen knowing xi_1 .. xi_t

==##############################################################################

using SDDP, JuMP, Clp, Base.Test

srand(100)

XI = collect(Base.product([linspace(0, 0.3, 3) for i in 1:3]...))[:]

m = SDDPModel(
                  sense = :Min,
                 stages = 5,
                 solver = ClpSolver(),
        objective_bound = -5
                                ) do sp, stage

    @state(sp, 0 <= stock[i=1:3] <= 1, stock0 == 0.5)

    @variable(sp, 0 <= control[i=1:3] <= 0.5)

    @constraint(sp, sum(control) - 0.5 * 3 <= 0)

    @noises(sp, xi = XI, begin
        stock[1] == stock0[1] + control[1] - xi[1]
        stock[2] == stock0[2] + control[2] - xi[2]
        stock[3] == stock0[3] + control[3] - xi[3]
    end)

    stageobjective!(sp, (sin(3 * stage) - 1) * sum(control))
end

@time status = SDDP.solve(m, max_iterations = 100,
print_level = 0)
@test isapprox(SDDP.getbound(m), -4.349, atol=0.01)

results = simulate(m, 5000)

@test isapprox(mean(r[:objective] for r in results), -4.349, atol=0.02)
