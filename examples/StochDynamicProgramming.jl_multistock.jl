#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    This example comes from github.com/JuliaOpt/StochDynamicProgramming.jl/tree/
        f68b9da541c2f811ce24fc76f6065803a0715c2f/examples/multistock-example.jl

    Vincent Leclere, Francois Pacaud and Henri Gerard

 Compare different ways of solving a stock problem :
 Min   E [\sum_{t=1}^TF c_t u_t]
 s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
         0 <= s_t <= 1
         u_min <= u_t <= u_max
         u_t choosen knowing xi_1 .. xi_t

=#

using SDDP, GLPK, Test

function test_multistock_example()
    model = SDDP.LinearPolicyGraph(
        stages = 5,
        lower_bound = -5.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, stage
        @variable(subproblem, 0 <= stock[i = 1:3] <= 1, SDDP.State, initial_value = 0.5)
        @variables(subproblem, begin
            0 <= control[i = 1:3] <= 0.5
            ξ[i = 1:3]  # Dummy for RHS noise.
        end)
        @constraints(
            subproblem,
            begin
                sum(control) - 0.5 * 3 <= 0
                [i = 1:3], stock[i].out == stock[i].in + control[i] - ξ[i]
            end
        )
        Ξ = collect(Base.product((0.0, 0.15, 0.3), (0.0, 0.15, 0.3), (0.0, 0.15, 0.3)))[:]
        SDDP.parameterize(subproblem, Ξ) do ω
            JuMP.fix.(ξ, ω)
        end
        @stageobjective(subproblem, (sin(3 * stage) - 1) * sum(control))
    end
    SDDP.train(model, iteration_limit = 100, print_level = 0, cut_type = SDDP.SINGLE_CUT)
    @test SDDP.calculate_bound(model) ≈ -4.349 atol = 0.01

    simulation_results = SDDP.simulate(model, 5000)
    @test length(simulation_results) == 5000
    @test SDDP.Statistics.mean(
        sum(data[:stage_objective] for data in simulation)
        for simulation in simulation_results
    ) ≈ -4.349 atol = 0.02
end

test_multistock_example()
