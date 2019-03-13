#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#=
    This example comes from github.com/JuliaOpt/StochDynamicProgramming.jl/tree/
        f68b9da541c2f811ce24fc76f6065803a0715c2f/examples/stock-example.jl

    Vincent Leclere, Francois Pacaud and Henri Gerard

 Compare different ways of solving a stock problem :
 Min   E [\sum_{t=1}^TF c_t u_t]
 s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
         0 <= s_t <= 1
         u_min <= u_t <= u_max
         u_t choosen knowing xi_1 .. xi_t

=#

using Kokako, GLPK, Test

function stock_example()
    model = Kokako.PolicyGraph(Kokako.LinearGraph(5),
                bellman_function = Kokako.AverageCut(lower_bound = -2),
                optimizer = with_optimizer(GLPK.Optimizer)
                                    ) do sp, stage

        # @state(sp, 0 <= state <= 1, state0 == 0.5)
        @variable(sp, 0 <= state <= 1, Kokako.State, initial_value = 0.5)
        @variable(sp, 0 <= control <= 0.5)
        @variable(sp, ξ)
        @constraint(sp, state.out == state.in - control + ξ)
        Kokako.parameterize(sp, 0.0:1/30:0.3) do ω
            JuMP.fix(ξ, ω)
        end
        @stageobjective(sp, (sin(3 * stage) - 1) * control)
    end

    Kokako.train(model, iteration_limit = 50, print_level = 0)
    @test Kokako.calculate_bound(model) ≈ -1.471 atol=0.001

    simulation_results = Kokako.simulate(model, 1_000)
    @test length(simulation_results) == 1_000
    @test Kokako.Statistics.mean(
        sum(data[:stage_objective] for data in simulation)
        for simulation in simulation_results
    ) ≈ -1.471 atol=0.05
end

stock_example()
