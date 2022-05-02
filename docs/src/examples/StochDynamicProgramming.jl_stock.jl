#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors.        #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # StochDynamicProgramming: the stock problem

# This example comes from [StochDynamicProgramming.jl](https://github.com/JuliaOpt/StochDynamicProgramming.jl/tree/f68b9da541c2f811ce24fc76f6065803a0715c2f/examples/stock-example.jl).

using SDDP, GLPK, Test

function stock_example()
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(5),
        lower_bound = -2,
        optimizer = GLPK.Optimizer,
    ) do sp, stage
        @variable(sp, 0 <= state <= 1, SDDP.State, initial_value = 0.5)
        @variable(sp, 0 <= control <= 0.5)
        @variable(sp, ξ)
        @constraint(sp, state.out == state.in - control + ξ)
        SDDP.parameterize(sp, 0.0:1/30:0.3) do ω
            return JuMP.fix(ξ, ω)
        end
        @stageobjective(sp, (sin(3 * stage) - 1) * control)
    end
    SDDP.train(model, iteration_limit = 50, log_frequency = 10)
    @test SDDP.calculate_bound(model) ≈ -1.471 atol = 0.001
    simulation_results = SDDP.simulate(model, 1_000)
    @test length(simulation_results) == 1_000
    μ = SDDP.Statistics.mean(
        sum(data[:stage_objective] for data in simulation) for
        simulation in simulation_results
    )
    @test μ ≈ -1.471 atol = 0.05
    return
end

stock_example()
