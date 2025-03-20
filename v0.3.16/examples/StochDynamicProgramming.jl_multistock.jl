#  Copyright 2017-21, Oscar Dowson.                                     #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # StochDynamicProgramming: the multistock problem

# This example comes from [StochDynamicProgramming.jl](https://github.com/JuliaOpt/StochDynamicProgramming.jl/tree/f68b9da541c2f811ce24fc76f6065803a0715c2f/examples/multistock-example.jl).

using SDDP, GLPK, Test

function test_multistock_example()
    model = SDDP.LinearPolicyGraph(
        stages = 5,
        lower_bound = -5.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, stage
        @variable(
            subproblem,
            0 <= stock[i = 1:3] <= 1,
            SDDP.State,
            initial_value = 0.5
        )
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
        Ξ = collect(
            Base.product((0.0, 0.15, 0.3), (0.0, 0.15, 0.3), (0.0, 0.15, 0.3)),
        )[:]
        SDDP.parameterize(subproblem, Ξ) do ω
            return JuMP.fix.(ξ, ω)
        end
        @stageobjective(subproblem, (sin(3 * stage) - 1) * sum(control))
    end
    SDDP.train(
        model,
        iteration_limit = 100,
        cut_type = SDDP.SINGLE_CUT,
        log_frequency = 10,
    )
    @test SDDP.calculate_bound(model) ≈ -4.349 atol = 0.01

    simulation_results = SDDP.simulate(model, 5000)
    @test length(simulation_results) == 5000
    μ = SDDP.Statistics.mean(
        sum(data[:stage_objective] for data in simulation) for
        simulation in simulation_results
    )
    @test μ ≈ -4.349 atol = 0.1
    return
end

test_multistock_example()
