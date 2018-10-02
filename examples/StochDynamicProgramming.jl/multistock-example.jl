#==
#  Copyright 2018, Oscar Dowson
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

using Kokako, GLPK, Test

function test_multistock_example()
    model = Kokako.PolicyGraph(Kokako.LinearGraph(5),
            optimizer = with_optimizer(GLPK.Optimizer),
            bellman_function = Kokako.AverageCut(lower_bound=-5)
                                    ) do subproblem, stage
        @variables(subproblem, begin
            stock[i=1:3]
            0 <= stock′[i=1:3] <= 1
            0 <= control[i=1:3] <= 0.5
            ξ[i=1:3]  # Dummy for RHS noise.
        end)
        @constraints(subproblem, begin
            sum(control) - 0.5 * 3 <= 0
            stock′ .== stock .+ control .- ξ
        end)
        Kokako.add_state_variable(subproblem, :stock_1, stock[1], stock′[1])
        Kokako.add_state_variable(subproblem, :stock_2, stock[2], stock′[2])
        Kokako.add_state_variable(subproblem, :stock_3, stock[3], stock′[3])
        Ξ = collect(Base.product(
                (0.0, 0.15, 0.3),
                (0.0, 0.15, 0.3),
                (0.0, 0.15, 0.3))
            )[:]
        Kokako.parameterize(subproblem, Ξ) do ω
            JuMP.fix.(ξ, ω)
        end
        Kokako.set_stage_objective(subproblem, :Min,
            (sin(3 * stage) - 1) * sum(control)
        )
    end
    initial_state = Dict(
        :stock_1 => 0.5,
        :stock_2 => 0.5,
        :stock_3 => 0.5
    )
    status = Kokako.train(model,
        iteration_limit = 100,
        print_level = 0,
        initial_state = initial_state
        )
    @test Kokako.calculate_bound(model, initial_state) ≈ -4.349 atol=0.01
    # results = simulate(m, 5000)
    # @test length(results) == 5000
    # @test isapprox(mean(r[:objective] for r in results), -4.349, atol=0.02)
end

test_multistock_example()
