#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Statistics, Test

function biobjective_hydro()
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(3),
        bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
        optimizer = GLPK.Optimizer,
    ) do subproblem, stage
        @variable(subproblem, 0 <= v <= 200, SDDP.State, initial_value = 50)

        SDDP.add_objective_state(
            subproblem,
            initial_value = 0.0,
            lower_bound = 0.0,
            upper_bound = 1.0,
            lipschitz = 1e6,
        ) do y, ω
            return y + ω.λ
        end

        @variables(subproblem, begin
            0 <= g[i = 1:2] <= 100
            0 <= u <= 150
            s >= 0
            a
            shortage_cost >= 0
            objective_1
            objective_2
        end)
        @constraints(subproblem, begin
            v.out == v.in - u - s + a
            demand_eq, g[1] + g[2] + u == 150
            objective_1 == g[1] + 10 * g[2]
            shortage_cost >= 40 - v.out
            shortage_cost >= 60 - 2 * v.out
            shortage_cost >= 80 - 4 * v.out
            objective_2 == shortage_cost
        end)
        price_noise_terms = (stage == 1) ? [0.04, 0.12, 0.22, 0.64] : [0.0]
        Ω = [(a = i, λ = j) for i = 0.0:5:50.0 for j in price_noise_terms]
        SDDP.parameterize(subproblem, Ω) do ω
            JuMP.fix(a, ω.a)
            # This *has* to be called from inside `SDDP.parameterize`,
            # otherwise it doesn't make sense.
            λ = SDDP.objective_state(subproblem)
            @stageobjective(subproblem, λ * objective_1 + (1 - λ) * objective_2)
        end
    end
    SDDP.train(model, iteration_limit = 50, print_level = 0)

    results = SDDP.simulate(model, 500)
    objectives = [sum(s[:stage_objective] for s in simulation) for simulation in results]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    println("Confidence_interval = $(sample_mean) ± $(sample_ci)")
    @test SDDP.calculate_bound(model) ≈ sample_mean atol = sample_ci
end

biobjective_hydro()
