#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

#=
    Example: newsvendor.

    This example is based on the classical newsvendor problem, but features
    an AR(1) spot-price.

    V(x[t-1], ω[t]) =         max p[t] × u[t]
                       subject to x[t] = x[t-1] - u[t] + ω[t]
                                  u[t] ∈ [0, 1]
                                  x[t] ≥ 0
                                  p[t] = p[t-1] + ϕ[t]

    x[0] = 2.0
    p[0] = 1.5
    ω[t] ~ {0, 0.05, 0.10, ..., 0.45, 0.5} with uniform probability.
    ϕ[t] ~ {-0.25, -0.125, 0.125, 0.25} with uniform probability.
=#

using Kokako, GLPK, Statistics, Test

function joint_distribution(; kwargs...)
    names = tuple([first(kw) for kw in kwargs]...)
    values = tuple([last(kw) for kw in kwargs]...)
    output_type = NamedTuple{names, Tuple{eltype.(values)...}}
    distribution = map(output_type, Base.product(values...))
    return distribution[:]
end

function newsvendor_example()
    model = Kokako.PolicyGraph(
            Kokako.LinearGraph(3),
            sense = :Max,
            bellman_function = Kokako.AverageCut(upper_bound = 50.0),
            optimizer = with_optimizer(GLPK.Optimizer)
            ) do subproblem, stage
        @variables(subproblem, begin
            x >= 0, (Kokako.State, initial_value = 2)
            0 <= u <= 1
            w
        end)
        @constraint(subproblem, x.out == x.in - u + w)
        Kokako.add_objective_state(subproblem, initial_value = 1.5,
            lower_bound = 0.75, upper_bound = 2.25, lipschitz = 100.0) do y, ω
                return (y[1] + ω.price_noise,)
        end
        noise_terms = joint_distribution(
            demand = 0:0.05:0.5, price_noise = [-0.25, -0.125, 0.125, 0.25])
        Kokako.parameterize(subproblem, noise_terms) do ω
            JuMP.fix(w, ω.demand)
            price = Kokako.objective_state(subproblem)
            @stageobjective(subproblem, price * u)
        end
    end
    return model
end

model = newsvendor_example()
Kokako.train(model, iteration_limit = 100, print_level = 0)

results = Kokako.simulate(model, 500)
objectives = [
    sum(s[:stage_objective] for s in simulation) for simulation in results
]
sample_mean = round(Statistics.mean(objectives); digits = 2)
sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
println("Confidence_interval = $(sample_mean) ± $(sample_ci)")
@test Kokako.calculate_bound(model) ≈ sample_mean atol = sample_ci
