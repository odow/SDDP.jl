#  Copyright 2017-20, Oscar Dowson
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

using SDDP, GLPK, Statistics, Test

function joint_distribution(; kwargs...)
    names = tuple([first(kw) for kw in kwargs]...)
    values = tuple([last(kw) for kw in kwargs]...)
    output_type = NamedTuple{names,Tuple{eltype.(values)...}}
    distribution = map(output_type, Base.product(values...))
    return distribution[:]
end

function newsvendor_example(; cut_type)
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(3),
        sense = :Max,
        upper_bound = 50.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, stage
        @variables(subproblem, begin
            x >= 0, (SDDP.State, initial_value = 2)
            0 <= u <= 1
            w
        end)
        @constraint(subproblem, x.out == x.in - u + w)
        SDDP.add_objective_state(
            subproblem,
            initial_value = 1.5,
            lower_bound = 0.75,
            upper_bound = 2.25,
            lipschitz = 100.0,
        ) do y, ω
            return y + ω.price_noise
        end
        noise_terms = joint_distribution(
            demand = 0:0.05:0.5,
            price_noise = [-0.25, -0.125, 0.125, 0.25],
        )
        SDDP.parameterize(subproblem, noise_terms) do ω
            JuMP.fix(w, ω.demand)
            price = SDDP.objective_state(subproblem)
            @stageobjective(subproblem, price * u)
        end
    end
    SDDP.train(
        model,
        iteration_limit = 100,
        print_level = 0,
        time_limit = 20.0,
        cut_type = cut_type,
    )
    @test SDDP.calculate_bound(model) ≈ 4.04 atol = 0.05
    results = SDDP.simulate(model, 500)
    objectives = [sum(s[:stage_objective] for s in simulation) for simulation in results]
    @test round(Statistics.mean(objectives); digits = 2) ≈ 4.04 atol = 0.1
end

newsvendor_example(cut_type = SDDP.SINGLE_CUT)
newsvendor_example(cut_type = SDDP.MULTI_CUT)
