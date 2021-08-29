#  Copyright 2017-21, Oscar Dowson.                                     #src
#  This Source Code Form is subject to the terms of the Mozilla Public  #src
#  License, v. 2.0. If a copy of the MPL was not distributed with this  #src
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.             #src

# # Biobjective hydro-thermal

using SDDP, GLPK, Statistics, Test

model = SDDP.LinearPolicyGraph(
    stages = 3,
    lower_bound = 0.0,
    optimizer = GLPK.Optimizer,
) do subproblem, stage
    @variable(subproblem, 0 <= v <= 200, SDDP.State, initial_value = 50)
    SDDP.add_objective_state(
        subproblem,
        initial_value = 0.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
        lipschitz = 1e6,
    ) do y, _
        return y
    end
    @variables(subproblem, begin
        0 <= g[i = 1:2] <= 100
        0 <= u <= 150
        s >= 0
        shortage_cost >= 0
    end)
    @expressions(subproblem, begin
        objective_1, g[1] + 10 * g[2]
        objective_2, shortage_cost
    end)
    @constraints(subproblem, begin
        inflow_constraint, v.out == v.in - u - s
        g[1] + g[2] + u == 150
        shortage_cost >= 40 - v.out
        shortage_cost >= 60 - 2 * v.out
        shortage_cost >= 80 - 4 * v.out
    end)
    Ω = 0.0:5:50.0
    SDDP.parameterize(subproblem, Ω) do ω
        JuMP.set_normalized_rhs(inflow_constraint, ω)
        ## This *has* to be called from inside `SDDP.parameterize`,
        ## otherwise it doesn't make sense.
        λ = SDDP.objective_state(subproblem)
        @stageobjective(subproblem, λ * objective_1 + (1 - λ) * objective_2)
    end
end

function set_trade_off_weight(model::SDDP.PolicyGraph, weight::Float64)
    @assert 0 <= weight <= 1
    for (_, node) in model.nodes
        node.objective_state.initial_value = (weight,)
    end
    return
end

function train_biobjective(
    model::SDDP.PolicyGraph;
    solution_limit::Int,
    kwargs...,
)
    solutions = Dict{Float64,Float64}()
    for weight in (0.0, 1.0)
        set_trade_off_weight(model, weight)
        SDDP.train(model; add_to_existing_cuts = true, kwargs...)
        solutions[weight] = SDDP.calculate_bound(model)
    end
    queue = Tuple{Float64,Float64}[(0.0, 1.0)]
    while length(queue) > 0 && length(solutions) < solution_limit
        (a, b) = popfirst!(queue)
        w = 0.5 * (a + b)
        set_trade_off_weight(model, w)
        SDDP.train(model; add_to_existing_cuts = true, kwargs...)
        bound = SDDP.calculate_bound(model)
        if !isapprox(0.5 * (solutions[a] + solutions[b]), bound; rtol=1e-4)
            push!(queue, (a, w))
            push!(queue, (w, b))
            solutions[w] = bound
        end
    end
    return solutions
end

pareto_weights = train_biobjective(
    model,
    solution_limit = 10,
    iteration_limit = 10,
)

@show pareto_weights
