#  Copyright 2017-21, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    set_biobjective_functions(subproblem, objective_1, objective_2)

Se the biobjective functions in `subproblem`.

This must be called from inside [`SDDP.parameterize`](@ref).

We recommend you define both objectives as JuMP expressions.

See also: [`initialize_biobjective_subproblem`](@ref).

!!! warning
    This function is experimental! It may change in any future release.
"""
function set_biobjective_functions(subproblem, objective_1, objective_2)
    λ = SDDP.objective_state(subproblem)
    @stageobjective(subproblem, λ * objective_1 + (1 - λ) * objective_2)
    return
end

"""
    initialize_biobjective_subproblem(subproblem)

Run some initialization code to setup a biobjective problem.

This must be called outside [`SDDP.parameterize`](@ref).

!!! warning
    This function is experimental! It may change in any future release.
"""
function initialize_biobjective_subproblem(subproblem)
    SDDP.add_objective_state(
        subproblem,
        initial_value = 0.0,
        lower_bound = 0.0,
        upper_bound = 1.0,
        lipschitz = 1e6,
    ) do y, _
        return y
    end
    return
end

"""
    set_trade_off_weight(model::SDDP.PolicyGraph, weight::Float64)

Set the trade-off weight of a bi-objective problem to `weight`.

!!! warning
    This function is experimental! It may change in any future release.
"""
function set_trade_off_weight(model::SDDP.PolicyGraph, weight::Float64)
    @assert 0 <= weight <= 1
    for (_, node) in model.nodes
        node.objective_state.initial_value = (weight,)
    end
    return
end

"""
    train_biobjective(
        model::SDDP.PolicyGraph;
        solution_limit::Int,
        kwargs...,
    )

Train a biobjective problem using a variation of the non-inferior set estimation
method.

## arguments

 * `solution_limit` is the maximum number of unique policies to return.
 * `kwargs` are passed to [`SDDP.train`](@ref) when solving the scalarized
   problems.

## Returns

Returns a dictionary mapping trade-off weights to their scalarized objective
value.

!!! warning
    This function is experimental! It may change in any future release.
"""
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
        if !isapprox(0.5 * (solutions[a] + solutions[b]), bound; rtol = 1e-4)
            push!(queue, (a, w))
            push!(queue, (w, b))
            solutions[w] = bound
        end
    end
    return solutions
end
