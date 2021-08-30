#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module BiObjectiveSDDP

# Include the Gurobi-specific versions of get_BinvA and get_basis.
include(joinpath(@__DIR__, "gurobi.jl"))

import Printf
import SDDP

const MOI = SDDP.MOI
const JuMP = SDDP.JuMP

### Utilities for bi-objective simplex problems. Assumes that the problem is in
### the standard-form:
###
###     min c'x
###     s.t Ax = b
###          x ≥ 0
###
### This code is based on the description of the bi-objective simplex method
### given by M. Ehrgott in his slides located at:
### https://www.lamsade.dauphine.fr/~projet_cost/ALGORITHMIC_DECISION_THEORY/pdf/Ehrgott/HanLecture2_ME.pdf

"""
    get_BinvA(model::MOI.ModelLike)

Return the matrix `B⁻¹A`, where `B` is the matrix formed by the columns of the
basic variables in the constraint matrix.

Note that this typically has `n + m` columns, where `n` is the number of
variables and `m` is the number of rows in the constraint matrix because Gurobi
adds slack variables to equality constraints constraints for some reason.
"""
function get_BinvA end

"""
    get_basis(model::MOI.ModelLike)

Return the 1-indexed columns that comprise the basis.
"""
function get_basis end

# Here is where we kick things off properly.

struct BiObjectiveModel{M<:MOI.ModelLike}
    model::M
    c_1::Vector{Float64}
    c_2::Vector{Float64}
end

function set_weighted_objective(bi_obj_model::BiObjectiveModel, lambda::Float64)
    model = bi_obj_model.model
    x = MOI.get(model, MOI.ListOfVariableIndices())
    c = lambda .* bi_obj_model.c_1 + (1 - lambda) .* bi_obj_model.c_2
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0),
    )
    return
end

"""
    phase_iii_biobjective_simplex_step(
        bi_obj_model::BiObjectiveModel,
        current_lambda::Float64,
    )

Perform a Phase III step of the bi-objective simplex method.
"""
function phase_iii_biobjective_simplex_step(
    bi_obj_model::BiObjectiveModel,
    current_lambda::Float64,
)
    set_weighted_objective(bi_obj_model, current_lambda)
    model = bi_obj_model.model
    MOI.optimize!(model)
    # Get the matrix B⁻¹A and sanity checks.
    BinvA = get_BinvA(model)
    m, np = size(BinvA)
    n = np - m
    # Get the list of basis and non-basic variables, ignoring the slack
    # variables.
    B = get_basis(model)
    N = setdiff(1:n, B)
    # Now compute the c_N - c_B' B⁻¹A, but expand the cost vectors with `0.0` to
    # account for the slack variables introduced by Gurobi. These are the
    # reduced costs w.r.t. each objective, and are the change in the objective
    # function that would occur if we were to bring the corresponding variable
    # into the basis.
    c_1 = copy(bi_obj_model.c_1)
    c_2 = copy(bi_obj_model.c_2)
    c_n = length(c_1)
    if c_n < np
        resize!(c_1, np)
        resize!(c_2, np)
        c_1[(c_n+1):end] .= 0.0
        c_2[(c_n+1):end] .= 0.0
    end
    cb1 = c_1' .- c_1[B]' * BinvA
    cb2 = c_2' .- c_2[B]' * BinvA
    @assert length(cb1) == length(c_1)
    @assert length(cb2) == length(c_2)
    # Wondering from where the next peice of the formula comes? Here is a
    # toy example to explain it.
    #
    #   f: min  1x + 0y + 0z
    #   g: min  0x + 1y + 0z
    #      s.t. 2x + 1y - 1z == 1
    #            x,   y,   z >= 0
    #
    #      Decision space          Objective space
    #  y                        g
    #  | |                      | |
    # 1+ x A:(0, 1, 0)         1+ x A:(0, 1)
    #  |   \                    |    \
    #  |     \ B:(0.5, 0, 0)    |       \ B:(0.5, 0)
    # 0+       x------         0+         x--
    #  |-+―――――+――――――+― x      |-+-------+- f
    #    0    0.5     1           0      0.5
    #
    # To move from point A to point B, we have the change in objectives of
    #
    #   df/dΔ = +0.5
    #   dg/dΔ = -1.0
    #
    # The objective facet is defined by the normal vector
    #
    #   u = [-dy/dΔ / (dx/dΔ - dy/dΔ), 1 + dy/dΔ / (dx/dΔ - dy/dΔ)]
    #
    # or, more concretely, u = [2 / 3, 1 / 3].
    #
    # Therefore, the first point is valid for λ ∈ (2 / 3, 1] and the second
    # for λ ∈ [0, 2 / 3).

    # Set the initial bounds [λ⁻, λ⁺] over which the basis remains optimal.
    λ⁻, λ⁺ = -1.0, 2.0
    for i in N
        tmp = -cb2[i] / (cb1[i] - cb2[i])
        if cb1[i] >= 0 && cb2[i] < 0
            # An entering variable that will decrease the second objective and
            # increase the first objective. Lift the lower limit of the bound.
            λ⁻ = max(λ⁻, tmp)
        elseif cb1[i] < 0 && cb2[i] >= 0
            # An entering variable that will decrease the first objective and
            # increase the second objective. Drop the upper limit of the bound.
            λ⁺ = min(λ⁺, tmp)
        end
    end
    if λ⁺ == 2.0
        # We couldn't find an entering variable that will decrease the first
        # objective while increasing the second objective.
        λ⁺ = current_lambda
    end
    if λ⁻ == -1.0
        # We couldn't find an entering variable that will decrease the second
        # objective while increasing the first objective.
        λ⁻ = 0.0
    end
    return λ⁻, λ⁺
end

function get_next_lambda(
    bi_obj_model::BiObjectiveModel,
    current_lambda::Float64;
    lambda_minimum_step::Float64,
    lambda_atol::Float64,
)
    if current_lambda == 0.0
        return 1.0
    elseif current_lambda < lambda_minimum_step
        return 0.0
    end
    λ⁻, _ = phase_iii_biobjective_simplex_step(bi_obj_model, current_lambda)
    if λ⁻ < current_lambda - lambda_minimum_step
        # The Phase III step found a new value for lambda.
        return λ⁻
    end
    # Okay, so we didn't find a new value, we have a  problem. The basis we
    # currently have after the solve is valid for all λ ∈ [λ⁻, lambda]. So,
    # if set the weight to λ⁻, then we're going to get the same basis :(.
    # Instead of manually performing a simplex step, or changing the basis,
    # we're going to do a bisection search in λ to find the next one. It's a
    # little slower, because there might be multiple solves, but it's
    # simpler for the solvers to implement.
    #
    # TODO(odow): improve this.
    candidate = 0.0
    while candidate < current_lambda - lambda_minimum_step
        lambda_bar = (candidate + current_lambda) / 2
        λ⁻, λ⁺ = phase_iii_biobjective_simplex_step(bi_obj_model, lambda_bar)
        if isapprox(λ⁺, current_lambda, atol = lambda_atol)
            # Success! We found a new lambda that is provably the next one
            # in the sequence. Take a step of at least lambda_minimum_step
            # unless we hit λ = 0.
            candidate = λ⁻
            break
        else
            candidate = λ⁺
        end
        if isapprox(candidate, 0.0, atol=lambda_atol)
            break
        end
    end
    # Before we leave, we need to reset the model to the state we found it in.
    set_weighted_objective(bi_obj_model, current_lambda)
    MOI.optimize!(bi_obj_model.model)
    return min(candidate, max(current_lambda - lambda_minimum_step, 0))
end

###
### Utilities for converting an arbitrary MOI model into the standard form:
###
###     min c'x
###     s.t Ax = b
###          x ≥ 0
###

# Here's a type that is used to bridge into Ax = b; x ∈ R₊.
# We can't do SingleVariable-in-GreaterThan because it won't added the required
# slacks, leaving us with x >= l.

MOI.Utilities.@model(
    VectorStandardForm,
    (),
    (MOI.EqualTo,),
    (MOI.Nonnegatives,),
    (),
    (),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    ()
)

# We dis-allow free variables to bridge x free into x = y⁺ - y⁻; y⁺, y⁻ >= 0.

function MOI.supports_constraint(
    ::VectorStandardForm,
    ::Type{MOI.VectorOfVariables},
    ::Type{MOI.Reals},
)
    return false
end

# We dis-allow SingleVariable-in-S constraints to force bridging into
# VectorOfVariables-in-Nonnegatives.

function MOI.supports_constraint(
    ::VectorStandardForm{Float64},
    ::Type{MOI.SingleVariable},
    ::Type{S},
) where {
    S<:Union{
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.Interval{Float64},
        MOI.EqualTo{Float64},
    },
}
    return false
end

function vec_bridged_terms(x::MOI.VariableIndex, bridged, index_map)
    variable = MOI.Bridges.bridged_variable_function(bridged, index_map[x])
    if typeof(variable) == MOI.SingleVariable
        return [MOI.ScalarAffineTerm(1.0, variable.variable)]
    end
    @assert typeof(variable) <: MOI.ScalarAffineFunction
    return variable.terms
end

function std_bridged_terms(
    terms::Vector{MOI.ScalarAffineTerm{Float64}},
    bridged,
    index_map,
)
    terms_out = MOI.ScalarAffineTerm{Float64}[]
    for term in terms
        var = index_map[term.variable_index]
        for new_term in vec_bridged_terms(var, bridged, index_map)
            push!(
                terms_out,
                MOI.ScalarAffineTerm(
                    term.coefficient * new_term.coefficient,
                    new_term.variable_index,
                ),
            )
        end
    end
    return terms_out
end

"""
    convert_to_standard_form(
        dest::MOI.ModelLike,
        src::MOI.ModelLike,
    )::Dict{MOI.VariableIndex,Vector{MOI.ScalarAffineTerm{Float64}}}

Convert `src` to an equivalent model in `dest`, where `dest` has the form
`min{c'x: Ax = b, x >= 0}`.

Return a dictionary that maps variables from `src`-space, to an equivalent
vector of `ScalarAffineTerm`'s (if they were summed) in `dest`-space.
"""
function convert_to_standard_form(dest::MOI.ModelLike, src::MOI.ModelLike)
    x = MOI.get(src, MOI.ListOfVariableIndices())
    # First step: convert model in to `min{c'x: a'x = b, x ∈ R₊}`.
    vec_std = VectorStandardForm{Float64}()
    vec_bridge = MOI.Bridges.full_bridge_optimizer(vec_std, Float64)
    vec_map = MOI.copy_to(vec_bridge, src)
    vec_terms =
        Dict(xi => vec_bridged_terms(xi, vec_bridge, vec_map) for xi in x)
    # Second step: shift constants to RHS.
    for index in MOI.get(
        vec_std,
        MOI.ListOfConstraintIndices{
            MOI.ScalarAffineFunction{Float64},
            MOI.EqualTo{Float64},
        }(),
    )
        constant = MOI.get(vec_std, MOI.ConstraintFunction(), index).constant
        if !iszero(constant)
            set = MOI.get(vec_std, MOI.ConstraintSet(), index)
            MOI.set(
                vec_std,
                MOI.ConstraintSet(),
                index,
                MOI.EqualTo(set.value - constant),
            )
            MOI.modify(vec_std, index, MOI.ScalarConstantChange(0.0))
        end
    end
    # Third step: convert model in to `min{c'x: a'x = b, xᵢ >= 0}`.
    dest_bridge = MOI.Bridges.full_bridge_optimizer(dest, Float64)
    dest_map = MOI.copy_to(dest_bridge, vec_std)
    # Fourth step: reconcile the variables.
    return Dict(
        xi => std_bridged_terms(terms, dest_bridge, dest_map) for
        (xi, terms) in vec_terms
    )
end

###
### Utilities for working with bi-objective SDDP problems.
###

"""
    get_next_lambda(
        model::SDDP.PolicyGraph{T},
        node_index::T,
        noise,
        λ::Float64,
        dest::MOI.ModelLike,
    ) where {T}

"""
function get_next_lambda(
    model::SDDP.PolicyGraph{T},
    node_index::T,
    noise::Any,
    λ::Float64,
    dest::MOI.ModelLike;
    lambda_minimum_step::Float64,
    lambda_atol::Float64,
) where {T}
    # Look-up the node that we want to compute the step at.
    node = model[node_index]
    # Convert `node.subproblem` into the standard form.
    dest_variables =
        convert_to_standard_form(dest, JuMP.backend(node.subproblem))
    x = MOI.get(dest, MOI.ListOfVariableIndices())
    function compute_objective_vector(λ)
        # Set the objective and noise inn `node.subproblem`.
        SDDP.set_trade_off_weight(model, λ)
        SDDP.parameterize(node, noise)
        # Initialize the objective vector `c`.
        c = fill(0.0, length(x))
        objective = JuMP.objective_function(node.subproblem)
        for (variable, coef) in objective.terms
            src_index = JuMP.index(variable)
            for term in dest_variables[src_index]
                # TODO: here, we assume that .value is the column! This is
                # true in Gurobi, but it might not be the case for all solvers.
                c[term.variable_index.value] = coef * term.coefficient
            end
        end
        return c
    end
    c_1 = compute_objective_vector(1.0)
    c_2 = compute_objective_vector(0.0)
    # Quickly optimize `dest` to obtain a basis. Note: for sanity sake, the
    # ObjectiveValue of `dest` after this should be identical to the objective
    # value of node.subproblem (modulo numerical tolerances).
    MOI.set(dest, MOI.Silent(), true)
    MOI.optimize!(dest)
    # Get the next lambda using MOI calls defined above.
    return get_next_lambda(
        BiObjectiveModel(dest, c_1, c_2),
        λ;
        lambda_minimum_step = lambda_minimum_step,
        lambda_atol = lambda_atol,
    )
end

"""
    surrogate_lower_bound(
        model::SDDP.PolicyGraph,
        optimizer;
        global_lower_bound::Real,
        lambda_minimum_step::Float64 = 1e-4,
        lambda_atol::Float64 = 1e-4,
    )

Compute the surrogate lower bound for `model` using `optimizer` to calculate
the weight update.

`global_lower_bound` must be a valid lower bound across the entire weight-space.
"""
function surrogate_lower_bound(
    model::SDDP.PolicyGraph,
    optimizer;
    global_lower_bound::Real,
    lambda_minimum_step::Float64 = 1e-4,
    lambda_atol::Float64 = 1e-4,
)
    key = model.root_children[1].term
    node = model[key]
    noise = node.noise_terms[1].term
    if (length(model.root_children) != 1) || (length(node.noise_terms) != 1)
        error("Need deterministic first-stage")
    end
    weights, bounds = Float64[], Float64[]
    λ = 1.0
    while true
        SDDP.set_trade_off_weight(model, λ)
        bound = SDDP.calculate_bound(model)
        push!(weights, λ)
        push!(bounds, bound)
        λ = get_next_lambda(
            model,
            key,
            noise,
            λ,
            optimizer();
            lambda_minimum_step = lambda_minimum_step,
            lambda_atol = lambda_atol,
        )
        if λ == 1.0
            break
        end
    end
    lower_bound = 0.0
    for i in 2:length(weights)
        d_bound = 0.5 * (bounds[i-1] + bounds[i]) - global_lower_bound
        lower_bound += d_bound * (weights[i-1] - weights[i])
    end
    return lower_bound, weights, bounds
end

function print_iteration_header(io)
    println(
        io,
        "------------------------------------------------------------------",
    )
    println(
        io,
        "          BI-OBJECTIVE SDDP.jl (c) Oscar Dowson, 2019-21          ",
    )
    println(io)
    println(io, "      Iterations")
    println(
        io,
        " Maj.    Min.    SDDP    Lower Bound       Weight       Time (s)  ",
    )
    println(
        io,
        "------------------------------------------------------------------",
    )
    return
end

function print_iteration(
    io::IO,
    major_iterations::Int,
    minor_iterations::Int,
    sddp_iterations::Int,
    lower::Float64,
    weight::Float64,
    time::Float64,
)
    println(
        io,
        print_value(major_iterations),
        "  ",
        print_value(minor_iterations),
        "  ",
        print_value(sddp_iterations),
        "  ",
        print_value(lower),
        "  ",
        print_value(weight),
        " ",
        print_value(time),
    )
    return
end

function print_footer(io::IO)
    println(
        io,
        "------------------------------------------------------------------",
    )
    return
end

print_value(x::Real) = lpad(Printf.@sprintf("%1.6e", x), 13)

print_value(x::Int) = lpad(Printf.@sprintf("%5d", x), 6)

abstract type AbstractLambdaUpdate end

struct RandomUpdate <: AbstractLambdaUpdate end

function lambda_update(
    ::RandomUpdate,
    model::SDDP.PolicyGraph,
    λ::Float64,
    optimizer;
    kwargs...,
)
    key, node = rand(model.nodes)
    noise = rand(node.noise_terms)
    return get_next_lambda(model, key, noise.term, λ, optimizer(); kwargs...)
end

struct MinimumUpdate <: AbstractLambdaUpdate end

function lambda_update(
    ::MinimumUpdate,
    model::SDDP.PolicyGraph,
    λ::Float64,
    optimizer;
    kwargs...,
)
    if length(model.nodes) != 2
        error("Minimum update only implemented for two-stage problems.")
    end
    weights = Float64[]
    for (key, node) in model.nodes
        if length(node.noise_terms) != 1
            error("Minimum update only implemented for deterministic problems.")
        end
        noise = node.noise_terms[1]
        push!(
            weights,
            get_next_lambda(model, key, noise.term, λ, optimizer(); kwargs...),
        )
    end
    return maximum(weights)
end

struct FirstStageUpdate <: AbstractLambdaUpdate end

function lambda_update(
    ::FirstStageUpdate,
    model::SDDP.PolicyGraph,
    λ::Float64,
    optimizer;
    kwargs...,
)
    key = model.root_children[1].term
    node = model[key]
    noise = node.noise_terms[1]
    return get_next_lambda(model, key, noise.term, λ, optimizer(); kwargs...)
end

function bi_objective_sddp(
    model::SDDP.PolicyGraph,
    optimizer;
    bi_objective_major_iteration_limit::Int = 100,
    bi_objective_minor_iteration_limit::Int = 1_000,
    bi_objective_sddp_iteration_limit::Int = 10_000,
    bi_objective_lower_bound::Float64,
    bi_objective_lambda_atol::Float64 = 1e-4,
    bi_objective_major_iteration_burn_in::Int = 10,
    bi_objective_lambda_update_method::AbstractLambdaUpdate = RandomUpdate(),
    bi_objective_post_train_callback::Union{Function,Nothing} = nothing,
    kwargs...,
)
    print_iteration_header(stdout)
    start_time = time()
    major_iterations = 0
    minor_iterations = 0
    sddp_iterations = 0
    λ = 1.0
    try
        while true
            if λ == 1.0
                major_iterations += 1
            end
            SDDP.set_trade_off_weight(model, λ)
            SDDP.train(
                model;
                run_numerical_stability_report = false,
                add_to_existing_cuts = true,
                kwargs...,
            )
            minor_iterations += 1
            sddp_iterations += length(model.most_recent_training_results.log)
            if bi_objective_post_train_callback !== nothing
                bi_objective_post_train_callback(model, λ)
            end
            lower_bound, _, _ = surrogate_lower_bound(
                model,
                optimizer;
                global_lower_bound = bi_objective_lower_bound,
                lambda_minimum_step = bi_objective_lambda_atol,
                lambda_atol = bi_objective_lambda_atol,
            )
            tmp_bi_objective_lambda_tol =
                max(
                    1,
                    bi_objective_major_iteration_burn_in - major_iterations,
                ) * (0.1 - bi_objective_lambda_atol) /
                bi_objective_major_iteration_burn_in
            λ′ = lambda_update(
                bi_objective_lambda_update_method,
                model,
                λ,
                optimizer;
                lambda_minimum_step = tmp_bi_objective_lambda_tol,
                lambda_atol = tmp_bi_objective_lambda_tol,
            )
            # Clean up the iteration.
            print_iteration(
                stdout,
                major_iterations,
                minor_iterations,
                sddp_iterations,
                lower_bound,
                λ,
                time() - start_time,
            )
            λ = λ′
            if major_iterations >= bi_objective_major_iteration_limit
                break
            elseif minor_iterations >= bi_objective_minor_iteration_limit
                break
            elseif sddp_iterations >= bi_objective_sddp_iteration_limit
                break
            end
        end
    catch ex
        if isa(ex, InterruptException)
            println("Terminating: solve interrupted by user.")
        else
            rethrow(ex)
        end
    end
    print_footer(stdout)
    return surrogate_lower_bound(
        model,
        optimizer;
        global_lower_bound = bi_objective_lower_bound,
        lambda_minimum_step = bi_objective_lambda_atol,
        lambda_atol = bi_objective_lambda_atol,
    )
end

struct Hyperplane
    λ::Float64
    y::Float64
    Δ::Float64
end

# TODO(odow): this function is a pretty inefficient O(N²) operation. We could
# make this more efficient by being clever.
function compute_vertices(hyperplanes::Vector{Hyperplane}, atol::Float64)
    N = length(hyperplanes)
    vertices = Tuple{Float64,Float64}[(0.0, hyperplanes[1].y)]
    x_min, i_min = 0.0, 1
    while x_min < 1.0
        # Here is our current λ co-ordinate
        x = vertices[end][1]
        # and our current intercept and slope.
        y, Δ = hyperplanes[i_min].y, hyperplanes[i_min].Δ
        # Now loop through every hyperplane and find the minimum x intersection
        # greater than our current x.
        x_min, i_min = Inf, 0
        for i in 1:N
            # Get intersection with current plane and the i'th one.
            yᵢ, Δᵢ = hyperplanes[i].y, hyperplanes[i].Δ
            x′ = if isapprox(y, yᵢ, atol = atol) && isapprox(Δ, Δᵢ; atol = atol)
                1.0
            else
                (y - yᵢ) / (Δᵢ - Δ)
            end
            if x < x′ < x_min
                x_min, i_min = x′, i
            end
        end
        if x_min ≈ 1.0
            y, Δ = hyperplanes[end].y, hyperplanes[end].Δ
        end
        push!(vertices, (x_min, y + Δ * x_min))
    end
    return vertices
end

struct ObjectiveSpacePoint
    f₁::Float64
    f₂::Float64
    λ::Float64
end

V(p::ObjectiveSpacePoint, λ::Float64) = λ * p.f₁ + (1 - λ) * p.f₂

function compute_upper_bound(
    points::Vector{ObjectiveSpacePoint};
    atol::Float64,
    lower_bound::Float64,
)
    # Each tuple in hyperplanes is `(x, intercept, slope)`.
    hyperplanes =
        [Hyperplane(p.λ, V(p, 0.0), V(p, 1.0) - V(p, 0.0)) for p in points]
    # Drop colinear hyperplanes.
    unique!(h -> (h.y, h.Δ), hyperplanes)
    # Sort first by smallest intercept, then by largest slope, then by x.
    sort!(hyperplanes, by = h -> (h.y, -h.Δ, h.λ))
    # Now compute the vertices.
    vertices = compute_vertices(hyperplanes, atol)
    # Now compute the area contained in the polygon
    area = 0.0
    for i in 2:length(vertices)
        xl, yl = vertices[i-1]
        xr, yr = vertices[i]
        area += (xr - xl) * ((yl + yr) / 2 - lower_bound)
    end
    return area
end

function surrogate_upper_bound(
    model::SDDP.PolicyGraph,
    optimizer;
    global_lower_bound::Real,
    lambda_minimum_step::Float64 = 1e-4,
    lambda_atol::Float64 = 1e-4,
)
    if length(model.nodes) != 2
        error("Upper bound only defined for two-stage problems")
    elseif any(length(node.noise_terms) > 1 for (key, node) in model.nodes)
        error("Upper bound only defined for deterministic problems")
    end
    points = ObjectiveSpacePoint[]
    λ = 1.0
    while true
        SDDP.set_trade_off_weight(model, λ)
        simulation = SDDP.simulate(model, 1, [:objective_1, :objective_2])
        push!(
            points,
            ObjectiveSpacePoint(
                sum(s[:objective_1] for s in simulation[1]),
                sum(s[:objective_2] for s in simulation[1]),
                λ,
            ),
        )
        key = model.root_children[1].term
        node = model[key]
        noise = node.noise_terms[1].term
        λ = get_next_lambda(
            model,
            key,
            noise,
            λ,
            optimizer();
            lambda_minimum_step = lambda_minimum_step,
            lambda_atol = lambda_atol,
        )
        if λ == 1.0
            break
        end
    end
    return compute_upper_bound(
        points;
        atol = lambda_atol,
        lower_bound = global_lower_bound,
    )
end

end
