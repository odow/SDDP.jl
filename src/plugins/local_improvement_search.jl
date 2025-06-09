#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors and contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module LocalImprovementSearch

import JuMP

_norm(x) = sqrt(sum(xi^2 for xi in x))

abstract type AbstractSearchMethod end

"""
    minimize(
        f::Function,
        [method::AbstractSearchMethod = BFGS(100)],
        x₀::Vector{Float64},
        lower_bound::Float64 = -Inf,
    )

Minimizes a convex function `f` using first-order information.

Compared to off-the-shelf implementations, it has a number of features tailored
to this purpose:

 * Infeasibility is indicated by the function returning `nothing`. No other
   constraint information is given.
 * Sub-optimal solutions are okay, so we should focus on improving the feasible
   starting point, instead of finding the global minimizer.
 * `f` can be piecewise-linear convex with non-differentiable points.

## Arguments

 * `f(::Vector{Float64})`: takes a vector `x` and returns one of the following:
   * `nothing` if `x` is infeasible
   * `(f, Δf)::Tuple{Float64,Vector{Float64}`:  a tuple of the function
     evaluation and first-order gradient information.
 * `x₀::Vector{Float64}`: a feasible starting point.

## Default method

The default algorithm is a modified version of BFGS, with a specialized
back-tracking inexact line-search.
"""
function minimize end

function minimize(f::Function, x₀::Vector{Float64}, lower_bound::Float64 = -Inf)
    return minimize(f, BFGS(100), x₀, lower_bound)
end

###
### BFGS
###

struct BFGS <: AbstractSearchMethod
    evaluation_limit::Int
end

function minimize(
    f::F,
    bfgs::BFGS,
    x₀::Vector{Float64},
    lower_bound::Float64 = -Inf,
) where {F<:Function}
    # Initial estimte for the Hessian matrix in BFGS
    B = zeros(length(x₀), length(x₀))
    for i in 1:size(B, 1)
        B[i, i] = 1.0
    end
    # We assume that the initial iterate is feasible
    xₖ = x₀
    fₖ, ∇fₖ = f(xₖ)::Tuple{Float64,Vector{Float64}}
    @show fₖ, ∇fₖ
    # Initial step-length
    αₖ = 1.0
    # Evaluation counter
    evals = Ref(bfgs.evaluation_limit)
    for _ in 1:bfgs.evaluation_limit
        # Search direction. We could be clever here and maintain B⁻¹, but we're
        # only ever going to be solving this for very small |x| << 100 problems,
        # so taking the linear solve every time is okay. (The MIP solve is much
        # more of a bottleneck.)
        pₖ = B \ -∇fₖ
        @show pₖ
        # Run line search in direction `pₖ`
        @show αₖ, fₖ₊₁, ∇fₖ₊₁ = _line_search(f, fₖ, ∇fₖ, xₖ, pₖ, αₖ, evals)
        @show _norm(αₖ * pₖ) / max(1.0, _norm(xₖ)) < 1e-3
        if _norm(αₖ * pₖ) / max(1.0, _norm(xₖ)) < 1e-3
            # Small steps! Probably at the edge of the feasible region.
            # Return the current iterate.
            #
            # Note that "1e-3" isn't thaaaat small. But we hit a very annoying
            # feature of solvers: their feasibility checks are only approximate.
            # This tolerance is needed to pass the `test_kelleys_ip_xxx` tests.
            # Decreasing the tolerance leads to a _worse_ estimate for the dual,
            # because we abuse the solvers feasibility tolerance, and end up
            # returning a solution that is on the edge of numerical dual
            # feasibility.
            return fₖ, xₖ
        elseif _norm(∇fₖ₊₁) < 1e-6
            # Zero(ish) gradient. Return what must be a local maxima.
            return fₖ₊₁, xₖ + αₖ * pₖ
        elseif evals[] <= 0
            # We have evaluated the function too many times. Return our current
            # best.
            return fₖ₊₁, xₖ + αₖ * pₖ
        end
        # BFGS update.
        sₖ = αₖ * pₖ
        yₖ = ∇fₖ₊₁ - ∇fₖ
        # A slight tweak to normal BFGS: because we're dealing with non-smooth
        # problems, ||yₖ|| might be 0.0, i.e., we just moved along a facet from
        # from an interior point to a vertex, so the gradient stays the same.
        if _norm(yₖ) > 1e-12
            B .=
                B .+ (yₖ * yₖ') / (yₖ' * sₖ) -
                (B * sₖ * sₖ' * B') / (sₖ' * B * sₖ)
        end
        fₖ, ∇fₖ, xₖ = fₖ₊₁, ∇fₖ₊₁, xₖ + sₖ
    end
end

function _line_search(
    f::F,
    fₖ::Float64,
    ∇fₖ::Vector{Float64},
    x::Vector{Float64},
    p::Vector{Float64},
    α::Float64,
    evals::Ref{Int},
) where {F<:Function}
    while _norm(α * p) > 1e-3 * max(1.0, _norm(x))
        @show α
        @show p
        xₖ = x + α * p
        ret = f(xₖ)
        @show ret
        evals[] -= 1
        if ret === nothing  # Infeasible. So take a smaller step
            α /= 2
            continue
        end
        fₖ₊₁, ∇fₖ₊₁ = ret
        if p' * ∇fₖ₊₁ < 1e-6
            # Still a descent direction, so take a step.
            return α, fₖ₊₁, ∇fₖ₊₁
        elseif isapprox(fₖ + α * p' * ∇fₖ, fₖ₊₁; atol = 1e-8)
            # Step is onto a kink
            return α, fₖ₊₁, ∇fₖ₊₁
        elseif evals[] <= 0
            # Too many iterations
            return α, fₖ₊₁, ∇fₖ₊₁
        end
        #  Step is an ascent, so use Newton's method to find the intersection
        α = (fₖ₊₁ - fₖ - p' * ∇fₖ₊₁ * α) / (p' * ∇fₖ - p' * ∇fₖ₊₁)
    end
    return 0.0, fₖ, ∇fₖ
end

###
### Cutting plane
###

struct OuterApproximation{O} <: AbstractSearchMethod
    optimizer::O
end

function minimize(
    f::F,
    method::OuterApproximation,
    x₀::Vector{Float64},
    lower_bound::Float64,
) where {F<:Function}
    model = JuMP.Model(method.optimizer)
    JuMP.set_silent(model)
    n = length(x₀)
    JuMP.@variable(model, x[i in 1:n], start = x₀[i])
    JuMP.@variable(model, θ >= lower_bound)
    JuMP.@objective(model, Min, θ)
    xₖ = x₀
    fₖ, ∇fₖ = f(xₖ)::Tuple{Float64,Vector{Float64}}
    upper_bound = fₖ
    JuMP.@constraint(model, θ >= fₖ + ∇fₖ' * (x - xₖ))
    evals = Ref(0)
    d_step = Inf
    while d_step > 1e-8 && evals[] < 20
        JuMP.optimize!(model)
        lower_bound, xₖ₊₁ = JuMP.value(θ), JuMP.value.(x)
        ret = f(xₖ₊₁)
        while ret === nothing
            # point is infeasible
            xₖ₊₁ = 0.5 * (xₖ + xₖ₊₁)
            ret = f(xₖ₊₁)
        end
        fₖ₊₁, ∇fₖ₊₁ = ret::Tuple{Float64,Vector{Float64}}
        evals[] += 1
        upper_bound = fₖ₊₁
        JuMP.@constraint(model, θ >= fₖ₊₁ + ∇fₖ₊₁' * (x - xₖ₊₁))
        d = xₖ₊₁ - xₖ
        d_step = _norm(d)
        if sign(∇fₖ' * d) != sign(∇fₖ₊₁' * d)
            # There is a kink between the x
            xₖ₊₂ = 0.5 * (xₖ + xₖ₊₁)
            fₖ₊₂, ∇fₖ₊₂ = f(xₖ₊₂)::Tuple{Float64,Vector{Float64}}
            evals[] += 1
            JuMP.@constraint(model, θ >= fₖ₊₂ + ∇fₖ₊₂' * (x - xₖ₊₂))
            fₖ, xₖ = fₖ₊₂, xₖ₊₂
        else
            fₖ, xₖ = fₖ₊₁, xₖ₊₁
        end
    end
    return fₖ, xₖ
end

end
