#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ========================== The Expectation Operator ======================== #

struct Expectation <: AbstractRiskMeasure end

function adjust_probability(measure::Expectation,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector,
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool)
    risk_adjusted_probability .= original_probability
    return
end

# ========================== The Worst Case Operator ========================= #

struct WorstCase <: AbstractRiskMeasure end

function adjust_probability(measure::WorstCase,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector,
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool)
    risk_adjusted_probability .= 0.0
    worst_index = 1
    worst_observation = is_minimization ? -Inf : Inf
    for (index, (probability, observation)) in enumerate(
            zip(original_probability, objective_realizations))
        if probability > 0.0
            if (is_minimization && observation > worst_observation) ||
                    (!is_minimization && observation < worst_observation)
                worst_index = index
                worst_observation = observation
            end
        end
    end
    risk_adjusted_probability[worst_index] = 1.0
    return
end

# =================================== AV@R =================================== #

struct AVaR <: AbstractRiskMeasure
    β::Float64
    function AVaR(β::Float64)
        if !(0 <= β <= 1)
            throw(ArgumentError(
                "Risk-quantile β must be in [0, 1]. Currently it is $(β)."))
        end
        return new(β)
    end
end

function adjust_probability(measure::AVaR,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector,
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool)
     if measure.β ≈ 0.0
         return adjust_probability(
            WorstCase(), risk_adjusted_probability, original_probability,
            noise_support, objective_realizations, is_minimization)
     elseif measure.β ≈ 1.0
         return adjust_probability(
            Expectation(), risk_adjusted_probability, original_probability,
            noise_support, objective_realizations, is_minimization)
     end
     risk_adjusted_probability .= 0.0
     quantile_collected = 0.0
     for i in sortperm(objective_realizations, rev = is_minimization)
         quantile_collected >= measure.β && break
         avar_prob = min(
                original_probability[i],
                measure.β - quantile_collected
            ) / measure.β
         risk_adjusted_probability[i] = avar_prob
         quantile_collected += avar_prob * measure.β
    end
    return
end

# ============================ ConvexCombination ============================= #

"""
    ConvexCombination((weight::Float64, measure::AbstractRiskMeasure)...)

Create a weighted combination of risk measures.

### Examples

    Kokako.ConvexCombination(
        (0.5, Kokako.Expectation()),
        (0.5, Kokako.AVaR(0.25))
    )

Convex combinations can also be constructed by adding weighted risk measures
together as follows:

    0.5 * Kokako.Expectation() + 0.5 * Kokako.AVaR(0.5)
"""
struct ConvexCombination{T<:Tuple} <: AbstractRiskMeasure
    measures::T
end
ConvexCombination(args::Tuple...) = ConvexCombination(args)

import Base: +, *

function Base.:+(a::ConvexCombination, b::ConvexCombination)
    return ConvexCombination(a.measures..., b.measures...)
end

function Base.:*(lhs::Float64, rhs::AbstractRiskMeasure)
    return ConvexCombination(((lhs, rhs),))
end

function adjust_probability(measure::ConvexCombination,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector,
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool)
    risk_adjusted_probability .= 0.0
    partial_distribution = similar(risk_adjusted_probability)
    for (weight, measure) in measure.measures
        partial_distribution .= 0.0
        adjust_probability(
            measure, partial_distribution, original_probability, noise_support,
            objective_realizations, is_minimization)
        risk_adjusted_probability .+= weight * partial_distribution
    end
end

# =================================== EAV@R ================================== #

"""
    EAVaR(;lambda=1.0, beta=1.0)

A risk measure that is a convex combination of Expectation and Average Value @
Risk (also called Conditional Value @ Risk).

        λ * E[x] + (1 - λ) * AV@R(1-β)[x]

### Keyword Arguments

* `lambda`: Convex weight on the expectation (`(1-lambda)` weight is put on the
  AV@R component. Inreasing values of `lambda` are less risk averse (more
  weight on expectation).

* `beta`: The quantile at which to calculate the Average Value @ Risk.
  Increasing values of `beta` are less risk averse. If `beta=0`, then the AV@R
  component is the worst case risk measure.
"""
function EAVaR(;lambda::Float64=1.0, beta::Float64=1.0)
    if !(0.0 <= lambda <= 1.0)
        error("Lambda must be in the range [0, 1]. Increasing values of " *
              "lambda are less risk averse. lambda=1 is identical to taking " *
              "the expectation.")
    end
    if !(0.0 <= beta <= 1.0)
        error("Beta must be in the range [0, 1]. Increasing values of beta " *
              "are less risk averse. lambda=1 is identical to taking the " *
              "expectation.")
    end
    return lambda * Expectation() + (1 - lambda) * AVaR(beta)
end

# ===================================== DRO ================================== #

#=
This code was contributed by Lea Kapelevich.

In a Distributionally Robust Optimization (DRO) approach, we modify the
probabilities we associate with all future scenarios so that the resulting
probability distribution is the "worst case" probability distribution, in some
sense.

In each backward pass we will compute a worst case probability distribution
vector ̃p. We compute ̃p so that:

̄p ∈ argmax{̃pᵀ̃z}
    ||̃p - ̃a||₂ ≤ r
    ∑̃p = 1
    ̃p ≥ 0

where

 1. ̃z is a vector of future costs. We assume that our aim is to minimize
    future cost pᵀ̃z. If we maximize reward, we would have ̃p ∈ argmin{̃pᵀ̃z}.
2. ̄a is the uniform distribution
3. r is a user specified radius - the larger the radius, the more conservative
   the policy.

Note: the largest radius that will work with S scenarios is sqrt((S-1)/S).
=#

"""
    DRO(radius::Float64)

The distributionally robust SDDP risk measure of

Philpott, A., de Matos, V., Kapelevich, L. Distributionally robust SDDP.
Computational Management Science (2018) 165:431-454.
"""
struct DRO <: AbstractRiskMeasure
    radius::Float64
end

function adjust_probability(measure::DRO,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector,
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool)
    if abs(measure.radius) < 1e-9 ||
            Statistics.std(objective_realizations, corrected=false) < 1e-9
        # Don't do any DRO reweighting if the radius is small or the variance
        # is too low. Use Expectation instead.
        return adjust_probability(
           Expectation(), risk_adjusted_probability, original_probability,
           noise_support, objective_realizations, is_minimization)
   end
   m = length(objective_realizations)
   if all(x-> x ≈ 1 / m, original_probability)
       uniform_dro(measure, risk_adjusted_probability, original_probability,
                   objective_realizations, is_minimization)
   else
       non_uniform_dro(measure, risk_adjusted_probability, original_probability,
                       objective_realizations, is_minimization)
   end
   return
end

"""
Algorithm (1) of Philpott et al. Assumes that the nominal distribution is _not_
uniform.
"""
function non_uniform_dro(
        measure::DRO,
        risk_adjusted_probability::Vector{Float64},
        original_probability::Vector{Float64},
        objective_realizations::Vector{Float64},
        is_minimization::Bool)
    error("Current implementation of DRO assumes that the nominal " *
          "distribution is uniform.")
end

"""
Algorithm (2) of Philpott et al. Assumes that the nominal distribution is
uniform.
"""
function uniform_dro(
        measure::DRO,
        risk_adjusted_probability::Vector{Float64},
        original_probability::Vector{Float64},
        objective_realizations::Vector{Float64},
        is_minimization::Bool)
    m = length(objective_realizations)
    # Take a permuted view of `risk_adjusted_probability` so we can refer to
    # `p[i]` instead of `risk_adjusted_probability[perm[i]]`.
    perm = sortperm(objective_realizations, rev = !is_minimization)
    p = view(risk_adjusted_probability, perm)
    z = view(objective_realizations, perm)
    # Compute the new probabilities according to Algorithm (2) of the Philpott
    # et al. paper.
    # Step (1):
    @inbounds for k in 0:m-2
        # Step (1a):
        z_bar = sum(z[i] for i in (k+1):m) / (m - k)
        s = sqrt(sum(z[i]^2 - z_bar^2 for i in (k+1):m) / (m - k))
        # Step (1b): note that we cache a couple of terms that don't depend on i
        #            to speed things up.
        term_1 = 1 / (m - k)
        term_2 = sqrt((m - k) * measure.radius^2 - k / m) / ((m - k) * s)
        # We really should set p[i] = 0 for i = 1, ..., k. But since we don't
        # touch p[k-1] again, we can just set the k'th element to 0.
        if k > 0
            p[k] = 0.0
        end
        if is_minimization
            @inbounds for i in (k+1):m
                p[i] = term_1 + term_2 * (z[i] - z_bar)
            end
        else
            # Okay, here's the rub: we should have converted
            # objective_realizations (rewards) into costs by negating them. This
            # would have required a copy. This means that z_bar is in fact the
            # -ve of what it should be. `s` is fine since it is a difference of
            # squares. Thus, all we have to do is negate both z[i] and z_bar
            # here.
            @inbounds for i in (k+1):m
                p[i] = term_1 + term_2 * (z_bar - z[i])
            end
        end
        # Step (1c)
        if p[k+1] >= 0.0
            return
        end
    end
    # Step (2):
    p[end] = 1.0
    return
end

# ================================= Wasserstein ============================== #

"""
    Wasserstein(solver_factory; alpha::Float64)

A distributionally-robust risk measure based on the Wasserstein distance.

As `alpha` increases, the measure becomes more risk-averse. When `alpha=0`, the
measure is equivalent to the expectation operator. As `alpha` increases, the
measure approaches the Worst-case risk measure.
"""
struct Wasserstein{T, F} <: AbstractRiskMeasure
    alpha::Float64
    solver_factory::T
    norm::F
    function Wasserstein(norm::Function, solver_factory; alpha::Float64)
        if alpha < 0.0
            error("alpha cannot be $(alpha) as it must be in the range [0, ∞).")
        end
        return new{typeof(solver_factory), typeof(norm)}(alpha, solver_factory, norm)
    end
end

function adjust_probability(measure::Wasserstein,
                            risk_adjusted_probability::Vector{Float64},
                            original_probability::Vector{Float64},
                            noise_support::Vector,
                            objective_realizations::Vector{Float64},
                            is_minimization::Bool)
    N = length(objective_realizations)
    wasserstein = JuMP.Model(measure.solver_factory)
    @variable(wasserstein, z[1:N, 1:N] >= 0)
    @variable(wasserstein, p[1:N] >= 0)
    for i in 1:N
        @constraint(wasserstein, sum(z[:, i]) == original_probability[i])
        @constraint(wasserstein, sum(z[i, :]) == p[i])
    end
    @constraint(
        wasserstein,
        sum(measure.norm(noise_support[i], noise_support[j]) * z[i, j]
            for i in 1:N, j in 1:N) <= measure.alpha
    )
    objective_sense = is_minimization ? MOI.MAX_SENSE : MOI.MIN_SENSE
    @objective(
        wasserstein,
        objective_sense,
        sum(objective_realizations[i] * p[i] for i in 1:N)
    )
    JuMP.optimize!(wasserstein)
    if JuMP.primal_status(wasserstein) != MOI.FEASIBLE_POINT
        error("Unable to solver Wasserstein subproblem. Status: ", JuMP.termination_status(wassserstein))
    end
    copyto!(risk_adjusted_probability, JuMP.value.(p))
    return
end
