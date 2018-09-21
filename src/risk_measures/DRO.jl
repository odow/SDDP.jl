#  Copyright 2017, Oscar Dowson, Lea Kapelevich
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export DRO

#==
This file relates to distributionally robust SDDP, see:
    Philpott, A., de Matos, V., Kapelevich, L. (2017).
    Distributionally robust SDDP.
    Tech Report. Electric Power Optimisation Centre. http://www.epoc.org.nz.

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

Here
[1] ̃z is a vector of future costs. (We assumed that our aim is to minimize
future cost pᵀ̃z. If we were  maximizing reward, we would have ̃p ∈ argmin{̃pᵀ̃z}).
[2] ̄a is the uniform distribution, i.e. ones(number_of_scenarios) / number_of_scenarios,
[3] r is a user specified radius- the larger the radius, the more conservative the policy.

Note: the largest radius that will work with S scenarios is sqrt((S-1)/S).
==#

"""
    DRO(radius::Float64)

The distributionally robust SDDP risk measure of

Philpott, A., de Matos, V., Kapelevich, L. (2017). Distributionally robust SDDP.
Tech Report. Electric Power Optimisation Centre. http://www.epoc.org.nz.

Note: requires that the initial probability distribution is uniform.
"""
struct DRO <: SDDP.AbstractRiskMeasure
    radius::Float64
end

# Helper to calculate population standard deviation, avoiding the penalty of
# using a keyword argument.
function popvar(x::Vector{Float64})::Float64
    return sum(x.^2) / length(x) - (sum(x) / length(x))^2
end

function popstd(x::Vector{Float64})
    return sqrt(popvar(x))
end

function is_dro_applicable(radius::Float64, observations::Vector{Float64})
    return abs(radius) >= 1e-9 && abs(popstd(observations)) >= 1e-9
end

function getconstfactor(S::Int, k::Int, radius::Float64, permuted_observations::Vector{Float64})
    stdz = popstd(permuted_observations[k+1:S])::Float64
    return sqrt((S-k) * radius^2 - k/S) / (stdz * (S-k))
end

function getconstadditive(S::Int, k::Int, const_factor, permuted_observations)
    avgz = mean(permuted_observations[k+1:S])
    return 1 / (S-k) + const_factor * avgz
end

function modify_probability(measure::DRO,
                            riskadjusted_distribution,
                            original_distribution::Vector{Float64},
                            observations::Vector{Float64},
                            model::SDDPModel,
                            subproblem::JuMP.Model)
    S = length(observations)  # Number of noise terms.
    # Don't do any DRO reweighting if we aren't distributionally robust or the
    # variance is too low.
    if !is_dro_applicable(measure.radius, observations)
        riskadjusted_distribution .= 1 / S
        return
    end
    # Sort future costs/rewards
    if getsense(subproblem) == :Min
        perm = sortperm(observations)
        permuted_observations = -observations[perm]
    else
        perm = sortperm(observations, rev=true)
        permuted_observations = observations[perm]
    end
    # Compute the new probabilities
    @inbounds for k in 0:S-2
        if k > 0
            riskadjusted_distribution[perm[k]] = 0.0
        end
        const_factor = getconstfactor(S, k, measure.radius, permuted_observations)
        const_additive = getconstadditive(S, k, const_factor, permuted_observations)
        @inbounds for i in k+1:S
            riskadjusted_distribution[perm[i]] = const_additive - const_factor * permuted_observations[i]::Float64
        end
        if riskadjusted_distribution[perm[k+1]] >= 0.0
            break
        end
    end
end
