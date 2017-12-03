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

The distributionally robust SDDP risk measure.
Constructs a DRO risk measure object that allows probabilities to deviate by
`radius` away from the uniform distribution.
"""
type DRO <: SDDP.AbstractRiskMeasure
    radius::Float64
end

# Helper to calculate population standard deviation, avoid type instability
function popvar(x::Vector{Float64})::Float64
    ninv = 1 / length(x)
    ninv * sum(x.^2) - (ninv * sum(x))^2
end

function popstd(x::Vector{Float64})
    sqrt(popvar(x))
end

function is_dro_applicable(radius::Float64, observations::Vector{Float64})
    if (abs(radius) < 1e-9)
        return false
    elseif abs(popstd(observations)) < 1e-9
        return false
    end
    return true
end

function getconstfactor(S::Int, k::Int, radius::Float64, permuted_observations::Vector{Float64})
    stdz = popstd(permuted_observations[k+1:S])::Float64
    sqrt((S-k) * radius^2 - k/S) / (stdz * (S-k))
end

function getconstadditive(S::Int, k::Int, const_factor, permuted_observations)
    avgz = mean(permuted_observations[k+1:S])
    1 / (S-k) + const_factor * avgz
end

function modifyprobability!(newprobabilities, dro::DRO, sense::Symbol, observations::Vector{Float64}, S::Int)

    # Don't do any DRO reweighting if we aren't distributionally robust or the variance is too low
    r = dro.radius
    if !is_dro_applicable(r, observations)
        newprobabilities .= 1 / S
        return nothing
    end

    # Sort future costs/rewards
    if sense == :Min
        perm = sortperm(observations)
        permuted_observations = -observations[perm]
    else
        perm = sortperm(observations, rev=true)
        permuted_observations = observations[perm]
    end

    # Compute the new probabilities
    @inbounds for k = 0:S-2
        if k > 0
            newprobabilities[perm[k]] = 0.0
        end
        const_factor   = getconstfactor(S, k, r, permuted_observations)
        const_additive = getconstadditive(S, k, const_factor, permuted_observations)
        @inbounds for i = k+1:S
            newprobabilities[perm[i]] = const_additive - const_factor * permuted_observations[i]::Float64
        end
        if newprobabilities[perm[k+1]] >= 0.0
            break
        end
    end

end

function modifyprobability!(
        measure::DRO,
        newprobabilities,
        original_distribution::Vector{Float64},
        observations::Vector{Float64},
        m::SDDPModel,
        sp::JuMP.Model
    )

     # Number of noises
    S = length(observations)

    modifyprobability!(newprobabilities, measure, getsense(sp), observations, S)
end
