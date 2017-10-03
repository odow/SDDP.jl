#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#   Average Value at Risk
#   λ * E[x] + (1 - λ) * AV@R(1-β)[x]
immutable NestedAVaR <: AbstractRiskMeasure
    lambda::Float64
    beta::Float64
    function NestedAVaR(lambda::Float64, beta::Float64)
        if lambda > 1.0 || lambda < 0.0
            error("Lambda must be in the range [0, 1]. Increasing values of lambda are less risk averse. lambda=1 is identical to taking the expectation.")
        end
        if beta > 1.0 || beta < 0.0
            error("Beta must be in the range [0, 1]. Increasing values of beta are less risk averse. beta=1 is identical to taking the expectation.")
        end
        new(lambda, beta)
    end
end

"""
    NestedAVaR(;lambda=1.0, beta=1.0)

# Description

A risk measure that is a convex combination of Expectation and Average Value @ Risk
(also called Conditional Value @ Risk).

    λ * E[x] + (1 - λ) * AV@R(1-β)[x]

# Keyword Arguments

 * `lambda`
Convex weight on the expectation (`(1-lambda)` weight is put on the AV@R component.
Inreasing values of `lambda` are less risk averse (more weight on expecattion)

 * `beta`
 The quantile at which to calculate the Average Value @ Risk. Increasing values
 of `beta` are less risk averse. If `beta=0`, then the AV@R component is the
 worst case risk measure.

# Returns

    m::NestedAVaR<:AbstractRiskMeasure
"""
NestedAVaR(;lambda=1.0, beta=1.0) = NestedAVaR(lambda, beta)

function modifyprobability!(measure::NestedAVaR,
    riskadjusted_distribution,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    ismax = getsense(sp) == :Max
    # sanity
    @assert length(riskadjusted_distribution) == length(original_distribution) == length(observations)
    if measure.beta < 1e-8
        # worst case
        positive_prob = original_distribution.>0
        if ismax
            positive_idx = indmin(observations[positive_prob])
        else
            positive_idx = indmax(observations[positive_prob])
        end
        idx = (1:length(observations))[positive_prob][positive_idx]
        riskadjusted_distribution .= measure.lambda * original_distribution
        riskadjusted_distribution[idx] += (1 - measure.lambda)
        return
    end
    # Quantile collected so far
    q = 0.0
    # For each noise in order
    for i in sortperm(observations, rev=!ismax)
        if q >=  measure.beta
            # We have collected the beta quantile, therefore
            # AV@R contribution is zero
            avar_prob  = 0.0
        else
            # AV@R
            avar_prob  = min(original_distribution[i], measure.beta-q) / measure.beta
        end
        # take the biggest proportion of the noise possible
        riskadjusted_distribution[i] = measure.lambda * original_distribution[i] + (1-measure.lambda) * avar_prob
        # Update total quantile collected
        q += avar_prob * measure.beta
    end
end
