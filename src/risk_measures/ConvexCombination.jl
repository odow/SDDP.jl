#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

struct ConvexCombination <: AbstractRiskMeasure
    measures::Vector{Tuple{Float64, AVaR}}
end

function modifyprobability!(riskmeasure::ConvexCombination,
    riskadjusted_distribution,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    riskadjusted_distribution .= 0.0
    y = similar(original_distribution)
    for (weight, measure) in riskmeasure.measures
        y .= 0.0
        modifyprobability!(measure, y, original_distribution, observations, m, sp)
        riskadjusted_distribution .+= weight * y
    end
    return
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
function NestedAVaR(;lambda::Float64=1.0, beta::Float64=1.0)
    if lambda > 1.0 || lambda < 0.0
        error("Lambda must be in the range [0, 1]. Increasing values of lambda are less risk averse. lambda=1 is identical to taking the expectation.")
    end
    if beta > 1.0 || beta < 0.0
        error("Beta must be in the range [0, 1]. Increasing values of beta are less risk averse. beta=1 is identical to taking the expectation.")
    end
    ConvexCombination([(lambda, AVaR(1.0)), ((1-lambda), AVaR(beta))])
end
