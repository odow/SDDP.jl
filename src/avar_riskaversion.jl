#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#   Average Value at Risk
#   (1 - λ) * E[x] + λ * AV@R(1-β)[x]
export NestedAVaR

immutable NestedAVaR <: AbstractRiskMeasure
    lambda::Float64
    beta::Float64
    function NestedAVaR(lambda::Float64, beta::Float64)
        if lambda > 1.0 || lambda < 0.0
            error("Lambda must be in the range [0, 1]. Increasing values of lambda are more risk averse. lambda=0 is identical to taking the expectation.")
        end
        if beta > 1.0 || beta < 0.0
            error("Beta must be in the range [0, 1]. Increasing values of beta are less risk averse. beta=1 is identical to taking the expectation.")
        end
        new(lambda, beta)
    end
end
NestedAVaR(;lambda=0.0, beta=0.0) = NestedAVaR(lambda, beta)

function modifyprobability!(
    measure::NestedAVaR,                    # risk measure to be overloaded
    newp,      # vector of new probabilities (to by modified in place)
    p::Vector{Float64},      # vector of old probabilities
    m::JuMP.Model,
    x::Vector{Float64},                     # vector of state values
    pi::Vector{Vector{Float64}},            # vector (for each outcome) of dual vectors (dual for each state)
    theta::Vector{Float64}                  # vector of future value/cost values
    )
    ismax = getsense(m) == :Max
    @assert length(newp) == length(p) == length(theta)  # sanity
    q = 0.0                                         # Quantile collected so far
    for i in sortperm(theta, rev=!ismax)                # For each scenario in order
        if q >=  measure.beta                               # We have collected the beta quantile
            riskaverse  = 0.0                       # riskaverse contribution is zero
        else
            riskaverse  = min(p[i], measure.beta-q) / measure.beta  # risk averse contribution
        end
        newp[i] = (1-measure.lambda) * p[i] + measure.lambda * riskaverse  # take the biggest proportion of the scenario possible
        q += riskaverse * measure.beta                      # Update total quantile collected
    end
end
