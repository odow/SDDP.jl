#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

struct AVaR <: AbstractRiskMeasure
    beta::Float64
    function AVaR(beta::Float64)
        if 0.0 <= beta <= 1.0
            return new(beta)
        else
            error("Beta must be in the range [0, 1]. Increasing values of beta are less risk averse. beta=1 is identical to taking the expectation.")
        end
    end
end

function modifyprobability!(measure::AVaR,
    riskadjusted_distribution,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    if measure.beta < 1e-8
        return modifyprobability!(WorstCase(), riskadjusted_distribution, original_distribution, observations, m, sp)
    elseif measure.beta > 1.0 - 1e-8
        return modifyprobability!(Expectation(), riskadjusted_distribution, original_distribution, observations, m, sp)
    else
        ismax = getsense(sp) == :Max
        riskadjusted_distribution .= 0.0
        q = 0.0 # Quantile collected so far
        for i in sortperm(observations, rev=!ismax) # For each noise in order
            if q >=  measure.beta # We have collected the beta quantile
                break
            end
            avar_prob  = min(original_distribution[i], measure.beta-q) / measure.beta
            # take the biggest proportion of the noise possible
            riskadjusted_distribution[i] = avar_prob
            # Update total quantile collected
            q += avar_prob * measure.beta
        end
    end
end
