#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export AVaR

"""
    AVaR(beta::Float64)

The Average Value @ Risk measure. When `beta=0`, the measure is the is
worst-case, when `beta=1` the measure is equivalent to expectation.
"""
struct AVaR <: AbstractRiskMeasure
    beta::Float64
    function AVaR(beta::Float64)
        if !(0.0 <= beta <= 1.0)
            error("Beta must be in the range [0, 1]. Increasing values of " *
                  "beta are less risk averse. beta=1 is identical to taking " *
                  "the expectation.")
        end
        return new(beta)
    end
end

function modify_probability(measure::AVaR,
                            riskadjusted_distribution,
                            original_distribution::Vector{Float64},
                            observations::Vector{Float64},
                            model::SDDPModel,
                            subproblem::JuMP.Model)
    TOLERANCE = 1e-8  # Tolerance for checking extreme points of beta.
    if measure.beta < TOLERANCE
        return modify_probability(WorstCase(), riskadjusted_distribution,
                                  original_distribution, observations, model,
                                  subproblem)
    elseif measure.beta > 1.0 - TOLERANCE
        return modify_probability(Expectation(), riskadjusted_distribution,
                                  original_distribution, observations, model,
                                  subproblem)
    end
    riskadjusted_distribution .= 0.0
    quantile_collected = 0.0
    for i in sortperm(observations, rev=getsense(subproblem) == :Min)
        if quantile_collected >=  measure.beta
            break
        end
        avar_prob = min(original_distribution[i], measure.beta -
                        quantile_collected) / measure.beta
        # take the biggest proportion of the noise possible
        riskadjusted_distribution[i] = avar_prob
        # Update total quantile collected
        quantile_collected += avar_prob * measure.beta
    end
end
