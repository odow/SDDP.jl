#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export WorstCase

"""
    WorstCase()

The worst-case risk measure.
"""
struct WorstCase <: AbstractRiskMeasure end

function modify_probability(measure::WorstCase,
                            riskadjusted_distribution,
                            original_distribution::Vector{Float64},
                            observations::Vector{Float64},
                            model::SDDPModel,
                            subproblem::JuMP.Model)
    sense = getsense(subproblem)
    riskadjusted_distribution .= 0.0
    worst_idx = 1
    worst_observation = -worstcase(sense)
    for (idx, (probability, observation)) in enumerate(zip(original_distribution, observations))
        if probability > 0.0 && dominates(sense, observation, worst_observation)
            worst_idx = idx
            worst_observation = observation
        end
    end
    riskadjusted_distribution[worst_idx] = 1.0
    return
end
