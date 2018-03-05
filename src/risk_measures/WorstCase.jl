#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    WorstCase()

The worstcase risk measure.
"""
struct WorstCase <: AbstractRiskMeasure end

function modifyprobability!(measure::WorstCase,
    riskadjusted_distribution,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    sense = getsense(sp)
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
#
# function modifyprobability!(measure::WorstCase,
#     riskadjusted_distribution,
#     original_distribution::Vector{Float64},
#     observations::Vector{Float64},
#     m::SDDPModel,
#     sp::JuMP.Model
#     )
#     riskadjusted_distribution .= 0.0
#     probabilities = original_distribution .> 0.0
#     if getsense(sp) == :Min
#         ind = indmax(observations[probabilities])
#     else
#         ind = indmin(observations[probabilities])
#     end
#     riskadjusted_distribution[probabilities][ind] = 1.0
#     return
# end
