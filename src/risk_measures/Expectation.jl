#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export Expectation

"""
    Expectation()

The expectation risk measure.
"""
struct Expectation <: AbstractRiskMeasure end

function modify_probability(::Expectation,
                            riskadjusted_distribution::AbstractVector,
                            original_distribution::Vector{Float64},
                            observations::Vector{Float64},
                            model::SDDPModel,
                            subproblem::JuMP.Model)
    riskadjusted_distribution .= original_distribution
    return
end
