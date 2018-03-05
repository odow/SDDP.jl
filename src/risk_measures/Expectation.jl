#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    Expectation()

# Description

The expectation risk measure.
"""
struct Expectation <: AbstractRiskMeasure end

function modifyprobability!(::Expectation,
    riskadjusted_distribution::AbstractVector,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    riskadjusted_distribution .= original_distribution
    return nothing
end
