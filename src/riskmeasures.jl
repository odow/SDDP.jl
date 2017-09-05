#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    modifyprobability!(measure::AbstractRiskMeasure,
            riskadjusted_distribution,
            original_distribution::Vector{Float64},
            observations::Vector{Float64},
            m::SDDPModel,
            sp::JuMP.Model
    )

# Description

Perform a 'change-of-probabilities' transformation.

# Arguments
 * `measure::AbstractRiskMeasure`
 The risk measure
 * `riskadjusted_distribution`
 * `original_distribution::Vector{Float64}`
 The original probability distribution.
 * `observations::Vector{Float64}`
 The vector of objective values from the next stage
 problems (one for each scenario).
 * `m::SDDPModel`
 The full SDDP model
 * `sp::JuMP.Model`
 The stage problem that the cut will be added to.
"""
function modifyprobability!(measure::AbstractRiskMeasure,
        riskadjusted_distribution,
        original_distribution::Vector{Float64},
        observations::Vector{Float64},
        m::SDDPModel,
        sp::JuMP.Model
    )
    error("You need to overload a `modifyprobability!` method for the measure of type $(typeof(measure)).")
end

# ==============================================================================
#   The Expectation risk measure:
#   In expectation, leave probabilities as they were

immutable Expectation <: AbstractRiskMeasure end

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
