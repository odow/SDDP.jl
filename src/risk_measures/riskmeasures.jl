#  Copyright 2018, Oscar Dowson
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

Calculate the risk-adjusted probability of each scenario using the
'change-of-probabilities' approach of Philpott, de Matos, and Finardi,(2013). On
solving multistage stochastic programs with coherent risk measures. Operations
Research 61(4), 957-970.

# Arguments
 * `measure::AbstractRiskMeasure`
 The risk measure
 * `riskadjusted_distribution`
 A new probability distribution
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

include("Expectation.jl")
include("WorstCase.jl")
include("AVaR.jl")
include("ConvexCombination.jl")
include("DRO.jl")
include("Wasserstein.jl")
