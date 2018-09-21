#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

"""
    modify_probability(measure::AbstractRiskMeasure,
            riskadjusted_distribution,
            original_distribution::Vector{Float64},
            observations::Vector{Float64},
            model::SDDPModel,
            subproblem::JuMP.Model)

Calculate the risk-adjusted probability of each scenario using the
'change-of-probabilities' approach of Philpott, de Matos, and Finardi,(2013). On
solving multistage stochastic programs with coherent risk measures. Operations
Research 61(4), 957-970.

The function should modify `riskadjusted_distribution` in-place based on the
original probabiltiy distribution (contained in `original_distribution`) and the
costs observed in each scenario (contained in `observations`).

The SDDPModel and subproblem are provided for advanced risk measures which may
want to make use of the information contained within.
"""
function modify_probability(measure::AbstractRiskMeasure,
        riskadjusted_distribution,
        original_distribution::Vector{Float64},
        observations::Vector{Float64},
        model::SDDPModel,
        subproblem::JuMP.Model
    )
    error("You need to overload a `modify_probability` method for the measure" *
          " of type $(typeof(measure)).")
end

include("Expectation.jl")
include("WorstCase.jl")
include("AVaR.jl")
include("ConvexCombination.jl")
include("DRO.jl")
include("Wasserstein.jl")

@deprecate modifyprobability! modify_probability
