#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

"""
    modifyprobability!(measure::AbstractRiskMeasure, newprobabilities, oldprobabilities, x)

    This function assembles a new cut using the following inputs
    + measure::AbstractRiskMeasure - used to dispatch
    + newprobabilities::Vector{Float64} - the probability support of the scenarios. Should sum to one
    + oldprobabilities::Vector{Float64} - the probability support of the scenarios. Should sum to one
    + x::Vector{Float64}                - objectives
"""
modifyprobability!(measure::AbstractRiskMeasure, newprobabilities::Vector{Float64}, oldprobabilities::Vector{Float64}, x::Vector{Float64}) = error("You need to overload a `modifyprobability` method for the measure of type $(typeof(measure)).")

# a more expansive method that can be overloaded
modifyprobability!(
    measure::AbstractRiskMeasure,           # risk measure to be overloaded
    newprobabilities::Vector{Float64},      # vector of new probabilities (to by modified in place)
    oldprobabilities::Vector{Float64},      # vector of old probabilities
    m::JuMP.Model,
    x::Vector{Float64},                     # vector of state values
    pi::Vector{Vector{Float64}},            # vector (for each outcome) of dual vectors (dual for each state)
    theta::Vector{Float64}                  # vector of future value/cost values
    ) = modifyprobability!(measure, newprobabilities, oldprobabilities, theta)


# ==============================================================================
#   The Expectation risk measure:
#   In expectation, leave probabilities as they were

struct Expectation <: AbstractRiskMeasure end

modifyprobability!(
    measure::Expectation,
    newprobabilities::Vector{Float64},
    oldprobabilities::Vector{Float64},
    x::Vector{Float64}
    ) = (newprobabilities .= oldprobabilities)

# ==============================================================================
#   Average Value at Risk
#   (1 - λ) * E[x] + λ * AV@R(1-β)[x]
struct NestedAVaR <: AbstractRiskMeasure
    lambda::Float64
    beta::Float64
    storage::Vector{Float64}
end
NestedAVaR(lamba, beta) = NestedAVaR(lambda, beta, Float64[])
NestedAVaR(;lamba=0.0, beta=0.0) = NestedAVaR(lambda, beta, Float64[])
