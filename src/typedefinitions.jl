#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

abstract type AbstractRiskMeasure end
abstract type AbstractCutOracle end
abstract type OptimisationSense end
struct Max <: OptimisationSense end
struct Min <: OptimisationSense end

const LinearConstraint=JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

struct Cut
    intercept::Float64
    coefficients::Vector{Float64}
end

struct State
    variable::JuMP.Variable
    constraint::LinearConstraint
end

struct Scenario
    # list of row indices
    constraints::Vector{LinearConstraint}
    # list of RHS values
    values::Vector{Float64}
end

abstract type AbstractPriceOracle end

# # mutable struct RibPriceOracle{T} <: AbstractPriceOracle
# #     pricetransition::Function  # ℜ² → ℜ
# #     pricescenarios::Vector{T}
# #     objective::Function        # ℜ → AffExpr
# #     ribs::Vector{Float64}
# #     thetas::Vector{JuMP.Variable}
# #     cutoracles::Vector{CutOracles}
# # end
# # PriceOracle() = PriceOracle((p)->p, Float64[], (p) -> AffExpr(p))

struct DefaultPriceOracle{T<:AbstractCutOracle} <: AbstractPriceOracle
    theta::JuMP.Variable
    cutoracle::T
end

struct SubproblemExt{R<:AbstractRiskMeasure, S<:OptimisationSense}
    stage::Int               # stage index
    markovstate::Int         # index of the subproblem by markov state

    states::Vector{State}                # a vector of states

    # priceoracle::AbstractPriceOracle
    # cuts::Vector{Cut}

    # vector of scenarios
    scenarios::Vector{Scenario}
    # probability[i] = probability of scenarios[i] occuring
    scenarioprobability::Vector{Float64}

    riskmeasure::R # A risk measure to use for the subproblem

    problembound::Float64
    sense::Type{S}
end
ext(m::JuMP.Model) = m.ext[:SDDP]::SubproblemExt

function Subproblem()
    m = Model()
    m.ext[:SDDP] = SDDP.SubproblemExt(
        1,
        1,
        State[],
        # ValueFunction[],
        # AbstractCutOracle[],
        # DefaultPriceOracle(),
        Scenario[],
        Float64[],
        Expectation(),
        0.0,
        Min
    )
    m
end

struct Stage
    # vector of subproblems in this stage
    subproblems::Vector{JuMP.Model}
    # transitionprobabilities[i] =
    # probability of transitioning to subproblem i in next stage
    transitionprobabilities::Vector{Float64}
end

struct SDDPModel
    stages::Vector{Stage}
end
