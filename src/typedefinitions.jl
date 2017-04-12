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

const LinearConstraint=JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

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

struct PriceOracle{T}
    pricetransition::Function  # ℜ → ℜ
    pricescenarios::Vector{T}
    objective::Function        # ℜ → AffExpr
end
PriceOracle() = PriceOracle((p)->p, Float64[], (p) -> AffExpr(p))

struct ValueFunction
    location::Float64
    theta::JuMP.Variable
    cuts::AbstractCutOracle # a cut oracle
end

struct SubproblemExt
    stage::Int               # stage index
    markovstate::Int         # index of the subproblem by markov state

    states::Vector{State}                # a vector of states

    valuefunctions::Vector{ValueFunction} # a value function
    priceoracle::PriceOracle

    # vector of scenarios
    scenarios::Vector{Scenario}
    # probability[i] = probability of scenarios[i] occuring
    scenarioprobability::Vector{Float64}

    riskmeasure::AbstractRiskMeasure # A risk measure to use for the subproblem
end
ext(m::JuMP.Model) = m.ext[:SDDP]::SubproblemExt

function Subproblem()
    m = Model()
    m.ext[:SDDP] = SDDP.SubproblemExt(
        1,
        1,
        SDDP.State[],
        SDDP.ValueFunction[],
        SDDP.PriceOracle(),
        SDDP.Scenario[],
        Float64[],
        Expectation()
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
