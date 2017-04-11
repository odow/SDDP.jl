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

const LinearConstraint=JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

struct State
    variable::JuMP.Variable
    constraint::LinearConstraint
end

struct Scenario
    # list of row indices
    row::Vector{LinearConstraint}
    # list of RHS values
    value::Vector{Float64}
end
Scenario() = Scenario(LinearConstraint[], Float64[])

struct Scenarios
    # vector of scenarios
    scenarios::Vector{Scenario}
    # probability[i] = probability of scenarios[i] occuring
    probability::Vector{Float64}
end
Scenarios() = Scenarios(Scenario[], Float64[])

import Base: length
Base.length(s::Scenarios) = length(s.scenarios)


struct ValueFunction
    theta::JuMP.Variable
    cuts#::AbstractCutOracle # a cut oracle
end

struct SubproblemExt
    stage::Int               # stage index
    markovstate::Int         # index of the subproblem by markov state

    states::Vector{State}                # a vector of states
    valuefunction::Vector{ValueFunction} # a value function

    scenarios::Scenarios

    riskmeasure::AbstractRiskMeasure # A risk measure to use for the subproblem
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
