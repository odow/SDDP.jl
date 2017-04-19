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
abstract type AbstractValueFunction end

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

struct SubproblemExt{S<:OptimisationSense, V<:AbstractValueFunction, R<:AbstractRiskMeasure}
    stage::Int              # stage index
    markovstate::Int        # index of the subproblem by markov state
    problembound::Float64   # objective bound
    sense::Type{S}          # optimisation sense (max or min)
    # a vector of states
    states::Vector{State}
    # an oracle to value function
    valueoracle::V
    # vector of scenarios
    scenarios::Vector{Scenario}
    # probability[i] = probability of scenarios[i] occuring
    scenarioprobability::Vector{Float64}
    # A risk measure to use for the subproblem
    riskmeasure::R
end
ext(m::JuMP.Model) = m.ext[:SDDP]::SubproblemExt
isext(m::JuMP.Model) = isa(m.ext[:SDDP], SubproblemExt)

function Subproblem(;stage=1, markov_state=1, sense=Min, bound=-1e6,
    risk_measure=Expectation(), cut_oracle=DefaultCutOracle(), value_function=DefaultValueFunction)

    m = Model()
    m.ext[:SDDP] = SDDP.SubproblemExt(
        stage,
        markov_state,
        bound,
        sense,
        State[],
        init!(value_function, m, sense, bound, cut_oracle),
        Scenario[],
        Float64[],
        risk_measure
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
Stage(transition=Float64[]) = Stage(JuMP.Model[], transition)

mutable struct BackwardPassStorage
    state::Vector{Float64}
    duals::Vector{Vector{Float64}}
    objective::Vector{Float64}
    probability::Vector{Float64}
    newprobability::Vector{Float64}
    n::Int
end
BackwardPassStorage() = BackwardPassStorage(Float64[], Vector{Float64}[], Float64[], Float64[], Float64[], 0)

struct ForwardPassStorage
    states::Vector{Vector{Float64}}
end

struct SDDPModel
    stages::Vector{Stage}
    storage::BackwardPassStorage
    forwardstorage::ForwardPassStorage
end
SDDPModel() = SDDPModel(Stage[], BackwardPassStorage())

struct Timer
    simulation::Float64
    cutting::Float64
    lpsolver::Float64
    riskmeasure::Float64
    cutselection::Float64
end
