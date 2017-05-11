#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

abstract type AbstractRiskMeasure end
abstract type AbstractCutOracle end
abstract type AbstractValueFunction end

abstract type OptimisationSense end
struct Max <: OptimisationSense end
struct Min <: OptimisationSense end

abstract type IterationDirection end
struct ForwardPass <: IterationDirection end
struct BackwardPass <: IterationDirection end


abstract type SDDPSolveType end
struct Serial <: SDDPSolveType end
Base.show(io::IO, async::Serial) = print(io, "Serial solver")


const LinearConstraint=JuMP.ConstraintRef{JuMP.Model, JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64, JuMP.Variable}}}

mutable struct CachedVector{T} <: AbstractArray{T, 1}
    data::Vector{T}
    n::Int
end

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
    finalstage::Bool        # if final stage
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

function Subproblem(;finalstage=false, stage=1, markov_state=1, sense=Min, bound=-1e6,
    risk_measure=Expectation(), cut_oracle=DefaultCutOracle(), value_function=DefaultValueFunction{typeof(cut_oracle)})

    m = Model()
    m.ext[:SDDP] = SDDP.SubproblemExt(
        finalstage,
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
    t::Int
    # vector of subproblems in this stage
    subproblems::Vector{JuMP.Model}
    # transitionprobabilities[i, j] =
    # probability of transitioning from subproblem i to subproblem j in next stage
    transitionprobabilities::Array{Float64, 2}
    # storage for state on forward pass
    state::Vector{Float64}
    # extension dictionary
    ext::Dict
end
Stage(t=1, transition=Array{Float64}(0,0)) = Stage(t, JuMP.Model[], transition, Float64[], Dict())

struct Storage
    state::Vector{Float64}
    scenario::CachedVector{Int}
    markov::CachedVector{Int}
    duals::CachedVector{Vector{Float64}}
    objective::CachedVector{Float64}
    probability::CachedVector{Float64}
    modifiedprobability::CachedVector{Float64}
end
Storage() = Storage(
    Float64[],
    CachedVector(Int),
    CachedVector(Int),
    CachedVector(Vector{Float64}),
    CachedVector(Float64),
    CachedVector(Float64),
    CachedVector(Float64)
)
function reset!(s::Storage)
    reset!(s.scenario)
    reset!(s.markov)
    reset!(s.duals)
    reset!(s.objective)
    reset!(s.probability)
    reset!(s.modifiedprobability)
end

struct SolutionLog
    iteration::Int
    bound::Float64
    lower_statistical_bound::Float64
    upper_statistical_bound::Float64
    timecuts::Float64
    simulations::Int
    timesimulations::Float64
    timetotal::Float64
end
SolutionLog() = SolutionLog(0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0)

struct SDDPModel{V<:AbstractValueFunction}
    sense::Symbol
    stages::Vector{Stage}
    storage::Storage
    log::Vector{SolutionLog}
    build!::Function
    lpsolver::JuMP.MathProgBase.AbstractMathProgSolver
    ext::Dict # extension dictionary
end
newSDDPModel(sense::Symbol, v::AbstractValueFunction, build!::Function, solver::JuMP.MathProgBase.AbstractMathProgSolver) = newSDDPModel(sense, typeof(v), build!, solver)
newSDDPModel{V<:AbstractValueFunction}(sense::Symbol, v::Type{V}, build!::Function, solver::JuMP.MathProgBase.AbstractMathProgSolver) = SDDPModel{V}(sense, Stage[], Storage(), SolutionLog[], build!, solver, Dict())

struct BoundConvergence
    iterations::Int
    rtol::Float64
    atol::Float64
end
BoundConvergence(;iterations=0,rtol=0.0,atol=0.0) = BoundConvergence(iterations,rtol,atol)
struct MonteCarloSimulation
    frequency::Int
    steps::Vector{Int}
    confidence::Float64
    termination::Bool
end
MonteCarloSimulation(;frequency=0,min=20,max=0,step=1,confidence=0.95,termination=false) = MonteCarloSimulation(frequency,collect(min:step:max),confidence,termination)

struct Settings
    max_iterations::Int
    time_limit::Float64
    simulation::MonteCarloSimulation
    bound_convergence::BoundConvergence
    cut_selection_frequency::Int
    print_level::Int
    log_file::String
    reduce_memory_footprint::Bool
    cut_output_file::String
end
Settings() = Settings(0,600.0, MonteCarloSimulation(), BoundConvergence(), 0,0,"", false, "")
