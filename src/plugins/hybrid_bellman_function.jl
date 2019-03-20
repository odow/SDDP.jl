#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

abstract type AbstractCutOracle end

const STATE = Dict{Symbol, Float64}

mutable struct Cut
    intercept::Float64
    coefficients::STATE
    constraint_ref::Union{Nothing, JuMP.ConstraintRef}
end

mutable struct SampledState
    state::STATE
    dominating_cut_index::Int
    best_objective::Float64
end

struct LevelOneOracle <: AbstractCutOracle
    cuts::Vector{Cut}
    states::Vector{SampledState}
    cuts_to_be_deleted::Vector{Cut}
    deletion_minimum::Int
end

struct ConvexApproximation
    θ::JuMP.VariableRef
    cut_oracle::LevelOneOracle
end

struct BellmanFunction <: AbstractBellmanFunction
    θ_global::JuMP.VariableRef
    θ_locals::Vector{JuMP.VariableRef}
end

function BellmanFunction(; lower_bound = -Inf, upper_bound = Inf,
                         deletion_minimum::Int = 1)
    return BellmanFactory{BellmanFunction}(lower_bound = lower_bound,
        upper_bound = upper_bound, deletion_minimum = deletion_minimum)
end

function initialize_bellman_function(
        factory::BellmanFactory{BellmanFunction}, model::PolicyGraph{T},
        node::Node{T}) where {T}
    lower_bound, upper_bound = -Inf, Inf
    deletion_minimum = 0
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in AverageCut.")
    end
    for (kw, value) in factory.kwargs
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :deletion_minimum
            deletion_minimum = value
        else
            error("Keyword $(kw) not recognised as argument to AverageCut.")
        end
    end
    if length(node.children) == 0
        lower_bound = upper_bound = 0.0
    end
    bellman_variable = if lower_bound == -Inf && upper_bound < Inf
        @variable(node.subproblem, upper_bound = upper_bound)
    elseif lower_bound > -Inf && upper_bound == Inf
        @variable(node.subproblem, lower_bound = lower_bound)
    elseif lower_bound > -Inf && upper_bound < Inf
        @variable(node.subproblem,
            lower_bound = lower_bound, upper_bound = upper_bound)
    else
        error("You must supply a non-infinite bound for the cost-to-go variable.")
    end
    return BellmanFunction(θ, JuMP.VariableRef[])
end

@enum(CutType, AVERAGE_CUT, MULTI_CUT)

function refine_approximation(
            model::PolicyGraph{T},
            node::Node{T},
            bellman_function::BellmanFunction,
            risk_measure::AbstractRiskMeasure,
            outgoing_state::STATE,
            dual_variables::Vector{STATE},
            noise_supports::Vector,
            nominal_probability::Vector{Float64},
            objective_realizations::Vector{Float64},
            cut_type::CutType
        ) where {T}
    # Sanity checks.
    @assert length(dual_variables) == length(noise_supports) ==
        length(nominal_probability) == length(objective_realizations)
    # Preliminaries that are common to all cut types.
    risk_adjusted_probability = similar(nominal_probability)
    adjust_probability(
        risk_measure, risk_adjusted_probability, nominal_probability,
        noise_supports, objective_realizations,
        model.objective_sense == MOI.MIN_SENSE)
    # The meat of the function.
    if cut_type == AVERAGE_CUT
        _add_average_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            duals)
    else  # Add a multi-cut
        @assert cut_type == MULTI_CUT
        _add_locals_if_necessary(bellman_function, length(dual_variables))
        _add_multi_cut(
            node,
            outgoing_state,
            risk_adjusted_probability,
            objective_realizations,
            dual_variables)
    end
end

# Add the average cut `Θ ≥ Eᵢ[θᵢ] + ⟨Eᵢ[πᵢ], x′ - xᵏ⟩`
function _add_average_cut(node::Node, outgoing_state::STATE,
                          risk_adjusted_probability::Vector{Float64},
                          objective_realizations::Vector{Float64},
                          duals::Vector{STATE})
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(duals)
    coefficients = STATE(key => 0.0 for key in keys(outgoing_state))
    intercept = 0.0
    for i in 1:length(objective_realizations)
        p = risk_adjusted_probability[i]
        intercept += p * objective_realizations[i]
        for (key, dual) in duals[i]
            coefficients[key] += p * dual
        end
    end
    for (key, dual) in coefficients
        intercept -= dual * outgoing_state[key]
    end
    θ = node.bellman_function.θ_global
    model = JuMP.owner_model(θ)
    terms = @expression(model,
        intercept + sum(coefficients[name] * x′ for (name, x′) in node.states))
    _add_bounding_constraint(model, θ, terms)
    return
end

function _add_bounding_constraint(model, θ, terms)
    if JuMP.objective_sense(model) == MOI.MIN_SENSE
        @constraint(model, θ >= terms)
    else
        @constraint(model, θ <= terms)
    end
end

# Add the cut `θ ≥ θᵏ + ⟨πᵏ, x′ - xᵏ⟩`.
function _add_local_cut(θ, x′, θᵏ, xᵏ, πᵏ)
    model = JuMP.owner_model(θ)
    terms = @expression(model, θᵏ + sum(πᵏ[i] * (x′[i] - xᵏ[i]) for i in keys(x′))
    _add_bounding_constraint(model, θ, terms)
    return
end

function _add_multi_cut(node::Node, outgoing_state::STATE,
                        risk_adjusted_probability::Vector{Float64},
                        objective_realizations::Vector{Float64},
                        duals::Vector{STATE})
    N = length(risk_adjusted_probability)
    @assert N == length(objective_realizations) == length(duals)
    bellman_function = node.bellman_function
    for i in 1:length(duals)
        _add_local_cut(
            bellman_function.θ_locals[i],
            node.states,
            objective_realizations[i],
            outgoing_state,
            dual_variables[i])
    end
    model = JuMP.owner_model(bellman.θ_global)
    terms = @expression(model, sum(
        risk_adjusted_probability[i] * bellman_function.θ_locals[i]
            for i in 1:length(risk_adjusted_probability)))
    _add_bounding_constraint(model, bellman_function.ϴ_global, terms)
    return
end

# If we are adding a multi-cut for the first time, then the local θ variables
# won't have been added.
# TODO(odow): a way to set different bounds for each variable in the multi-cut.
function _add_locals_if_necessary(bellman::BellmanFunction, N::Int)
    num_local_thetas = length(bellman.θ_locals)
    if num_local_thetas == N
        # Do nothing. Already initialized.
    elseif num_local_thetas == 0
        model = JuMP.owner_model(bellman.θ_global)
        θ = @variable(model, [1:N])
        if JuMP.has_lower_bound(bellman.θ_global)
            JuMP.set_lower_bound.(θ, JuMP.lower_bound(bellman.θ_global))
        end
        if JuMP.has_upper_bound(bellman.θ_global)
            JuMP.set_upper_bound.(θ, JuMP.upper_bound(bellman.θ_global))
        end
        resize!(bellman.θ_locals, N)
        copyto!(bellman.θ_locals, θ)
    else
        error("Expected $(N) local θ variables but there were $(num_local_thetas).")
    end
    return
end
