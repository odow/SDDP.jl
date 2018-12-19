#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ============================== SDDP.AverageCut ===============================

struct AverageCut <: AbstractBellmanFunction
    variable::JuMP.VariableRef
    cut_improvement_tolerance::Float64
end

struct BellmanFactory{T}
    args
    kwargs
    BellmanFactory{T}(args...; kwargs...) where T = new{T}(args, kwargs)
end

"""
    AverageCut(; lower_bound = -Inf, upper_bound = Inf)

The AverageCut Bellman function. Provide a lower_bound if minimizing, or an
upper_bound if maximizing.
"""
function AverageCut(; kwargs...)
    return BellmanFactory{AverageCut}(; kwargs...)
end

function initialize_bellman_function(factory::BellmanFactory{AverageCut},
                                     graph::PolicyGraph{T},
                                     node::Node{T}) where T
    lower_bound, upper_bound = -Inf, Inf
    cut_improvement_tolerance = 0.0
    if length(factory.args) > 0
        error("Positional arguments $(factory.args) ignored in AverageCut.")
    end
    for (kw, value) in factory.kwargs
        if kw == :lower_bound
            lower_bound = value
        elseif kw == :upper_bound
            upper_bound = value
        elseif kw == :cut_improvement_tolerance
            if value < 0
                error("Cut cut_improvement_tolerance must be > 0.")
            end
            cut_improvement_tolerance = value
        else
            error("Keyword $(kw) not recognised as argument to AverageCut.")
        end
    end
    bellman_variable = if length(node.children) > 0
        @variable(node.subproblem,
                  lower_bound = lower_bound, upper_bound = upper_bound)
    else
        @variable(node.subproblem, lower_bound = 0, upper_bound = 0)
    end
    # Initialize bounds for the objective states. If objective_state==nothing,
    # this check will be skipped by dispatch.
    add_initial_bounds(node.objective_state, node.subproblem, bellman_variable)
    return AverageCut(bellman_variable, cut_improvement_tolerance)
end

# Internal function: helper used in add_objective_state_constraint.
function _dot(y::NTuple{N, Float64}, μ::NTuple{N, JuMP.VariableRef}) where {N}
    return sum(y[i] * μ[i] for i in 1:N)
end

# Internal function: helper used in add_initial_bounds.
function add_objective_state_constraint(subproblem, y, μ, θ)
    lower_bound = JuMP.lower_bound(θ)
    upper_bound = JuMP.upper_bound(θ)
    if lower_bound > -Inf
        @constraint(subproblem, _dot(y, μ) + θ >= lower_bound)
    end
    if upper_bound < Inf
        @constraint(subproblem, _dot(y, μ) + θ <= upper_bound)
    end
    if lower_bound ≈ upper_bound ≈ 0.0
        @constraint(subproblem, [i=1:length(μ)], μ[i] == 0.0)
    end
    return
end

# Internal function: When created, θ has bounds of [-M, M], but, since we are
# adding these μ terms, we really want to bound <y, μ> + θ ∈ [-M, M]. We need to
# consider all possible values for `y`. Because the domain of `y` is
# rectangular, we want to add a constraint at each extreme point. This involves
# adding 2^N constraints where N = |μ|. This is only feasible for
# low-dimensional problems, e.g., N < 5.
add_initial_bounds(obj_state::Nothing, subproblem, θ) = nothing
function add_initial_bounds(obj_state::ObjectiveState, subproblem, θ)
    if length(obj_state.μ) < 5
        for y in Base.product(zip(obj_state.lower_bound, obj_state.upper_bound)...)
            add_objective_state_constraint(subproblem, y, obj_state.μ, θ)
        end
    else
        add_objective_state_constraint(
            subproblem, obj_state.lower_bound, obj_state.μ, θ)
        add_objective_state_constraint(
            subproblem, obj_state.upper_bound, obj_state.μ, θ)
    end
end

bellman_term(bellman::AverageCut) = bellman.variable

function refine_bellman_function(graph::PolicyGraph{T},
                                 node::Node{T},
                                 bellman_function::AverageCut,
                                 risk_measure::AbstractRiskMeasure,
                                 outgoing_state::Dict{Symbol, Float64},
                                 dual_variables::Vector{Dict{Symbol, Float64}},
                                 noise_supports::Vector,
                                 original_probability::Vector{Float64},
                                 objective_realizations::Vector{Float64}
                                     ) where T
    is_minimization = graph.objective_sense == MOI.MinSense
    risk_adjusted_probability = similar(original_probability)
    adjust_probability(risk_measure,
                       risk_adjusted_probability,
                       original_probability,
                       noise_supports,
                       objective_realizations,
                       is_minimization)
    # Initialize average cut coefficients.
    intercept = 0.0
    coefficients = Dict{Symbol, Float64}()
    for state in keys(outgoing_state)
        coefficients[state] = 0.0
    end
    # Gather up coefficients for cut calculation.
    # β = F[λ]
    # α = F[θ] - βᵀ ̄x'
    # θ ≥ α + βᵀ x'
    for (objective, dual, prob) in zip(objective_realizations, dual_variables,
                                       risk_adjusted_probability)
        intercept += prob * objective
        for state in keys(outgoing_state)
            coefficients[state] += prob * dual[state]
        end
    end
    # Height of the cut at outgoing_state. We cache the value here for the
    # tolerance check that happens later.
    current_height = intercept
    # Calculate the intercept of the cut.
    for (name, value) in outgoing_state
        intercept -= coefficients[name] * value
    end

    # Coefficients in the objective state dimension.
    objective_state_component = get_objective_state_component(node)

    # Test whether we should add the new cut to the subproblem. We do this now
    # before collating the intercept to avoid twice the work.
    cut_is_an_improvement = if bellman_function.cut_improvement_tolerance > 0.0
        abs(JuMP.objective_value(node.subproblem) - current_height) >
            bellman_function.cut_improvement_tolerance
    else
        true
    end

    if cut_is_an_improvement
        if is_minimization
            @constraint(node.subproblem,
                bellman_function.variable + objective_state_component >=
                    intercept + sum(coefficients[name] * state.out
                        for (name, state) in node.states))
        else
            @constraint(node.subproblem,
                bellman_function.variable + objective_state_component <=
                    intercept + sum(coefficients[name] * state.out
                        for (name, state) in node.states))
        end
    end
    return
end
