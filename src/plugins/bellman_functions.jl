#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

# ============================== SDDP.AverageCut ===============================

mutable struct Cut
    # The intercept of the cut.
    intercept::Float64
    # The coefficients on the state variables.
    coefficients::Dict{Symbol, Float64}
    # A count of the number of points at which the cut is non-dominated.
    non_dominated_count::Int
    # The constraint reference (if it is currently in the model).
    constraint_ref::Union{Nothing, JuMP.ConstraintRef}  # TODO(odow): improve type.
end

mutable struct SampledState
    # A sampled point in state-space.
    state::Dict{Symbol, Float64}
    # Current dominating cut.
    dominating_cut::Cut
    # Current evaluation of the dominating cut at `state`.
    best_objective::Float64
end

struct AverageCut <: AbstractBellmanFunction
    # The cost-to-go variable.
    variable::JuMP.VariableRef
    # Data for Level-One cut selection.
    cuts::Vector{Cut}
    states::Vector{SampledState}
    # Storage for cuts that can be purged.
    cuts_to_be_deleted::Vector{Cut}
    # Minimum number of cuts to purge before activating JuMP.delete.
    deletion_minimum::Int
end

# Internal struct: this struct is just a cache for arguments until we can build
# an actual instance of the type T at a later point.
struct BellmanFactory{T}
    args
    kwargs
    BellmanFactory{T}(args...; kwargs...) where {T} = new{T}(args, kwargs)
end

"""
    AverageCut(; lower_bound = -Inf, upper_bound = Inf, deletion_minimum = 1)

The AverageCut Bellman function. Provide a `lower_bound` if minimizing, or an
`upper_bound` if maximizing.

This Bellman function also implements Level-One cut selection. This requires a
solver that supports constraint deletion. If the solver doesn't support
constraint deletion (e.g., Ipopt), JuMP must be used in Automatic mode instead
of direct mode. In automatic mode, constraint deletion incurs a high cost
because a new model is copied to the solver with every change. Thus, it is
possible to cache constraint deletions until `deletion_minimum` are ready to
improve performance. Cut selection can be "turned off" by setting
`deletion_minimum` to a very large positive integer.
"""
function AverageCut(; lower_bound = -Inf, upper_bound = Inf,
                    deletion_minimum::Int = 1)
    return BellmanFactory{AverageCut}(lower_bound = lower_bound,
        upper_bound = upper_bound, deletion_minimum = deletion_minimum)
end

function initialize_bellman_function(factory::BellmanFactory{AverageCut},
                                     graph::PolicyGraph{T},
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
    bellman_variable = @variable(node.subproblem, lower_bound = lower_bound,
                                 upper_bound = upper_bound)
    # Initialize bounds for the objective states. If objective_state==nothing,
    # this check will be skipped by dispatch.
    add_initial_bounds(node.objective_state, node.subproblem, bellman_variable)
    return AverageCut(
        bellman_variable, Cut[], SampledState[], Cut[], deletion_minimum
    )
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
                                     ) where {T}
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
    # Initialize the cut struct. It gets initialized with a non_dominated_count
    # of 1 because if there is cut selection, it dominates at the most recent
    # point (`outgoing_state`).
    cut = Cut(intercept, coefficients, 1, nothing)
    # No cut selection if there is objective-state interpolation :(.
    if objective_state_component == JuMP.AffExpr(0.0)
        levelone_update(bellman_function, cut, outgoing_state, is_minimization)
        purge_cuts(node, bellman_function)
    end
    add_new_cut(node, bellman_function.variable, cut, objective_state_component,
                is_minimization)
    return
end

# Internal function: delete dominated cuts.
function purge_cuts(node::Node, bellman::AverageCut)
    if length(bellman.cuts_to_be_deleted) >= bellman.deletion_minimum
        for existing_cut in bellman.cuts_to_be_deleted
            JuMP.delete(node.subproblem, existing_cut.constraint_ref)
        end
        empty!(bellman.cuts_to_be_deleted)
    end
    return
end

# Internal function: add a new cut to the model.
function add_new_cut(node::Node, theta::JuMP.VariableRef, cut::Cut,
                     objective_state, is_minimization::Bool)
    if is_minimization
        cut.constraint_ref = @constraint(node.subproblem,
            theta + objective_state >= cut.intercept +
                sum(cut.coefficients[name] * state.out
                    for (name, state) in node.states))
    else
        cut.constraint_ref = @constraint(node.subproblem,
            theta + objective_state <= cut.intercept +
                sum(cut.coefficients[name] * state.out
                    for (name, state) in node.states))
    end
    return
end

# Internal function: calculate the height of `cut` evaluated at `state`.
function eval_height(cut::Cut, state::Dict{Symbol, Float64})
    height = cut.intercept
    for (key, value) in cut.coefficients
        height += value * state[key]
    end
    return height
end

# Internal function: check if the candidate point dominates the incumbent.
function dominates(candidate, incumbent, minimization::Bool)
    return minimization ? candidate > incumbent : candidate < incumbent
end

# Internal function: update the Level-One datastructures inside
# `bellman_function`.
function levelone_update(bellman_function::AverageCut, cut::Cut,
                         state::Dict{Symbol, Float64}, is_minimization)
    sampled_state = SampledState(state, cut, eval_height(cut, state))
    # Loop through previously sampled states and compare the height of the most
    # recent cut against the current best. If this cut is an improvement, store
    # this one instead.
    for old_state in bellman_function.states
        height = eval_height(cut, old_state.state)
        if dominates(height, old_state.best_objective, is_minimization)
            old_state.dominating_cut.non_dominated_count -= 1
            if old_state.dominating_cut.non_dominated_count <= 0
                push!(bellman_function.cuts_to_be_deleted,
                    old_state.dominating_cut)
            end
            cut.non_dominated_count += 1
            old_state.dominating_cut = cut
            old_state.best_objective = height
        end
    end
    # Now loop through previously discovered cuts and compare their height at
    # the most recent sampled point in the state-space. If this cut is an
    # improvement, store this one instead. Note that we have to do this because
    # we might have previously thrown out a cut that is now relevant.
    for old_cut in bellman_function.cuts
        # If the constriant ref is not nothing, this cut is already in the
        # model, so it can't be better than the one we just found.
        if old_cut.constraint_ref !== nothing
            continue
        elseif !JuMP.is_valid(old_cut.constraint_ref)
            old_cut.constraint_ref = nothing
            continue
        end
        height = eval_height(old_cut, state)
        if dominates(height, sampled_state.best_objective, is_minimization)
            sampled_state.dominating_cut.non_dominated_count -= 1
            if sampled_state.dominating_cut.non_dominated_count <= 0
                push!(bellman_function.cuts_to_be_deleted,
                    sampled_state.dominating_cut)
            end
            old_cut.non_dominated_count += 1
            sampled_state.dominating_cut = old_cut
            sampled_state.best_objective = height
        end
    end
    push!(bellman_function.cuts, cut)
    push!(bellman_function.states, sampled_state)
    return
end
