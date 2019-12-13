#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    ValueFunction

A representation of the value function. SDDP.jl uses the following unique representation of
the value function that is undocumented in the literature.

It supports three types of state variables:

 1) x - convex "resource" states
 2) b - concave "belief" states
 3) y - concave "objective" states

In addition, we have three types of cuts:

 1) Single-cuts (also called "average" cuts in the literature), which involve the
    risk-adjusted expectation of the cost-to-go.
 2) Multi-cuts, which use a different cost-to-go term for each realization w.
 3) Risk-cuts, which correspond to the facets of the dual interpretation of a coherent risk
    measure.

Therefore, ValueFunction returns a JuMP model of the following form:

    V(x, b, y) = min: μᵀb + νᵀy + θ
                 s.t. # "Single" / "Average" cuts
                      μᵀb(j) + νᵀy(j) + θ >= α(j) + xᵀβ(j), ∀ j ∈ J
                      # "Multi" cuts
                      μᵀb(k) + νᵀy(k) + φ(w) >= α(k, w) + xᵀβ(k, w), ∀w ∈ Ω, k ∈ K
                      # "Risk-set" cuts
                      θ ≥ Σ{p(k, w) * φ(w)}_w - μᵀb(k) - νᵀy(k), ∀ k ∈ K
"""
struct ValueFunction{
    O<:Union{Nothing,NTuple{N,JuMP.VariableRef} where {N}},
    B<:Union{Nothing,Dict{T,JuMP.VariableRef} where {T}},
}
    model::JuMP.Model
    theta::JuMP.VariableRef
    states::Dict{Symbol,JuMP.VariableRef}
    objective_state::O
    belief_state::B
end

function _add_to_value_function(
    model::JuMP.Model,
    states::Dict{Symbol,JuMP.VariableRef},
    objective_state,
    belief_state,
    convex_approximation::ConvexApproximation,
    theta_name::String,
)
    theta = @variable(model, base_name = theta_name)
    if objective_sense(model) == MOI.MIN_SENSE
        set_lower_bound(theta, lower_bound(convex_approximation.theta))
    else
        set_upper_bound(theta, upper_bound(convex_approximation.theta))
    end
    for cut in convex_approximation.cut_oracle.cuts
        cut_expr = AffExpr(cut.intercept)
        for (key, coef) in cut.coefficients
            if !haskey(states, key)
                states[key] = @variable(model, base_name = "$(key)")
            end
            add_to_expression!(cut_expr, coef, states[key])
        end
        if objective_state !== nothing
            @assert cut.obj_y !== nothing
            for (y, μ) in zip(cut.obj_y, objective_state)
                add_to_expression!(cut_expr, -y, μ)
            end
        end
        if belief_state !== nothing
            @assert cut.belief_y !== nothing
            for (key, μ) in belief_state
                add_to_expression!(cut_expr, -cut.belief_y[key], μ)
            end
        end
        if objective_sense(model) == MOI.MIN_SENSE
            @constraint(model, theta >= cut_expr)
        else
            @constraint(model, theta <= cut_expr)
        end
    end
    return theta
end

function ValueFunction(node::Node{T}, optimizer) where {T}
    b = node.bellman_function
    sense = objective_sense(node.subproblem)
    model = Model(optimizer)
    set_objective_sense(model, sense)
    states = Dict{Symbol,VariableRef}()
    objective_state = if node.objective_state === nothing
        nothing
    else
        VariableRef[@variable(model, lower_bound = lower_bound(μ), upper_bound = upper_bound(μ)) for μ in node.objective_state.μ]
    end
    belief_state = if node.belief_state === nothing
        nothing
    else
        Dict{T,VariableRef}(
            key => @variable(
                model,
                lower_bound = lower_bound(μ),
                upper_bound = upper_bound(μ),
                base_name = "_belief_$(key)"
            )
            for (key, μ) in node.belief_state.μ
        )
    end
    global_theta = _add_to_value_function(
        model,
        states,
        objective_state,
        belief_state,
        b.global_theta,
        "V",
    )
    local_thetas = VariableRef[_add_to_value_function(
        model,
        states,
        belief_state,
        objective_state,
        l,
        "v$(i)",
    ) for (i, l) in enumerate(b.local_thetas)]
    for risk_set in b.risk_set_cuts
        expr = @expression(model, sum(p * v for (p, v) in zip(risk_set, local_thetas)))
        if sense == MOI.MIN_SENSE
            @constraint(model, global_theta >= expr)
        else
            @constraint(model, global_theta <= expr)
        end
    end
    return ValueFunction(model, global_theta, states, objective_state, belief_state)
end

"""
    evaluate(
        V::ValueFunction,
        point::Dict{Symbol, Float64},
        objective_state = nothing,
        belief_state = nothing
    )

Evaluate the value function `V` at `point` in the state-space.

Returns a tuple containing the height of the function, and the subgradient
w.r.t. the convex state-variables.
"""
function evaluate(
    V::ValueFunction,
    point::Dict{Symbol,Float64};
    objective_state = nothing,
    belief_state = nothing,
)
    for (state, val) in point
        fix(V.states[state], val, force = true)
    end
    saddle = AffExpr(0.0)
    if V.objective_state !== nothing
        @assert objective_state !== nothing
        for (y, x) in zip(objective_state, V.objective_state)
            add_to_expression!(saddle, y, x)
        end
    end
    if V.belief_state !== nothing
        @assert belief_state !== nothing
        for (key, x) in V.belief_state
            add_to_expression!(saddle, belief_state[key], x)
        end
    end
    @objective(V.model, objective_sense(V.model), V.theta + saddle)
    optimize!(V.model)
    obj = objective_value(V.model)
    duals = Dict{Symbol,Float64}()
    for (key, var) in V.states
        duals[key] = dual(FixRef(var))
    end
    return obj, duals
end
