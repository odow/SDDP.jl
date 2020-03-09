#  Copyright 2017-20, Oscar Dowson.
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

Base.show(io::IO, v::ValueFunction) = "An SDDP value function"

JuMP.set_optimizer(v::ValueFunction, optimizer) = set_optimizer(v.model, optimizer)

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
        cut_expr = @expression(
            model,
            cut.intercept + sum(coef * states[key] for (key, coef) in cut.coefficients)
        )
        if objective_state !== nothing
            @assert cut.obj_y !== nothing
            cut_expr = @expression(
                model,
                cut_expr - sum(y * μ for (y, μ) in zip(cut.obj_y, objective_state))
            )
        end
        if belief_state !== nothing
            @assert cut.belief_y !== nothing
            cut_expr = @expression(
                model,
                cut_expr - sum(cut.belief_y[key] * μ for (key, μ) in belief_state)
            )
        end
        if objective_sense(model) == MOI.MIN_SENSE
            @constraint(model, theta >= cut_expr)
        else
            @constraint(model, theta <= cut_expr)
        end
    end
    return theta
end

function ValueFunction(node::Node{T}) where {T}
    b = node.bellman_function
    sense = objective_sense(node.subproblem)
    model = Model()
    if node.optimizer !== nothing
        set_optimizer(model, node.optimizer)
    end
    set_objective_sense(model, sense)
    states = Dict{Symbol,VariableRef}(
        key => @variable(model, base_name = "$(key)") for (key, x) in node.states
    )
    objective_state = if node.objective_state === nothing
        nothing
    else
        tuple(VariableRef[
            @variable(
                model,
                lower_bound = lower_bound(μ),
                upper_bound = upper_bound(μ),
                base_name = "_objective_state_$(i)"
            ) for (i, μ) in enumerate(node.objective_state.μ)
        ]...,)
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
            ) for (key, μ) in node.belief_state.μ
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
    local_thetas = VariableRef[
        _add_to_value_function(model, states, belief_state, objective_state, l, "v$(i)") for (i, l) in enumerate(b.local_thetas)
    ]
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
    sign = objective_sense(V.model) == MOI.MIN_SENSE ? 1.0 : -1.0
    for (key, var) in V.states
        duals[key] = sign * dual(FixRef(var))
    end
    return obj, duals
end

"""
    evalute(V::ValueFunction{Nothing, Nothing}; kwargs...)

Evalute the value function `V` at the point in the state-space specified by `kwargs`.

### Example

    evaluate(V; volume = 1)
"""
function evaluate(V::ValueFunction{Nothing,Nothing}; kwargs...)
    return evaluate(V, Dict(k => float(v) for (k, v) in kwargs))
end

struct Point{Y,B}
    x::Dict{Symbol,Float64}
    y::Y
    b::B
end
Point(x::Dict{Symbol,Float64}) = Point(x, nothing, nothing)

function height(V::ValueFunction{Y,B}, x::Point{Y,B}) where {Y,B}
    return evaluate(V, x.x; objective_state = x.y, belief_state = x.b)[1]
end

function get_axis(x::Vector{Dict{K,V}}) where {K,V}
    @assert length(x) >= 2
    changing_key = nothing
    for (key, val) in x[1]
        if val == x[2][key]
            continue
        elseif changing_key !== nothing
            error("Too many elements are changing")
        end
        changing_key = key
    end
    return changing_key === nothing ? nothing : [xi[changing_key] for xi in x]
end

function get_axis(x::Vector{NTuple{N,T}}) where {N,T}
    @assert length(x) >= 2
    changing_index = nothing
    for i = 1:N
        if x[1][i] == x[2][i]
            continue
        elseif changing_index !== nothing
            error("Too many elements are changing")
        end
        changing_index = i
    end
    return changing_index === nothing ? nothing : [xi[changing_index] for xi in x]
end

get_axis(x::Vector{Nothing}) = nothing

function get_axis(X::Vector{Point{Y,B}}) where {Y,B}
    for f in [x -> x.x, x -> x.y, x -> x.b]
        x = get_axis(f.(X))
        x !== nothing && return x
    end
    return nothing
end

function get_data(V::ValueFunction{Y,B}, X::Vector{Point{Y,B}}) where {Y,B}
    x = get_axis(X)
    if x === nothing
        error("Unable to detect changing dimension")
    end
    y = height.(Ref(V), X)
    return x, y, Float64[]
end

function get_data(V::ValueFunction{Y,B}, X::Matrix{Point{Y,B}}) where {Y,B}
    x = get_axis(collect(X[:, 1]))
    if x === nothing
        error("Unable to detect changing row")
    end
    y = get_axis(collect(X[1, :]))
    if y === nothing
        error("Unable to detect changing column")
    end
    z = height.(Ref(V), X)
    return [i for _ in y for i in x], [i for i in y for _ in x], vec(z)
end

function plot(
    V::ValueFunction{Y,B},
    X::Array{Point{Y,B}};
    filename::String = joinpath(tempdir(), string(Random.randstring(), ".html")),
    open::Bool = true,
) where {Y,B}
    x, y, z = get_data(V, X)
    fill_template(
        filename,
        "<!--X-->" => JSON.json(x),
        "<!--Y-->" => JSON.json(y),
        "<!--Z-->" => JSON.json(z);
        template = joinpath(@__DIR__, "value_functions.html"),
        launch = open,
    )
    return
end


function plot(
    V::ValueFunction{Nothing,Nothing};
    filename::String = joinpath(tempdir(), string(Random.randstring(), ".html")),
    open::Bool = true,
    kwargs...,
)
    d = Dict{Symbol,Float64}()
    variables = Symbol[]
    for (key, val) in kwargs
        if isa(val, AbstractVector)
            push!(variables, key)
        else
            d[key] = float(val)
        end
    end
    if length(variables) == 1
        points = Point{Nothing,Nothing}[]
        key = variables[1]
        for val in kwargs[key]
            d2 = copy(d)
            d2[key] = val
            push!(points, Point(d2))
        end
        return plot(V, points; filename = filename, open = open)
    elseif length(variables) == 2
        k1, k2 = variables
        N1, N2 = length(kwargs[k1]), length(kwargs[k2])
        points = Array{Point{Nothing,Nothing},2}(undef, N1, N2)
        for i = 1:N1
            for j = 1:N2
                d2 = copy(d)
                d2[k1] = kwargs[k1][i]
                d2[k2] = kwargs[k2][j]
                points[i, j] = Point(d2)
            end
        end
        return plot(V, points; filename = filename, open = open)
    end
    error(
        "Can only plot 1- or 2-dimensional value functions. You provided " *
        "$(length(variables)).",
    )
end
