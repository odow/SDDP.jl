struct ValueFunction
    model::JuMP.Model
    theta::JuMP.VariableRef
    states::Dict{Symbol, JuMP.VariableRef}
end

function _add_to_value_function(
    model::JuMP.Model,
    states::Dict{Symbol, JuMP.VariableRef},
    convex_approximation::ConvexApproximation
)
    theta = @variable(model)
    if JuMP.objective_sense(model) == MOI.MIN_SENSE
        JuMP.set_lower_bound(theta, JuMP.lower_bound(convex_approximation.theta))
    else
        JuMP.set_upper_bound(theta, JuMP.upper_bound(convex_approximation.theta))
    end
    for state in keys(convex_approximation.states)
        states[state] = @variable(model)
    end
    for cut in convex_approximation.cut_oracle.cuts
        cut = JuMP.AffExpr(cut.intercept)
        for (key, coef) in cut.coefficients
            if !haskey(x,  key)
                x[key] = JuMP.@variable(model)
            end
            JuMP.add_to_expression!(cut, coef * x[key])
        end
        if JuMP.objective_sense(model) == MOI.MIN_SENSE
            JuMP.@constraint(model, theta >= cut)
        else
            JuMP.@constraint(model, theta <= cut)
        end
    end
    return theta
end

function ValueFunction(b::BellmanFunction, sense::MOI.OptimizationSense)
    model = JuMP.Model()
    JuMP.set_objective_sense(model, sense)
    states = Dict{Symbol, VariableRef}()
    global_theta = _add_to_value_function(model, states, b.global_theta)
    local_thetas = VariableRef[
        _add_to_value_function(model, states, l) for l in b.local_thetas
    ]
    for risk_set in b.risk_set_cuts
        expr = JuMP.@expression(
            model, sum(p * v for (p, v) in zip(risk_set, local_thetas)
        )
        if sense == MOI.MIN_SENSE
            JuMP.@constraint(model, global_theta >= expr)
        else
            JuMP.@constraint(model, global_theta <= expr)
        end
    end
    JuMP.set_objective_function(global_theta)
    return ValueFunction(model, global_theta, states)
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
w.r.t. the state-variables.
"""
function evaluate(
    V::ValueFunction,
    point::Dict{Symbol, Float64},
    objective_state = nothing,
    belief_state = nothing
)
    for (state, val) in point
        JuMP.fix(V.model.states[state], val, force = true)
    end
    JuMP.optimize!(model)
    obj = JuMP.objective_value(model)
    duals = Dict{Symbol, Float64}()
    for (key, var) in V.states
        duals[key] = JuMP.dual(JuMP.FixRef(var))
    end
    return obj, duals
end
