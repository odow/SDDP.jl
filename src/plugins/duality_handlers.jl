#  Copyright 2017-21, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# ========================= General methods ================================== #

function get_duality_handler(subproblem::JuMP.Model)
    return get_node(subproblem).duality_handler
end

function relax_integrality(model::PolicyGraph)
    undo = Function[]
    for (_, node) in model.nodes
        if node.has_integrality
            push!(undo, relax_integrality(node, node.duality_handler))
        end
    end
    function undo_relax()
        for f in undo
            f()
        end
        return
    end
    return undo_relax
end

# ========================= Continuous relaxation ============================ #

"""
    ConicDuality()

Obtain dual variables in the backward pass using conic duality.
"""
struct ConicDuality <: AbstractDualityHandler end

function get_dual_solution(node::Node, ::ConicDuality)
    if JuMP.dual_status(node.subproblem) != JuMP.MOI.FEASIBLE_POINT
        write_subproblem_to_file(
            node,
            "subproblem.mof.json",
            throw_error = true,
        )
    end
    # Note: due to JuMP's dual convention, we need to flip the sign for
    # maximization problems.
    dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1
    λ = Dict{Symbol,Float64}(
        name => dual_sign * JuMP.dual(JuMP.FixRef(state.in)) for
        (name, state) in node.states
    )
    return objective_value(node.subproblem), λ
end

function relax_integrality(node::Node, ::ConicDuality)
    return JuMP.relax_integrality(node.subproblem)
end

# =========================== LagrangianDuality ============================== #

"""
    LagrangianDuality(;
        iteration_limit::Int = 100,
        atol::Float64 = 1e-8,
        rtol::Float64 = 1e-8,
    )

Obtain dual variables in the backward pass using Lagrangian duality.

Kelley's method is used to compute Lagrange multipliers.

`iteration_limit` controls the maximum number of iterations
`atol` and `rtol` are the absolute and relative tolerances used in the
termination criteria.

All state variables are assumed to take nonnegative values only.
"""
mutable struct LagrangianDuality <: AbstractDualityHandler
    iteration_limit::Int
    atol::Float64
    rtol::Float64
    optimizer::Any

    function LagrangianDuality(;
        iteration_limit::Int = 100,
        atol::Float64 = 1e-8,
        rtol::Float64 = 1e-8,
        optimizer = nothing,
    )
        return new(iteration_limit, atol, rtol, optimizer)
    end
end

"""
    get_dual_solution(node::Node, lagrange::LagrangianDuality)

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w)
    x̄ - x == 0          [λ]
```
Return the dual solution `λ` computed using Lagrangian duality.

That is, solve the problem:
```
max L(λ)
```
where:
```
L(λ) = min Cᵢ(x̄, u, w) + θᵢ - λ' h(x̄)
        st (x̄, x′, u) in Xᵢ(w)
```
where `h(x̄) = x̄ - x`.

In the maximization case, the optimization senses are reversed, but the sign of
λ stays the same.
"""
function get_dual_solution(node::Node, lagrange::LagrangianDuality)
    # Assume the model has been solved. Solving the MIP is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal
    # duals.
    @assert JuMP.termination_status(node.subproblem) == MOI.OPTIMAL
    # Query the current MIP solution  here. For an optimal dual, we must have
    # equal objective  values. See the check below.
    primal_obj = JuMP.objective_value(node.subproblem)
    # TODO(odow): check the dual objective value is equal to the primal
    # objective value.
    _, λ_vector = _solve_lagrange_with_kelleys(node, lagrange, primal_obj)
    λ = Dict{Symbol,Float64}(
        name => λ_vector[i] for (i, name) in enumerate(keys(node.states))
    )
    return primal_obj, λ
end

"""
    _solve_primal_problem(
        model::JuMP.Model,
        λ::Vector{Float64},
        h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
        h_k::Vector{Float64},
    )

Solve the problem:
```
L(λ_k) = min/max Cᵢ(x̄, u, w) + θᵢ - λ_k' (x̄ - x_k)
              st (x̄, x′, u) in Xᵢ(w)
```
where `h_expr = x̄ - x_k` and store the solution of `x̄ - x_k` in `h_k`.

By inspection, subgradient of `L(λ_k)` is the slack term `-(x̄ - x_k)`.

!!! warning
    The dual of λ_k is independent of the optimization sense!
"""
function _solve_primal_problem(
    model::JuMP.Model,
    λ::Vector{Float64},
    h_expr::Vector{GenericAffExpr{Float64,VariableRef}},
    h_k::Vector{Float64},
)
    primal_obj = JuMP.objective_function(model)
    JuMP.set_objective_function(
        model,
        @expression(model, primal_obj - λ' * h_expr),
    )
    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) == MOI.OPTIMAL
    h_k .= -JuMP.value.(h_expr)
    L_λ = JuMP.objective_value(model)
    JuMP.set_objective_function(model, primal_obj)
    return L_λ
end

"""
    _solve_lagrange_with_kelleys(
        node::Node,
        lagrange::LagrangianDuality,
        initial_bound::Float64,
    )

We solve the Lagrangian dual problem using Kelley's cutting plane algorithm.

For
```
L(λ_k) = max t
          st t <= t′ + h(x̄)' * (λ - λ_k)
             t <= initial_bound
```

```
L(λ) >= min t
         st t >= t′ + h(x̄)' * (λ - λ_k)
            t >= initial_bound
```
"""
function _solve_lagrange_with_kelleys(
    node::Node,
    lagrange::LagrangianDuality,
    initial_bound::Float64,
)
    num_states = length(node.states)
    incoming_state_value = zeros(num_states)     # The original value of x.
    λ_k = zeros(num_states)                      # The current estimate for λ
    λ_star = zeros(num_states)                   # The best estimate for λ
    h_expr = Vector{AffExpr}(undef, num_states)  # The expression for x̄ - x
    h_k = zeros(num_states)                      # The value of x̄_k - x
    for (i, (_, state)) in enumerate(node.states)
        incoming_state_value[i] = JuMP.fix_value(state.in)
        h_expr[i] =
            @expression(node.subproblem, state.in - incoming_state_value[i])
        JuMP.unfix(state.in)
        JuMP.set_lower_bound(state.in, -1e9)
        JuMP.set_upper_bound(state.in, 1e9)
    end
    primal_sense = JuMP.objective_sense(node.subproblem)
    model = if lagrange.optimizer === nothing
        JuMP.Model(node.optimizer)
    else
        JuMP.Model(lagrange.optimizer)
    end
    @variable(model, λ[1:num_states])
    if primal_sense == MOI.MIN_SENSE
        @variable(model, t <= initial_bound)
        @objective(model, Max, t)
        L_best, t_k = -Inf, initial_bound
    else
        @variable(model, t >= initial_bound)
        @objective(model, Min, t)
        L_best, t_k = Inf, initial_bound
    end
    iter = 0
    while !isapprox(L_best, t_k, atol = lagrange.atol, rtol = lagrange.rtol)
        iter += 1
        if iter > lagrange.iteration_limit
            error("Iteration limit exceeded in Lagrangian subproblem.")
        end
        JuMP.optimize!(model)
        @assert JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(model)
        λ_k .= value.(λ)
        L_k = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
        if primal_sense == MOI.MIN_SENSE
            JuMP.@constraint(model, t <= L_k + h_k' * (λ .- λ_k))
            if L_k >= L_best
                L_best = L_k
                λ_star .= λ_k
            end
        else
            JuMP.@constraint(model, t >= L_k + h_k' * (λ .- λ_k))
            if L_k <= L_best
                L_best = L_k
                λ_star .= λ_k
            end
        end
    end
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, incoming_state_value[i], force = true)
    end
    return L_best, λ_star
end
