#  Copyright 2017-21, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

function _deprecate_integrality_handler()
    return error(
        """
SDDP.jl v0.4.0 introduced a number of breaking changes in how we deal with
binary and integer variables.

## Breaking changes

 * We have renamed `SDDiP` to `LagrangianDuality`.
 * We have renamed `ContinuousRelaxation` to `ConicDuality`.
 * Instead of passing the argument to `PolicyGraph`, you now pass it to
   `train`, e.g., `SDDP.train(model; duality_handler = SDDP.LagrangianDuality())`
 * We no longer turn continuous and integer states into a binary expansion. If
   you want to binarize your states, do it manually.

## Why did we do this?

SDDiP (the algorithm presented in the paper) is really two parts:

 1. If you have an integer program, you can compute the dual of the fishing
    constraint using Lagrangian duality; and
 2. If you have pure binary state variables, then cuts constructed from the
    Lagrangian duals result in an optimal policy.

However, these two points are quite independent. If you have integer or
continuous state variables, you can still use Lagrangian duality!

The new system is more flexible because the duality handler is a property of the
solution process, not the model. This allows us to use Lagrangian duality to
solve any dual problem, and it leaves the decision of binarizing the state
variables up to the user. (Hint: we don't think you should do it!)

## Other additions

We also added support for strengthened Benders cuts, which we call
`SDDP.StrengthenedConicDuality()`.

## Future plans

We have a number of future plans in the works, including better Lagrangian
solution methods and better ways of integrating the different types of duality
handlers (e.g., start with ConicDuality, then shift to StrengthenedConicDuality,
then LagrangianDuality).

If these sorts of things interest you, the code is now much more hackable, so
please reach out or read https://github.com/odow/SDDP.jl/issues/246.

Alternatively, if you have interesting examples using SDDiP that you find are
too slow, please send me the examples so we can use them as benchmarks in future
improvements.
    """,
    )
end

SDDiP(args...; kwargs...) = _deprecate_integrality_handler()

ContinuousRelaxation(args...; kwargs...) = _deprecate_integrality_handler()

function relax_integrality(
    model::PolicyGraph,
    duality_handler::AbstractDualityHandler,
)
    undo = Function[]
    for (_, node) in model.nodes
        if node.has_integrality
            push!(undo, relax_integrality(node, duality_handler))
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

function get_dual_solution(node::Node, ::Nothing)
    return JuMP.objective_value(node.subproblem), Dict{Symbol,Float64}()
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
    if !node.has_integrality
        return () -> nothing
    end
    return JuMP.relax_integrality(node.subproblem)
end

# =========================== LagrangianDuality ============================== #

"""
    LagrangianDuality(;
        iteration_limit::Int = 100,
        atol::Float64 = 1e-8,
        rtol::Float64 = 1e-8,
        optimizer = nothing
    )

Obtain dual variables in the backward pass using Lagrangian duality and Kelley's
cutting plane method.

## Arguments

 * `iteration_limit` controls the maximum number of iterations
 * `atol` and `rtol` are the absolute and relative tolerances used in the
   termination criteria
 * If `optimizer` is `nothing`, use the same solver as the main PolicyGraph.
   Otherwise, pass a new optimizer factory, potentially with different
   tolerances to ensure tighter convergence.
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

For primal minimization problems, we solve:
```
L(λ_k) = max t
          st t <= t′ + h(x̄)' * (λ - λ_k)
             t <= initial_bound
```
For primal maximization problems, we solve:
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

# ======================= StrengthenedConicDuality =========================== #

"""
    StrengthenedConicDuality()

Obtain dual variables in the backward pass using strengthened conic duality.
"""
mutable struct StrengthenedConicDuality <: AbstractDualityHandler end

"""
    get_dual_solution(node::Node, ::StrengthenedConicDuality)

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w)
    x̄ - x == 0          [λ]
```
Return the dual solution `λ` computed using strengthened conic duality.

That is, compute λ using conic duality (ignoring integrality), then evaluate the
Lagrangian function:
```
L(λ) = min Cᵢ(x̄, u, w) + θᵢ - λ' h(x̄)
        st (x̄, x′, u) in Xᵢ(w)
```
where `h(x̄) = x̄ - x` to obtain a better estimate of the intercept.
"""
function get_dual_solution(node::Node, ::StrengthenedConicDuality)
    undo_relax = relax_integrality(node, ConicDuality())
    optimize!(node.subproblem)
    conic_obj, conic_dual = get_dual_solution(node, ConicDuality())
    undo_relax()
    if !node.has_integrality
        return conic_obj, conic_dual  # If we're linear, return this!
    end
    num_states = length(node.states)
    λ_k, h_k, x = zeros(num_states), zeros(num_states), zeros(num_states)
    h_expr = Vector{AffExpr}(undef, num_states)
    for (i, (key, state)) in enumerate(node.states)
        x[i] = JuMP.fix_value(state.in)
        h_expr[i] = @expression(node.subproblem, state.in - x[i])
        JuMP.unfix(state.in)
        JuMP.set_lower_bound(state.in, -1e9)
        JuMP.set_upper_bound(state.in, 1e9)
        λ_k[i] = conic_dual[key]
    end
    lagrangian_obj = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, x[i], force = true)
    end
    return lagrangian_obj, conic_dual
end
