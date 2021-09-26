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
 * We have renamed `ContinuousRelaxation` to `ContinuousConicDuality`.
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
handlers (e.g., start with ContinuousConicDuality, then shift to
StrengthenedConicDuality, then LagrangianDuality).

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

function prepare_backward_pass(
    model::PolicyGraph,
    duality_handler::AbstractDualityHandler,
    options::Options,
)
    undo = Function[]
    for (_, node) in model.nodes
        push!(undo, prepare_backward_pass(node, duality_handler, options))
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
    ContinuousConicDuality()

Compute dual variables in the backward pass using conic duality, relaxing any
binary or integer restrictions as necessary.

## Theory

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w) ∪ S
    x̄ - x == 0          [λ]
```
where `S ⊆ ℝ×ℤ`, we relax integrality and using conic duality to solve for `λ`
in the problem:
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w)
    x̄ - x == 0          [λ]
```
"""
struct ContinuousConicDuality <: AbstractDualityHandler end

function get_dual_solution(node::Node, ::ContinuousConicDuality)
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

function _relax_integrality(node::Node)
    if !node.has_integrality
        return () -> nothing
    end
    return JuMP.relax_integrality(node.subproblem)
end

function prepare_backward_pass(node::Node, ::ContinuousConicDuality, ::Options)
    return _relax_integrality(node)
end

# =========================== LagrangianDuality ============================== #

"""
    LagrangianDuality(;
        iteration_limit::Int = 100,
        atol::Float64 = 1e-8,
        rtol::Float64 = 1e-8,
        optimizer = nothing,
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

## Theory

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w) ∪ S
    x̄ - x == 0          [λ]
```
where `S ⊆ ℝ×ℤ`, we solve the problem `max L(λ)`, where:
```
L(λ) = min Cᵢ(x̄, u, w) + θᵢ - λ' h(x̄)
        st (x̄, x′, u) in Xᵢ(w) ∪ S
```
and where `h(x̄) = x̄ - x`.

In the maximization case, the optimization senses are reversed, but the sign of
λ stays the same.

The Lagrangian problem is computed using Kelleys cutting plane method.

We solve:
```
L(λ) <= max t
         st t <= s * (L(λ_k) + h(x̄_k)' * (λ - λ_k)) for k=1,...
            t <= s * initial_bound
```
where `s=1` for primal minimization problems and `s=-1` for primal maximization
problems.

To generate a new cut we solve the cutting plane problem to generate an estimate
`λ_k`. Then we solve the primal problem:
```
L(λ_k) = min/max Cᵢ(x̄, u, w) + θᵢ - λ_k' (x̄ - x_k)
              st (x̄, x′, u) in Xᵢ(w) ∪ S
```
using our estimate `λ_k`.

By inspection, subgradient of `L(λ_k)` is the slack term `-(x̄ - x_k)`.

We converge once `L(λ_k) ≈ t` to the tolerance given by `atol` and `rtol`.

This gives us an optimal dual solution `λ_k`. However, since this will be used
in a cut for SDDP, we can go a step further and attempt to find a dual solution
that is "flat" by solving:
```
min e
 st e >= ‖λ‖
    t <= s * (L(λ_k) + h(x̄_k)' * (λ - λ_k)) for k=1,...
    t <= s * initial_bound
    t >= t_star
```
The L1-norm is implemented as the standard abs-value reformulation:
```
min Σλ⁺ + Σλ⁻
 st t <= s * (L(λ_k) + h(x̄_k)' * ([λ⁺ - λ⁻] - λ_k)) for k=1,...
    t <= s * initial_bound
    t >= t_star
    λ⁺, λ⁻ >= 0
```

If we hit the iteration limit, then we terminate because something probably went
wrong.
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

function get_dual_solution(node::Node, lagrange::LagrangianDuality)
    # Assume the model has been solved. Solving the MIP is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal
    # duals.
    @assert JuMP.termination_status(node.subproblem) == MOI.OPTIMAL
    # Query the current MIP solution  here. This is used as a bound for the
    # cutting plane method.
    primal_obj = JuMP.objective_value(node.subproblem)
    # A sign bit that is used to avoid if-statements in the models.
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1 : -1
    # Storage for the cutting plane method.
    num_states = length(node.states)
    x_in_value = zeros(num_states)               # The original value of x.
    λ_k = zeros(num_states)                      # The current estimate for λ
    λ_star = zeros(num_states)                   # The best estimate for λ
    # The best estimate for the dual objective value, ignoring optimization
    # sense (bigger is better).
    L_star = -Inf
    h_expr = Vector{AffExpr}(undef, num_states)  # The expression for x̄ - x
    h_k = zeros(num_states)                      # The value of x̄_k - x
    # Start by relaxing the fishing constraint.
    for (i, (_, state)) in enumerate(node.states)
        # We're going to need this value later when we reset things.
        x_in_value[i] = JuMP.fix_value(state.in)
        h_expr[i] = @expression(node.subproblem, state.in - x_in_value[i])
        # Relax the constraint from the problem.
        JuMP.unfix(state.in)
        # We need bounds to ensure that the dual problem is feasible. However,
        # they can't be too tight. Let's use 1e9 as a default...
        lb = has_lower_bound(state.out) ? lower_bound(state.out) : -1e9
        ub = has_upper_bound(state.out) ? upper_bound(state.out) : 1e9
        JuMP.set_lower_bound(state.in, lb)
        JuMP.set_upper_bound(state.in, ub)
    end
    # Create the model for the cutting plane algorithm
    model = JuMP.Model(something(lagrange.optimizer, node.optimizer))
    @variable(model, λ⁺[1:num_states] >= 0)
    @variable(model, λ⁻[1:num_states] >= 0)
    @variable(model, t <= s * primal_obj)
    @expression(model, λ, λ⁺ .- λ⁻)
    @objective(model, Max, t)
    # Step 1: find an optimal dual solution and corresponding objective value.
    iter, t_k = 0, s * primal_obj
    while !isapprox(L_star, t_k, atol = lagrange.atol, rtol = lagrange.rtol)
        iter += 1
        if iter > lagrange.iteration_limit
            error("Iteration limit exceeded in Lagrangian subproblem.")
        end
        JuMP.optimize!(model)
        @assert JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        t_k = JuMP.objective_value(model)
        λ_k .= value.(λ)
        L_k = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
        JuMP.@constraint(model, t <= s * (L_k + h_k' * (λ .- λ_k)))
        if s * L_k >= L_star
            L_star = s * L_k
            λ_star .= λ_k
        end
    end
    # Step 2: given the optimal dual objective value, try to find the optimal
    # dual solution with the smallest L1-norm ‖λ‖.
    @objective(model, Min, sum(λ⁺) + sum(λ⁻))
    set_lower_bound(t, t_k)
    # The worst-case scenario in this for-loop is that we run through the
    # iterations without finding a new dual solution. However if that happens
    # we can just keep our current λ_star.
    for _ in (iter+1):lagrange.iteration_limit
        JuMP.optimize!(model)
        @assert JuMP.termination_status(model) == JuMP.MOI.OPTIMAL
        λ_k .= value.(λ)
        L_k = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
        if isapprox(L_star, L_k, atol = lagrange.atol, rtol = lagrange.rtol)
            # At this point we tried the smallest ‖λ‖ from the cutting plane
            # problem, and it returned the optimal dual objective value. No
            # other optimal dual vector can have a smaller norm.
            λ_star = λ_k
            break
        end
        JuMP.@constraint(model, t <= s * (L_k + h_k' * (λ .- λ_k)))
    end
    # Restore the fishing constraint x.in == x_in_value
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, x_in_value[i], force = true)
    end
    λ_solution = Dict{Symbol,Float64}(
        name => λ_star[i] for (i, name) in enumerate(keys(node.states))
    )
    # Remember to correct the sign of the optimal dual objective value.
    return s * L_star, λ_solution
end

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

# ==================== StrengthenedConicDuality ==================== #

"""
    StrengthenedConicDuality()

Obtain dual variables in the backward pass using strengthened conic duality.

## Theory

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w) ∪ S
    x̄ - x == 0          [λ]
```
we first obtain an estiamte for `λ` using [`ContinuousConicDuality`](@ref).

Then, we evaluate the Lagrangian function:
```
L(λ) = min Cᵢ(x̄, u, w) + θᵢ - λ' (x̄ - x`)
        st (x̄, x′, u) in Xᵢ(w) ∪ S
```
to obtain a better estimate of the intercept.
"""
mutable struct StrengthenedConicDuality <: AbstractDualityHandler end

function get_dual_solution(node::Node, ::StrengthenedConicDuality)
    undo_relax = _relax_integrality(node)
    optimize!(node.subproblem)
    conic_obj, conic_dual = get_dual_solution(node, ContinuousConicDuality())
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
        # We need bounds to ensure that the dual problem is feasible. However,
        # they can't be too tight. Let's use 1e9 as a default...
        lb = has_lower_bound(state.out) ? lower_bound(state.out) : -1e9
        ub = has_upper_bound(state.out) ? upper_bound(state.out) : 1e9
        JuMP.set_lower_bound(state.in, lb)
        JuMP.set_upper_bound(state.in, ub)
        λ_k[i] = conic_dual[key]
    end
    lagrangian_obj = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, x[i], force = true)
    end
    return lagrangian_obj, conic_dual
end

# ============================== BanditDuality =============================== #

mutable struct _BanditArm{T}
    handler::T
    rewards::Vector{Float64}
end

"""
    BanditDuality()

Formulates the problem of choosing a duality handler as a multi-armed bandit
problem. The arms to choose between are:

 * [`ContinuousConicDuality`](@ref)
 * [`StrengthenedConicDuality`](@ref)
 * [`LagrangianDuality`](@ref)

Our problem isn't a typical multi-armed bandit for a two reasons:

 1. The reward distribution is non-stationary (each arm converges to 0 as it
    keeps getting pulled.
 2. The distribution of rewards is dependent on the history of the arms that
    were chosen.

We choose a very simple heuristic: pick the arm with the best mean + 1 standard
deviation. That should ensure we consistently pick the arm with the best
likelihood of improving the value function.

In future, we should consider discounting the rewards of earlier iterations, and
focus more on the more-recent rewards.
"""
mutable struct BanditDuality <: AbstractDualityHandler
    arms::Vector{_BanditArm}
    last_arm_index::Int
    iteration::Int
    function BanditDuality(;kwargs...)
        return new(
            _BanditArm[
                _BanditArm(ContinuousConicDuality(), Float64[]),
                _BanditArm(StrengthenedConicDuality(), Float64[]),
                # _BanditArm(LagrangianDuality(; kwargs...), Float64[]),
            ],
            1,
            1,
        )
    end
end

function _reward(arm::_BanditArm)
    return Statistics.mean(arm.rewards) + Statistics.std(arm.rewards)
end

function _update_best_arm(handler::BanditDuality, options::Options)
    # The bound is monotonic, so instead of worring about whether we are
    # maximizing or minimizing, let's just compute:
    #          |bound_t - bound_{t-1}|
    # reward = -----------------------
    #            time_t - time_{t-1}
    t, t′ = options.log[end], options.log[end-1]
    reward = abs(t.bound - t′.bound) / (t.time - t′.time)
    # This check is needed because we should probably keep using the first
    # handler until we start to improve the bound. This can take quite a few
    # iterations in some models. (Until we start to improve, the reward will be
    # zero, so we'd never revisit it.
    const_bound =
        isapprox(options.log[1].bound, options.log[end].bound; atol = 1e-6)
    # To start with, we should add the reward to all arms to construct a prior
    # distribution for the arms. The 10 is somewhat arbitrary.
    if length(options.log) < 10 || const_bound
        for arm in handler.arms
            push!(arm.rewards, reward)
        end
    else
        push!(handler.arms[handler.last_arm_index].rewards, reward)
    end
    # Now pick the best arm ...
    _, index = findmax([_reward(arm) for arm in handler.arms])
    handler.last_arm_index = index
    return
end

function prepare_backward_pass(
    node::Node,
    handler::BanditDuality,
    options::Options,
)
    if handler.iteration < length(options.log)
        _update_best_arm(handler, options)
        handler.iteration += 1
    end
    return prepare_backward_pass(
        node,
        handler.arms[handler.last_arm_index].handler,
        options,
    )
end

function get_dual_solution(node::Node, handler::BanditDuality)
    return get_dual_solution(
        node,
        handler.arms[handler.last_arm_index].handler,
    )
end
