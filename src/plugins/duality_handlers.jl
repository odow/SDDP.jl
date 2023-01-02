#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
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
 st (x̄, x′, u) in Xᵢ(w) ∩ S
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
        # Attempt to recover by resetting the optimizer and re-solving.
        if JuMP.mode(node.subproblem) != JuMP.DIRECT
            MOI.Utilities.reset_optimizer(node.subproblem)
            optimize!(node.subproblem)
        end
    end
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

duality_log_key(::ContinuousConicDuality) = " "

# =========================== LagrangianDuality ============================== #

"""
    LagrangianDuality(;
        method::LocalImprovementSearch.AbstractSearchMethod =
            LocalImprovementSearch.BFGS(100),
    )

Obtain dual variables in the backward pass using Lagrangian duality.

## Arguments

 * `method`: the `LocalImprovementSearch` method for maximizing the Lagrangian
   dual problem.

## Theory

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w) ∩ S
    x̄ - x == 0          [λ]
```
where `S ⊆ ℝ×ℤ`, we solve the problem `max L(λ)`, where:
```
L(λ) = min Cᵢ(x̄, u, w) + θᵢ - λ' h(x̄)
        st (x̄, x′, u) in Xᵢ(w) ∩ S
```
and where `h(x̄) = x̄ - x`.
"""
mutable struct LagrangianDuality <: AbstractDualityHandler
    method::LocalImprovementSearch.AbstractSearchMethod

    function LagrangianDuality(;
        method = LocalImprovementSearch.BFGS(100),
        kwargs...,
    )
        if length(kwargs) > 0
            @warn(
                "Keyword arguments to LagrangianDuality have changed. " *
                "See the documentation for details.",
            )
        end
        return new(method)
    end
end

function get_dual_solution(node::Node, lagrange::LagrangianDuality)
    undo_relax = _relax_integrality(node)
    optimize!(node.subproblem)
    conic_obj, conic_dual = get_dual_solution(node, ContinuousConicDuality())
    undo_relax()
    s = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? -1 : 1
    N = length(node.states)
    x_in_value = zeros(N)
    λ_star, h_expr, h_k = zeros(N), Vector{AffExpr}(undef, N), zeros(N)
    for (i, (key, state)) in enumerate(node.states)
        x_in_value[i] = JuMP.fix_value(state.in)
        h_expr[i] = @expression(node.subproblem, state.in - x_in_value[i])
        JuMP.unfix(state.in)
        λ_star[i] = conic_dual[key]
    end
    # Check that the conic dual is feasible for the subproblem. Sometimes it
    # isn't if the LP dual solution is slightly infeasible due to numerical
    # issues.
    L_k = _solve_primal_problem(node.subproblem, λ_star, h_expr, h_k)
    if L_k === nothing
        return conic_obj, conic_dual
    end
    L_star, λ_star =
        LocalImprovementSearch.minimize(lagrange.method, λ_star) do x
            L_k = _solve_primal_problem(node.subproblem, x, h_expr, h_k)
            return L_k === nothing ? nothing : (s * L_k, s * h_k)
        end
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, x_in_value[i], force = true)
    end
    λ_solution = Dict{Symbol,Float64}(
        name => λ_star[i] for (i, name) in enumerate(keys(node.states))
    )
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
    if JuMP.termination_status(model) != MOI.OPTIMAL
        JuMP.set_objective_function(model, primal_obj)
        return nothing
    end
    h_k .= -JuMP.value.(h_expr)
    L_λ = JuMP.objective_value(model)
    JuMP.set_objective_function(model, primal_obj)
    return L_λ
end

duality_log_key(::LagrangianDuality) = "L"

# ==================== StrengthenedConicDuality ==================== #

"""
    StrengthenedConicDuality()

Obtain dual variables in the backward pass using strengthened conic duality.

## Theory

Given the problem
```
min Cᵢ(x̄, u, w) + θᵢ
 st (x̄, x′, u) in Xᵢ(w) ∩ S
    x̄ - x == 0          [λ]
```
we first obtain an estimate for `λ` using [`ContinuousConicDuality`](@ref).

Then, we evaluate the Lagrangian function:
```
L(λ) = min Cᵢ(x̄, u, w) + θᵢ - λ' (x̄ - x`)
        st (x̄, x′, u) in Xᵢ(w) ∩ S
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
        λ_k[i] = conic_dual[key]
    end
    lagrangian_obj = _solve_primal_problem(node.subproblem, λ_k, h_expr, h_k)
    for (i, (_, state)) in enumerate(node.states)
        JuMP.fix(state.in, x[i], force = true)
    end
    # If lagrangian_obj is `nothing`, then the primal problem didn't solve
    # correctly, probably because it was unbounded (i.e., the dual was
    # infeasible.) But we got the dual from solving the LP relaxation so it must
    # be feasible! Sometimes however, the dual from the LP solver might be
    # numerically infeasible when solved in the primal. That's a shame :(
    # If so, return the conic_obj instead.
    return something(lagrangian_obj, conic_obj), conic_dual
end

duality_log_key(::StrengthenedConicDuality) = "S"

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
    function BanditDuality(args::AbstractDualityHandler...)
        return new(_BanditArm[_BanditArm(arg, Float64[]) for arg in args], 1)
    end
end

function Base.show(io::IO, handler::BanditDuality)
    print(io, "BanditDuality with arms:")
    for arm in handler.arms
        print(io, "\n * ", arm.handler)
    end
    return
end

function BanditDuality()
    return BanditDuality(ContinuousConicDuality(), StrengthenedConicDuality())
end

function _choose_best_arm(handler::BanditDuality)
    _, index = findmax(
        map(handler.arms) do arm
            return Statistics.mean(arm.rewards) + Statistics.std(arm.rewards)
        end,
    )
    handler.last_arm_index = index
    return handler.arms[index]
end

function _update_rewards(handler::BanditDuality, log::Vector{Log})
    # The bound is monotonic, so instead of worring about whether we are
    # maximizing or minimizing, let's just compute:
    #          |bound_t - bound_{t-1}|
    # reward = -----------------------
    #            time_t - time_{t-1}
    t, t′ = log[end], log[end-1]
    reward = abs(t.bound - t′.bound) / (t.time - t′.time)
    # This check is needed because we should probably keep using the first
    # handler until we start to improve the bound. This can take quite a few
    # iterations in some models. (Until we start to improve, the reward will be
    # zero, so we'd never revisit it.
    const_bound = isapprox(log[1].bound, log[end].bound; atol = 1e-6)
    # To start with, we should add the reward to all arms to construct a prior
    # distribution for the arms. The 10 is somewhat arbitrary.
    if length(log) < 10 || const_bound
        for arm in handler.arms
            push!(arm.rewards, reward)
        end
    else
        push!(handler.arms[handler.last_arm_index].rewards, reward)
    end
    return
end

function prepare_backward_pass(
    model::PolicyGraph,
    handler::BanditDuality,
    options::Options,
)
    if length(options.log) > 1
        _update_rewards(handler, options.log)
    end
    arm = _choose_best_arm(handler)
    return prepare_backward_pass(model, arm.handler, options)
end

function get_dual_solution(node::Node, handler::BanditDuality)
    return get_dual_solution(node, handler.arms[handler.last_arm_index].handler)
end

function duality_log_key(handler::BanditDuality)
    return duality_log_key(handler.arms[handler.last_arm_index].handler)
end
