#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
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

function get_dual_solution(node::Node, ::Nothing)
    return JuMP.objective_value(node.subproblem), Dict{Symbol,Float64}()
end

# ========================= Continuous relaxation ============================ #

"""
    ContinuousConicDuality(optimizer = nothing)

Compute dual variables in the backward pass using conic duality, relaxing any
binary or integer restrictions as necessary.

## Arguments

 * `optimizer`: if  specified, SDDP.jl will call
   `JuMP.set_optimizer(subproblem, optimizer)` before solving problems on the
   backward pass. Use this option only if your default optimizer does not
   support returning a dual solution after the integrality has been relaxed.

## Example

Train a model using `ContinuousConicDuality` by passing it to the
`duality_handler` keyword argument of [`SDDP.train`](@ref):

```jldoctest
julia> import SDDP, HiGHS, Ipopt

julia> duality_handler = SDDP.ContinuousConicDuality(Ipopt.Optimizer)
SDDP.ContinuousConicDuality{DataType}(Ipopt.Optimizer)

julia> model = SDDP.LinearPolicyGraph(;
           stages = 2,
           lower_bound = 0.0,
           optimizer = HiGHS.Optimizer,
       ) do sp, t
           @variable(sp, x, SDDP.State, Int, initial_value = 0)
           @constraint(sp, x.out >= x.in + 0.5)
           @stageobjective(sp, x.out)
       end;

julia> SDDP.train(
           model;
           duality_handler = SDDP.ContinuousConicDuality(),
           print_level = 0,
       )
```
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
struct ContinuousConicDuality{O} <: AbstractDualityHandler
    optimizer::O

    function ContinuousConicDuality(optimizer = nothing)
        return new{typeof(optimizer)}(optimizer)
    end
end

function get_dual_solution(node::Node, ::ContinuousConicDuality)
    if !_has_dual_solution(node)
        model = node.subproblem.ext[:sddp_policy_graph]
        attempt_numerical_recovery(model, node; require_dual = true)
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

function _relax_integrality(node::Node, optimizer)
    if !node.has_integrality
        return () -> nothing
    elseif optimizer === nothing
        return JuMP.relax_integrality(node.subproblem)
    end
    undo_relax = JuMP.relax_integrality(node.subproblem)
    JuMP.set_optimizer(node.subproblem, optimizer)
    return () -> begin
        JuMP.set_optimizer(node.subproblem, node.optimizer)
        undo_relax()
        return
    end
end

function prepare_backward_pass(
    node::Node,
    handler::ContinuousConicDuality,
    ::Options,
)
    return _relax_integrality(node, handler.optimizer)
end

duality_log_key(::ContinuousConicDuality) = " "

# =========================== LagrangianDuality ============================== #

"""
    LagrangianDuality(
        optimizer = nothing;
        method::LocalImprovementSearch.AbstractSearchMethod =
            LocalImprovementSearch.BFGS(100),
    )

Obtain dual variables in the backward pass using Lagrangian duality.

## Arguments

 * `optimizer`: if  specified, SDDP.jl will call
   `JuMP.set_optimizer(subproblem, optimizer)` before solving problems on the
   backward pass. Use this option only if your default optimizer does not
   support returning a dual solution after the integrality has been relaxed.

 * `method`: the `LocalImprovementSearch` method for maximizing the Lagrangian
   dual problem.

## Example

Train a model using `LagrangianDuality` by passing it to the `duality_handler`
keyword argument of [`SDDP.train`](@ref):

```jldoctest
julia> import SDDP, HiGHS, Ipopt

julia> duality_handler = SDDP.LagrangianDuality(Ipopt.Optimizer)
SDDP.LagrangianDuality{DataType}(SDDP.LocalImprovementSearch.BFGS(100), Ipopt.Optimizer)

julia> model = SDDP.LinearPolicyGraph(;
           stages = 2,
           lower_bound = 0.0,
           optimizer = HiGHS.Optimizer,
       ) do sp, t
           @variable(sp, x, SDDP.State, Int, initial_value = 0)
           @constraint(sp, x.out >= x.in + 0.5)
           @stageobjective(sp, x.out)
       end;

julia> SDDP.train(
           model;
           duality_handler = SDDP.LagrangianDuality(),
           print_level = 0,
       )
```

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
mutable struct LagrangianDuality{O} <: AbstractDualityHandler
    method::LocalImprovementSearch.AbstractSearchMethod
    optimizer::O

    function LagrangianDuality(
        optimizer = nothing;
        method = LocalImprovementSearch.BFGS(100),
        kwargs...,
    )
        if length(kwargs) > 0
            @warn(
                "Keyword arguments to LagrangianDuality have changed. " *
                "See the documentation for details.",
            )
        end
        return new{typeof(optimizer)}(method, optimizer)
    end
end

_sparsify(x::Float64) = ifelse(abs(x) < 1e-15, 0.0, x)

function get_dual_solution(node::Node, lagrange::LagrangianDuality)
    if isempty(node.incoming_state_bounds)
        error(
            "LagrangianDuality requires incoming state bounds to be set. " *
            "Please set the `save_incoming_state_bounds` keyword argument to " *
            "`true` when constructing the `PolicyGraph`.",
        )
    end
    undo_relax = _relax_integrality(node, lagrange.optimizer)
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
        l, u, is_integer = node.incoming_state_bounds[key]
        if l > -Inf
            JuMP.set_lower_bound(state.in, l)
        end
        if u < Inf
            JuMP.set_upper_bound(state.in, u)
        end
        if is_integer
            JuMP.set_integer(state.in)
        end
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
        LocalImprovementSearch.minimize(lagrange.method, λ_star, conic_obj) do x
            L_k = _solve_primal_problem(node.subproblem, x, h_expr, h_k)
            return L_k === nothing ? nothing : (s * L_k, s * h_k)
        end
    for (i, (_, state)) in enumerate(node.states)
        if JuMP.is_integer(state.in)
            JuMP.unset_integer(state.in)
        end
        JuMP.fix(state.in, x_in_value[i]; force = true)
    end
    λ_solution = Dict{Symbol,Float64}(
        k => _sparsify(λ_star[i]) for (i, k) in enumerate(keys(node.states))
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
    StrengthenedConicDuality(optimizer = nothing)

Obtain dual variables in the backward pass using strengthened conic duality.

This method is also known in the literature as Strengthened Benders.

## Arguments

 * `optimizer`: if  specified, SDDP.jl will call
   `JuMP.set_optimizer(subproblem, optimizer)` before solving problems on the
   backward pass. Use this option only if your default optimizer does not
   support returning a dual solution after the integrality has been relaxed.

## Example

Train a model using `StrengthenedConicDuality` by passing it to the
`duality_handler` keyword argument of [`SDDP.train`](@ref):

```jldoctest
julia> import SDDP, HiGHS, Ipopt

julia> duality_handler = SDDP.StrengthenedConicDuality(Ipopt.Optimizer)
SDDP.StrengthenedConicDuality{DataType}(Ipopt.Optimizer)

julia> model = SDDP.LinearPolicyGraph(;
           stages = 2,
           lower_bound = 0.0,
           optimizer = HiGHS.Optimizer,
       ) do sp, t
           @variable(sp, x, SDDP.State, Int, initial_value = 0)
           @constraint(sp, x.out >= x.in + 0.5)
           @stageobjective(sp, x.out)
       end;

julia> SDDP.train(
           model;
           duality_handler = SDDP.StrengthenedConicDuality(),
           print_level = 0,
       )
```

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
mutable struct StrengthenedConicDuality{O} <: AbstractDualityHandler
    optimizer::O

    function StrengthenedConicDuality(optimizer = nothing)
        return new{typeof(optimizer)}(optimizer)
    end
end

function get_dual_solution(node::Node, handler::StrengthenedConicDuality)
    undo_relax = _relax_integrality(node, handler.optimizer)
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
        JuMP.fix(state.in, x[i]; force = true)
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
    BanditDuality(args::AbstractDualityHandler...)

Formulates the problem of choosing a duality handler as a multi-armed bandit
problem. The arms to choose between are given by `args`.

    BanditDuality(optimizer::Any = nothing)

The default implementation of `BanditDuality` that picks between the arms:

 * [`ContinuousConicDuality`](@ref)
 * [`StrengthenedConicDuality`](@ref)

If `optimizer` is specified, SDDP.jl will call
`JuMP.set_optimizer(subproblem, optimizer)` before solving problems on the
backward pass. Use this option only if your default optimizer does not support
returning a dual solution after the integrality has been relaxed.

## Example

Train a model using `BanditDuality` by passing it to the `duality_handler`
keyword argument of [`SDDP.train`](@ref):

```jldoctest
julia> import SDDP, HiGHS, Ipopt

julia> SDDP.BanditDuality(Ipopt.Optimizer)
BanditDuality with arms:
 * SDDP.ContinuousConicDuality{DataType}(Ipopt.Optimizer)
 * SDDP.StrengthenedConicDuality{DataType}(Ipopt.Optimizer)

julia> duality_handler = SDDP.BanditDuality(
           SDDP.ContinuousConicDuality(),
           SDDP.StrengthenedConicDuality(),
           SDDP.LagrangianDuality(),
       )
BanditDuality with arms:
 * SDDP.ContinuousConicDuality{Nothing}(nothing)
 * SDDP.StrengthenedConicDuality{Nothing}(nothing)
 * SDDP.LagrangianDuality{Nothing}(SDDP.LocalImprovementSearch.BFGS(100), nothing)

julia> model = SDDP.LinearPolicyGraph(;
           stages = 2,
           lower_bound = 0.0,
           optimizer = HiGHS.Optimizer,
       ) do sp, t
           @variable(sp, x, SDDP.State, Int, initial_value = 0)
           @constraint(sp, x.out >= x.in + 0.5)
           @stageobjective(sp, x.out)
       end;

julia> SDDP.train(model; duality_handler = duality_handler, print_level = 0)
```

## Theory

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
    logs_seen::Int

    function BanditDuality(args::AbstractDualityHandler...)
        return new(_BanditArm[_BanditArm(arg, Float64[]) for arg in args], 1, 1)
    end
end

function Base.empty!(handler::BanditDuality)
    for arm in handler.arms
        empty!(arm.rewards)
    end
    handler.last_arm_index = 1
    handler.logs_seen = 1
    return
end

function Base.show(io::IO, handler::BanditDuality)
    print(io, "BanditDuality with arms:")
    for arm in handler.arms
        print(io, "\n * ", arm.handler)
    end
    return
end

function BanditDuality(optimizer = nothing)
    return BanditDuality(
        ContinuousConicDuality(optimizer),
        StrengthenedConicDuality(optimizer),
    )
end

function _update_arm(handler::BanditDuality)
    scores = map(handler.arms) do arm
        μ, σ = Statistics.mean(arm.rewards), Statistics.std(arm.rewards)
        # σ may be NaN if there are 0 or 1 observations, or if all observations
        # are the same.
        if isnan(σ)
            return μ
        end
        return μ + σ
    end
    if any(isnan, scores)
        # Some scores may be NaN if there are no observations. Pick an arm
        # randomly.
        index = rand(findall(isnan.(scores)))
        handler.last_arm_index = index
        return
    end
    # Compute softmax
    z = exp.(scores .- maximum(scores))
    z ./= sum(z)
    # Sample arm from softmax
    r = rand()
    index = length(z)
    for i in 1:length(z)
        r -= z[i]
        if r <= 0
            index = i
            break
        end
    end
    handler.last_arm_index = index
    return
end

function _update_rewards(handler::BanditDuality, log::Vector{Log})
    # The bound is monotonic, so instead of worrying about whether we are
    # maximizing or minimizing, let's just compute:
    #          |bound_t - bound_{t-1}|
    # reward = -----------------------
    #            time_t - time_{t-1}
    t, t′ = log[end], log[end-1]
    reward = abs(t.bound - t′.bound) / max(t.time - t′.time, 0.1)
    push!(handler.arms[handler.last_arm_index].rewards, reward)
    return
end

function _is_no_progress(log; kwargs...)
    return length(log) >= 2 &&
           isapprox(log[end].bound, log[end-1].bound; kwargs...)
end

function prepare_backward_pass(
    node::Node,
    handler::BanditDuality,
    options::Options,
)
    log = options.log
    if isempty(log)
        empty!(handler)
    end
    if length(log) > handler.logs_seen
        _update_rewards(handler, log)
        handler.logs_seen = length(log)
        if _is_no_progress(log; atol = 1e-6)
            # The last iteration made no progress. Try the next arm.
            handler.last_arm_index =
                mod1(handler.last_arm_index + 1, length(handler.arms))
        else
            _update_arm(handler)
        end
    end
    arm = handler.arms[handler.last_arm_index]
    return prepare_backward_pass(node, arm.handler, options)
end

function get_dual_solution(node::Node, handler::BanditDuality)
    return get_dual_solution(node, handler.arms[handler.last_arm_index].handler)
end

function duality_log_key(handler::BanditDuality)
    return duality_log_key(handler.arms[handler.last_arm_index].handler)
end

# ==================== FixedDiscreteDuality ==================== #

"""
    FixedDiscreteDuality()

Obtain dual variables in the backward pass by solving the MIP, fixing the
integrality, and then solving the continuous relaxation. The intercept is then
adjusted by evaluating the Lagrangian objective value.

## Example

Train a model using `FixedDiscreteDuality` by passing it to the
`duality_handler` keyword argument of [`SDDP.train`](@ref):

```jldoctest
julia> using SDDP

julia> import HiGHS, Ipopt

julia> model = SDDP.LinearPolicyGraph(;
           stages = 2,
           lower_bound = 0.0,
           optimizer = HiGHS.Optimizer,
       ) do sp, t
           @variable(sp, x, SDDP.State, Int, initial_value = 0)
           @constraint(sp, x.out >= x.in + 0.5)
           @stageobjective(sp, x.out)
       end;

julia> SDDP.train(
           model;
           duality_handler = SDDP.FixedDiscreteDuality(),
           print_level = 0,
       )
```
"""
mutable struct FixedDiscreteDuality{O} <: AbstractDualityHandler
    optimizer::O

    function FixedDiscreteDuality(optimizer = nothing)
        return new{typeof(optimizer)}(optimizer)
    end
end

function _fix_integrality(node::Node, ::Nothing)
    return JuMP.fix_discrete_variables(node.subproblem)
end

function _fix_integrality(node::Node, optimizer)
    undo_fix = JuMP.fix_discrete_variables(node.subproblem)
    JuMP.set_optimizer(node.subproblem, optimizer)
    return () -> begin
        JuMP.set_optimizer(node.subproblem, node.optimizer)
        undo_fix()
        return
    end
end

function get_dual_solution(node::Node, handler::FixedDiscreteDuality)
    if !node.has_integrality
        return get_dual_solution(node, ContinuousConicDuality())
    end
    undo_fix = _fix_integrality(node, handler.optimizer)
    optimize!(node.subproblem)
    _, conic_dual = get_dual_solution(node, ContinuousConicDuality())
    undo_fix()
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
        JuMP.fix(state.in, x[i]; force = true)
    end
    if lagrangian_obj !== nothing
        return lagrangian_obj, conic_dual
    end
    # The conic_dual is infeasible for the Lagrangian. We need to bail with
    # _something_, but we can't use the conic_obj because it might cut off part
    # of the true value function.
    undo_relax = relax_integrality(node.subproblem)
    optimize!(node.subproblem)
    ret = get_dual_solution(node, ContinuousConicDuality())
    undo_relax()
    return ret
end

duality_log_key(::FixedDiscreteDuality) = "F"
