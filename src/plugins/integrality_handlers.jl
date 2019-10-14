#  Copyright 2017-19, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# ========================= General methods ================================== #

"""
    enforce_integrality(
        binaries::Vector{Tuple{JuMP.VariableRef, Float64, Float64}},
        integers::Vector{VariableRef})

Set all variables in `binaries` to `SingleVariable-in-ZeroOne()`, and all
variables in `integers` to `SingleVariable-in-Integer()`.

See also [`relax_integrality`](@ref).
"""
function enforce_integrality(
    binaries::Vector{Tuple{JuMP.VariableRef, Float64, Float64}},
    integers::Vector{VariableRef}
)
    JuMP.set_integer.(integers)
    for (x, x_lb, x_ub) in binaries
        if isnan(x_lb)
            JuMP.delete_lower_bound(x)
        else
            JuMP.set_lower_bound(x, x_lb)
        end
        if isnan(x_ub)
            JuMP.delete_upper_bound(x)
        else
            JuMP.set_upper_bound(x, x_ub)
        end
        JuMP.set_binary(x)
    end
    return
end

function get_integrality_handler(subproblem::JuMP.Model)
    return get_node(subproblem).integrality_handler
end

# ========================= Continuous relaxation ============================ #

"""
    ContinuousRelaxation()

The continuous relaxation integrality handler. Duals are obtained in the
backward pass by solving a continuous relaxation of each subproblem.
Integrality constraints are retained in policy simulation.
"""
struct ContinuousRelaxation <: AbstractIntegralityHandler end

function setup_state(
    subproblem::JuMP.Model,
    state::State,
    state_info::StateInfo,
    name::String,
    ::ContinuousRelaxation
)
    node = get_node(subproblem)
    sym_name = Symbol(name)
    @assert !haskey(node.states, sym_name)  # JuMP prevents duplicate names.
    node.states[sym_name] = state
    graph = get_policy_graph(subproblem)
    graph.initial_root_state[sym_name] = state_info.initial_value
    return
end

# Requires node.subproblem to have been solved with DualStatus == FeasiblePoint
function get_dual_variables(node::Node, ::ContinuousRelaxation)
    # Note: due to JuMP's dual convention, we need to flip the sign for
    # maximization problems.
    dual_values = Dict{Symbol, Float64}()
    if JuMP.dual_status(node.subproblem) != MOI.FEASIBLE_POINT
        write_subproblem_to_file(node, "subproblem", throw_error = true)
    end
    dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1.0 : -1.0
    for (name, state) in node.states
        ref = JuMP.FixRef(state.in)
        dual_values[name] = dual_sign * JuMP.dual(ref)
    end
    return dual_values
end

function relax_integrality(model::PolicyGraph, ::ContinuousRelaxation)
    binaries = Tuple{JuMP.VariableRef, Float64, Float64}[]
    integers = JuMP.VariableRef[]
    for (key, node) in model.nodes
        bins, ints = _relax_integrality(node.subproblem)
        append!(binaries, bins)
        append!(integers, ints)
    end
    return binaries, integers
end

# Relax all binary and integer constraints in `model`. Returns two vectors:
# the first containing a list of binary variables and previous bounds,
# and the second containing a list of integer variables.
function _relax_integrality(m::JuMP.Model)
    # Note: if integrality restriction is added via @constraint then this loop doesn't catch it.
    binaries = Tuple{JuMP.VariableRef, Float64, Float64}[]
    integers = JuMP.VariableRef[]
    # Run through all variables on model and unset integrality
    for x in JuMP.all_variables(m)
        if JuMP.is_binary(x)
            x_lb, x_ub = NaN, NaN
            JuMP.unset_binary(x)
            # Store upper and lower bounds
            if JuMP.has_lower_bound(x)
                x_lb = JuMP.lower_bound(x)
                JuMP.set_lower_bound(x, max(x_lb, 0.0))
            else
                JuMP.set_lower_bound(x, 0.0)
            end
            if JuMP.has_upper_bound(x)
                x_ub = JuMP.upper_bound(x)
                JuMP.set_upper_bound(x, min(x_ub, 1.0))
            else
                JuMP.set_upper_bound(x, 1.0)
            end
            push!(binaries, (x, x_lb, x_ub))
        elseif JuMP.is_integer(x)
            JuMP.unset_integer(x)
            push!(integers, x)
        end
    end
    return binaries, integers
end

# =========================== SDDiP ========================================== #

"""
    SDDiP(; iteration_limit::Int = 100, atol::Float64, rtol::Float64)

The SDDiP integrality handler introduced by Zou, J., Ahmed, S. & Sun, X.A.
Math. Program. (2019) 175: 461. Stochastic dual dynamic integer programming.
https://doi.org/10.1007/s10107-018-1249-5.

Calculates duals by solving the Lagrangian dual for each subproblem. Kelley's
method is used to compute Lagrange multipliers. `iteration_limit` controls the
maximum number of iterations, and `atol` and `rtol` are the absolute and
relative tolerances used in the termination criteria.

All state variables are assumed to take nonnegative values only.
"""
mutable struct SDDiP <: AbstractIntegralityHandler
    iteration_limit::Int
    atol::Float64
    rtol::Float64
    regularizing_weight::Float64
    # Storage for use in Kelley's to avoid excessive memory allocation.
    optimizer::JuMP.OptimizerFactory
    subgradients::Vector{Float64}
    old_rhs::Vector{Float64}
    best_mult::Vector{Float64}
    slacks::Vector{Tuple{VariableRef, Float64}}

    function SDDiP(
        ;
        iteration_limit::Int = 100,
        atol::Float64 = 1e-8,
        rtol::Float64 = 1e-8,
        regularizing_weight::Float64 = 0.0
    )
        integrality_handler = new()
        integrality_handler.iteration_limit = iteration_limit
        integrality_handler.atol = atol
        integrality_handler.rtol = rtol
        integrality_handler.regularizing_weight = regularizing_weight
        return integrality_handler
    end
end

function update_integrality_handler!(
    integrality_handler::SDDiP,
    optimizer::JuMP.OptimizerFactory,
    num_states::Int
)
    integrality_handler.optimizer = optimizer
    integrality_handler.subgradients = Vector{Float64}(undef, num_states)
    integrality_handler.old_rhs = Vector{Float64}(undef, num_states)
    integrality_handler.best_mult = Vector{Float64}(undef, num_states)
    integrality_handler.slacks = Vector{Tuple{VariableRef, Float64}}(undef, num_states)
    return integrality_handler
end

function setup_state(
    subproblem::JuMP.Model,
    state::State,
    state_info::StateInfo,
    name::String,
    ::SDDiP
)
    if state_info.out.binary
        # Only in this case we treat `state` as a real state variable
        setup_state(subproblem, state, state_info, name, ContinuousRelaxation())
        return
    end
    if !isfinite(state_info.out.upper_bound)
        error("When using SDDiP, state variables require an upper bound.")
    end
    if state_info.out.integer
        epsilon = 1.0
        initial_value = binexpand(
            Int(state_info.initial_value),
            floor(Int, state_info.out.upper_bound)
        )
    else
        epsilon = get(state_info.kwargs, :epsilon, 0.1)
        initial_value = binexpand(
            float(state_info.initial_value),
            float(state_info.out.upper_bound),
            epsilon
        )
    end
    num_vars = length(initial_value)
    binary_vars = JuMP.@variable(
        subproblem,
        [i in 1:num_vars],
        SDDP.State,
        Bin,
        initial_value = initial_value[i],
        base_name = "_bin_" * name
    )
    JuMP.@constraint(
        subproblem,
        state.in == bincontract([binary_vars[i].in for i in 1:num_vars], epsilon)
    )
    JuMP.@constraint(
        subproblem,
        state.out == bincontract([binary_vars[i].out for i in 1:num_vars], epsilon)
    )
    return
end

function get_dual_variables(node::Node, integrality_handler::SDDiP)
    dual_values = Dict{Symbol, Float64}()
    # TODO: implement smart choice for initial duals.
    # One possibility is to start with infeasible duals so that Kelley's works
    # towards a feasible one.
    dual_vars = zeros(length(node.states))
    solver_obj = JuMP.objective_value(node.subproblem)
    try
        kelley_obj = _kelley(node, dual_vars, integrality_handler)::Float64
        if !isapprox(solver_obj, kelley_obj, atol = 1e-4, rtol = 1e-3)
            @warn("""
            SDDiP: objective of dual problem is far from the primal problem.
            Consider tighting the tolerances in SDDiP(; atol, rtol).
            Primal objective: $(solver_obj)
            Dual objective: $(kelley_obj)
            """)
        end
    catch e
        write_subproblem_to_file(node, "subproblem", throw_error = false)
        rethrow(e)
    end
    for (i, name) in enumerate(keys(node.states))
        dual_values[name] = -dual_vars[i]
    end
    return dual_values
end

function relax_integrality(::PolicyGraph, ::SDDiP)
    return Tuple{JuMP.VariableRef, Float64, Float64}[], JuMP.VariableRef[]
end

function _solve_primal(
    subgradients::Vector{Float64},
    node::Node,
    dual_vars::Vector{Float64},
    slacks::Vector{Tuple{VariableRef, Float64}},
    old_obj
)
    N = length(subgradients)
    @assert N == length(dual_vars) == length(slacks)
    sense = JuMP.objective_sense(node.subproblem)
    sign = sense == JuMP.MOI.MIN_SENSE ? 1 : -1
    @objective(
        node.subproblem,
        sense,
        old_obj + sign * sum(
            dual_vars[i] * (slacks[i][1] - slacks[i][2]) for i in 1:N
        )
    )
    JuMP.optimize!(node.subproblem)
    lagrangian_obj = JuMP.objective_value(node.subproblem)
    for i in 1:N
        subgradients[i] = sign * (JuMP.value(slacks[i][1]) - slacks[i][2])
    end
    return lagrangian_obj
end

function _kelley(
    node::Node, dual_vars::Vector{Float64}, sddip::SDDiP
)
    N = length(dual_vars)
    @assert N == length(node.states)
    sddip.best_mult .= NaN
    sddip.subgradients .= NaN

    # Check that the MIP has already been solved. This is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal
    # duals.
    @assert JuMP.termination_status(node.subproblem) == MOI.OPTIMAL

    # Cache some information about the model.
    dual_sense = if JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE
        MOI.MAX_SENSE
    else
        MOI.MIN_SENSE
    end
    sgn = dual_sense == MOI.MIN_SENSE ? 1.0 : -1.0
    primal_bound = JuMP.objective_value(node.subproblem)
    old_objective = JuMP.objective_function(node.subproblem)

    # Collect information about the state variables.
    for (i, (name, state)) in enumerate(node.states)
        sddip.old_rhs[i] = JuMP.fix_value(state.in)
        sddip.slacks[i] = (state.in, sddip.old_rhs[i])
        JuMP.unfix(state.in)
        JuMP.set_lower_bound(state.in, 0.0)
        JuMP.set_upper_bound(state.in, 1.0)
    end

    # Approximation of Lagrangian dual as a function of the multipliers.
    approx_model = JuMP.Model(sddip.optimizer)
    @variables(approx_model, begin
        theta
        pi[1:N]
        pi_p[1:N] >= 0
        pi_n[1:N] >= 0
    end)
    @constraints(approx_model, begin
        reg_p[i = 1:N], pi_p[i] - pi[i] >= 0.0
        reg_n[i = 1:N], pi_n[i] + pi[i] >= 0.0
    end)
    @objective(
        approx_model,
        dual_sense,
        theta + sgn * sddip.regularizing_weight * sum(
                pi_p[i] + pi_n[i] for i in 1:N
            ) / N
        )
    if dual_sense == MOI.MIN_SENSE
        JuMP.set_lower_bound(theta, primal_bound)
        (best_actual, f_actual, f_approx) = (Inf, Inf, -Inf)
    else
        JuMP.set_upper_bound(theta, primal_bound)
        (best_actual, f_actual, f_approx) = (-Inf, -Inf, Inf)
    end

    for iter in 1:sddip.iteration_limit
        # Evaluate the real function and a subgradient. This will update the
        # subgradients in `sddip.subgradients`!
        f_actual = _solve_primal(
            sddip.subgradients, node, dual_vars, sddip.slacks, old_objective
        )

        # Update the model and update best function value so far.
        if dual_sense == MOI.MIN_SENSE
            @constraint(
                approx_model,
                theta >= f_actual + sum(
                    sddip.subgradients[j] * (pi[j] - dual_vars[j]) for j in 1:N
                )
            )
            if f_actual <= best_actual
                best_actual = f_actual
                sddip.best_mult .= dual_vars
            end
        else
            @constraint(
                approx_model,
                theta <= f_actual + sum(
                    sddip.subgradients[j] * (pi[j] - dual_vars[j]) for j in 1:N
                )
            )
            if f_actual >= best_actual
                best_actual = f_actual
                sddip.best_mult .= dual_vars
            end
        end

        # Get a bound from the approximate model.
        # TODO(odow): regularizing around the best multiplier doesn't seem to
        # work. In fact, with the air_conditioning problem, it hurts even when
        # the regularizer_weight=0. Investigate why.
        # for i = 1:N
        #     # n.b.: even on the first iteration, we should have a non-NaN value
        #     # for the best multipliers.
        #     JuMP.set_normalized_rhs(reg_p[i], -sddip.best_mult[i])
        #     JuMP.set_normalized_rhs(reg_n[i], sddip.best_mult[i])
        # end
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == MOI.OPTIMAL
        f_approx = JuMP.value(theta)

        # More reliable than checking whether subgradient is zero, check that
        # the approximated objective is the same as the actual objective.
        if isapprox(best_actual, f_approx, atol = sddip.atol, rtol = sddip.rtol)
            # Before returning, clean up the problem.
            dual_vars .= -sgn .* sddip.best_mult
            for (i, (name, state)) in enumerate(node.states)
                JuMP.fix(state.in, sddip.old_rhs[i], force = true)
            end
            JuMP.set_objective_function(node.subproblem, old_objective)
            return best_actual
        end

        # Next iterate.
        dual_vars .= value.(pi)
    end

    # Make sure to restore the objective before throwing an error to prevent
    # corrupting the user's model.
    for (i, (name, state)) in enumerate(node.states)
        JuMP.fix(state.in, sddip.old_rhs[i], force = true)
    end
    JuMP.set_objective_function(node.subproblem, old_objective)

    error("""
    Could not solve for Lagrangian duals. Iteration limit exceeded. A likely
    cause for this is degeneracy. Consider reformulating the model, or
    increasing the iteration limit in SDDiP(; iteration_limit).
    """)
end
