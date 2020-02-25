#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
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
    binaries::Vector{Tuple{JuMP.VariableRef,Float64,Float64}},
    integers::Vector{VariableRef},
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

get_integrality_handler(subproblem::JuMP.Model) = get_node(subproblem).integrality_handler

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
    ::ContinuousRelaxation,
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
    dual_values = Dict{Symbol,Float64}()
    if JuMP.dual_status(node.subproblem) != JuMP.MOI.FEASIBLE_POINT
        write_subproblem_to_file(node, "subproblem.mof.json", throw_error = true)
    end
    dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1.0 : -1.0
    for (name, state) in node.states
        ref = JuMP.FixRef(state.in)
        dual_values[name] = dual_sign * JuMP.dual(ref)
    end
    return dual_values
end

function relax_integrality(model::PolicyGraph, ::ContinuousRelaxation)
    binaries = Tuple{JuMP.VariableRef,Float64,Float64}[]
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
    binaries = Tuple{JuMP.VariableRef,Float64,Float64}[]
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
    optimizer::Any
    subgradients::Vector{Float64}
    old_rhs::Vector{Float64}
    best_mult::Vector{Float64}
    slacks::Vector{GenericAffExpr{Float64,VariableRef}}
    atol::Float64
    rtol::Float64

    function SDDiP(; iteration_limit::Int = 100, atol::Float64 = 1e-8, rtol::Float64 = 1e-8)
        integrality_handler = new()
        integrality_handler.iteration_limit = iteration_limit
        integrality_handler.atol = atol
        integrality_handler.rtol = rtol
        return integrality_handler
    end
end

function update_integrality_handler!(
    integrality_handler::SDDiP,
    optimizer::Any,
    num_states::Int,
)
    integrality_handler.optimizer = optimizer
    integrality_handler.subgradients = Vector{Float64}(undef, num_states)
    integrality_handler.old_rhs = similar(integrality_handler.subgradients)
    integrality_handler.best_mult = similar(integrality_handler.subgradients)
    integrality_handler.slacks =
        Vector{GenericAffExpr{Float64,VariableRef}}(undef, num_states)
    return integrality_handler
end

function setup_state(
    subproblem::JuMP.Model,
    state::State,
    state_info::StateInfo,
    name::String,
    ::SDDiP,
)
    if state_info.out.binary
        # Only in this case we treat `state` as a real state variable
        setup_state(subproblem, state, state_info, name, ContinuousRelaxation())
    else
        if !isfinite(state_info.out.upper_bound)
            error("When using SDDiP, state variables require an upper bound.")
        end

        if state_info.out.integer
            # Initial value must be integral
            initial_value = binexpand(
                Int(state_info.initial_value),
                floor(Int, state_info.out.upper_bound),
            )
            num_vars = length(initial_value)

            binary_vars = JuMP.@variable(
                subproblem,
                [i in 1:num_vars],
                base_name = "_bin_" * name,
                SDDP.State,
                Bin,
                initial_value = initial_value[i]
            )

            JuMP.@constraint(
                subproblem,
                state.in == bincontract([binary_vars[i].in for i = 1:num_vars])
            )
            JuMP.@constraint(
                subproblem,
                state.out == bincontract([binary_vars[i].out for i = 1:num_vars])
            )
        else
            epsilon =
                (
                    haskey(state_info.kwargs, :epsilon) ? state_info.kwargs[:epsilon] : 0.1
                )::Float64
            initial_value = binexpand(
                float(state_info.initial_value),
                float(state_info.out.upper_bound),
                epsilon,
            )
            num_vars = length(initial_value)

            binary_vars = JuMP.@variable(
                subproblem,
                [i in 1:num_vars],
                base_name = "_bin_" * name,
                SDDP.State,
                Bin,
                initial_value = initial_value[i]
            )

            JuMP.@constraint(
                subproblem,
                state.in == bincontract([binary_vars[i].in for i = 1:num_vars], epsilon)
            )
            JuMP.@constraint(
                subproblem,
                state.out == bincontract([binary_vars[i].out for i = 1:num_vars], epsilon)
            )
        end
    end
    return
end

function get_dual_variables(node::Node, integrality_handler::SDDiP)
    dual_values = Dict{Symbol,Float64}()
    # TODO implement smart choice for initial duals
    dual_vars = zeros(length(node.states))
    solver_obj = JuMP.objective_value(node.subproblem)
    try
        kelley_obj = _kelley(node, dual_vars, integrality_handler)::Float64
        @assert isapprox(solver_obj, kelley_obj, atol = 1e-8, rtol = 1e-8)
    catch e
        write_subproblem_to_file(node, "subproblem.mof.json", throw_error = false)
        rethrow(e)
    end
    for (i, name) in enumerate(keys(node.states))
        # TODO (maybe) change dual signs inside kelley to match LP duals
        dual_values[name] = -dual_vars[i]
    end
    return dual_values
end

relax_integrality(::PolicyGraph, ::SDDiP) =
    Tuple{JuMP.VariableRef,Float64,Float64}[], JuMP.VariableRef[]

function _solve_primal!(
    subgradients::Vector{Float64},
    node::Node,
    dual_vars::Vector{Float64},
    slacks,
)
    model = node.subproblem
    old_obj = JuMP.objective_function(model)
    # Set the Lagrangian the objective in the primal model
    fact = (JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE ? 1 : -1)
    new_obj = old_obj + fact * LinearAlgebra.dot(dual_vars, slacks)
    JuMP.set_objective_function(model, new_obj)
    JuMP.optimize!(model)
    lagrangian_obj = JuMP.objective_value(model)

    # Reset old objective, update subgradients using slack values
    JuMP.set_objective_function(model, old_obj)
    subgradients .= fact .* JuMP.value.(slacks)
    return lagrangian_obj
end

function _kelley(node::Node, dual_vars::Vector{Float64}, integrality_handler::SDDiP)
    atol = integrality_handler.atol
    rtol = integrality_handler.rtol
    model = node.subproblem
    # Assume the model has been solved. Solving the MIP is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal duals
    @assert JuMP.termination_status(model) == MOI.OPTIMAL
    obj = JuMP.objective_value(model)

    for (i, (name, state)) in enumerate(node.states)
        integrality_handler.old_rhs[i] = JuMP.fix_value(state.in)
        integrality_handler.slacks[i] = state.in - integrality_handler.old_rhs[i]
        JuMP.unfix(state.in)
        JuMP.set_lower_bound(state.in, 0)
        JuMP.set_upper_bound(state.in, 1)
    end

    # Subgradient at current solution
    subgradients = integrality_handler.subgradients
    # Best multipliers found so far
    best_mult = integrality_handler.best_mult
    # Dual problem has the opposite sense to the primal
    dualsense = (
        JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE ? JuMP.MOI.MAX_SENSE :
            JuMP.MOI.MIN_SENSE
    )

    # Approximation of Lagrangian dual as a function of the multipliers
    approx_model = JuMP.Model(integrality_handler.optimizer)

    # Objective estimate and Lagrangian duals
    @variables approx_model begin
        θ
        x[1:length(dual_vars)]
    end
    JuMP.@objective(approx_model, dualsense, θ)

    if dualsense == MOI.MIN_SENSE
        JuMP.set_lower_bound(θ, obj)
        (best_actual, f_actual, f_approx) = (Inf, Inf, -Inf)
    else
        JuMP.set_upper_bound(θ, obj)
        (best_actual, f_actual, f_approx) = (-Inf, -Inf, Inf)
    end

    iter = 0
    while iter < integrality_handler.iteration_limit
        iter += 1
        # Evaluate the real function and a subgradient
        f_actual = _solve_primal!(subgradients, node, dual_vars, integrality_handler.slacks)

        # Update the model and update best function value so far
        if dualsense == MOI.MIN_SENSE
            JuMP.@constraint(
                approx_model,
                θ >= f_actual + LinearAlgebra.dot(subgradients, x - dual_vars)
            )
            if f_actual <= best_actual
                best_actual = f_actual
                best_mult .= dual_vars
            end
        else
            JuMP.@constraint(
                approx_model,
                θ <= f_actual + LinearAlgebra.dot(subgradients, x - dual_vars)
            )
            if f_actual >= best_actual
                best_actual = f_actual
                best_mult .= dual_vars
            end
        end
        # Get a bound from the approximate model
        JuMP.optimize!(approx_model)
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        f_approx = JuMP.objective_value(approx_model)

        # More reliable than checking whether subgradient is zero
        if isapprox(best_actual, f_approx, atol = atol, rtol = rtol)
            dual_vars .= best_mult
            if dualsense == JuMP.MOI.MIN_SENSE
                dual_vars .*= -1
            end
            for (i, (name, state)) in enumerate(node.states)
                JuMP.fix(state.in, integrality_handler.old_rhs[i], force = true)
            end

            return best_actual
        end
        # Next iterate
        dual_vars .= value.(x)
    end
    error("Could not solve for Lagrangian duals. Iteration limit exceeded.")
end
