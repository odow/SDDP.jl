#  Copyright 2017-19, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import LinearAlgebra.dot

# ========================= Fallbacks ======================================== #
update_integrality_handler!(integrality_handler::AbstractIntegralityHandler, ::JuMP.OptimizerFactory, ::Int) = integrality_handler

# ========================= Continuous relaxation ============================ #

"""
    ContinuousRelaxation()

The continuous relaxation integrality handler. Duals are obtained in the
backward pass by solving a continuous relaxation of each subproblem.
Integrality constraints are retained in policy simulation.
"""
struct ContinuousRelaxation <: AbstractIntegralityHandler end

# Requires node.subproblem to have been solved with DualStatus == FeasiblePoint
function get_dual_variables(node::Node, ::ContinuousRelaxation)
    # Note: due to JuMP's dual convention, we need to flip the sign for
    # maximization problems.
    dual_values = Dict{Symbol, Float64}()
    if JuMP.dual_status(node.subproblem) != JuMP.MOI.FEASIBLE_POINT
        write_subproblem_to_file(node, "subproblem", throw_error = true)
    end
    dual_sign = JuMP.objective_sense(node.subproblem) == MOI.MIN_SENSE ? 1.0 : -1.0
    for (name, state) in node.states
        ref = JuMP.FixRef(state.in)
        dual_values[name] = dual_sign * JuMP.dual(ref)
    end
    return dual_values
end

# =========================== SDDiP ========================================== #

"""
    SDDiP(; max_iter::Int = 100)

The SDDiP integrality handler introduced by Zhou, J., Ahmed, S., Sun, X.A. in
Nested Decomposition of Multistage Stochastic Integer Programs with Binary State
Variables (2016).

Calculates duals by solving the Lagrangian dual for each subproblem.
"""
mutable struct SDDiP <: AbstractIntegralityHandler
    max_iter::Int
    optimizer::JuMP.OptimizerFactory
    subgradients::Vector{Float64}
    old_rhs::Vector{Float64}
    best_mult::Vector{Float64}
    slacks::Vector{GenericAffExpr{Float64, VariableRef}}

    function SDDiP(; max_iter::Int = 100)
        integrality_handler = new()
        integrality_handler.max_iter = max_iter
        return integrality_handler
    end
end


function update_integrality_handler!(
    integrality_handler::SDDiP, optimizer::JuMP.OptimizerFactory,
    num_states::Int)
    integrality_handler.optimizer = optimizer
    integrality_handler.subgradients = Vector{Float64}(undef, num_states)
    integrality_handler.old_rhs = similar(integrality_handler.subgradients)
    integrality_handler.best_mult = similar(integrality_handler.subgradients)
    integrality_handler.slacks = Vector{GenericAffExpr{Float64, VariableRef}}(undef, num_states)
    return integrality_handler
end

function get_dual_variables(node::Node, integrality_handler::SDDiP)
    dual_values = Dict{Symbol, Float64}()
    # TODO implement smart choice for initial duals
    dual_vars = zeros(length(node.states))
    solver_obj = JuMP.objective_value(node.subproblem)
    try
        kelley_obj = _kelley(node, dual_vars, integrality_handler)::Float64
        @assert isapprox(solver_obj, kelley_obj, atol = 1e-8, rtol = 1e-8)
    catch e
        write_subproblem_to_file(node, "subproblem", throw_error = false)
        rethrow(e)
    end
    for (i, name) in enumerate(keys(node.states))
        dual_values[name] = -dual_vars[i]
    end
    return dual_values
end

function _solve_primal!(subgradients::Vector{Float64}, node::Node, dual_vars::Vector{Float64}, slacks)
    model = node.subproblem
    old_obj = JuMP.objective_function(model)
    # Set the Lagrangian the objective in the primal model
    fact = (JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE ? 1 : -1)
    new_obj = old_obj + fact * dot(dual_vars, slacks)
    JuMP.set_objective_function(model, new_obj)
    JuMP.optimize!(model)
    lagrangian_obj = JuMP.objective_value(model) # TODO moi/solver issue? gurobi/cplex can't get objective after modification

    # Reset old objective, update subgradients using slack values
    JuMP.set_objective_function(model, old_obj)
    subgradients .= fact .* JuMP.value.(slacks)
    return lagrangian_obj
end

function _kelley(node::Node, dual_vars::Vector{Float64}, integrality_handler::SDDiP)
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
    dualsense = (JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE ? JuMP.MOI.MAX_SENSE : JuMP.MOI.MIN_SENSE)

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
    while iter < integrality_handler.max_iter
        iter += 1
        # Evaluate the real function and a subgradient
        f_actual = _solve_primal!(subgradients, node, dual_vars, integrality_handler.slacks)

        # Update the model and update best function value so far
        if dualsense == MOI.MIN_SENSE
            JuMP.@constraint(approx_model, θ >= f_actual + dot(subgradients, x - dual_vars))
            if f_actual <= best_actual
                best_actual = f_actual
                best_mult .= dual_vars
            end
        else
            JuMP.@constraint(approx_model, θ <= f_actual + dot(subgradients, x - dual_vars))
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
        if isapprox(best_actual, f_approx, atol = 1e-8, rtol = 1e-8)
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
    error("Could not solve for Lagrangian duals.")
end
