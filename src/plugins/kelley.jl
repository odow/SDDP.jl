using GLPK
using LinearAlgebra

function _solve_primal(node::Node, dual_vars::Vector{Float64}, slacks)
    model = node.subproblem
    # Set the Lagrangian the objective in the primal model
    old_obj = JuMP.objective_function(model)

    if JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE
        subgradient = slacks
    else
        subgradient = -slacks
    end

    new_obj = old_obj + dot(dual_vars, subgradient)
    JuMP.set_objective_function(model, new_obj)
    JuMP.optimize!(model)
    # println(model)
    # @show value.(model[:x][1].in)
    # @show value.(model[:x][2].in)
    # @show objective_value(model)

    JuMP.set_objective_function(model, old_obj)
    return (JuMP.objective_value(model), JuMP.value.(subgradient))
end

function _kelley(node::Node, dual_vars::Vector{Float64}, mip_solver::SDDiP)
    model = node.subproblem
    # Assume the model has been solved. Solving the MIP is usually very quick
    # relative to solving for the Lagrangian duals, so we cheat and use the
    # solved model's objective as our bound while searching for the optimal duals
    @assert JuMP.termination_status(model) == MOI.OPTIMAL
    obj = JuMP.objective_value(model)

    old_rhs = Dict()
    slacks = []
    for (name, state) in node.states
        old_rhs[name] = JuMP.fix_value(state.in)
        push!(slacks, state.in - old_rhs[name])
        # push!(old_rhs, JuMP.fix_value(state.in))
        JuMP.unfix(state.in)
        # TODO 0, 1 bounds should just always be there
        JuMP.set_lower_bound(state.in, 0)
        JuMP.set_upper_bound(state.in, 1)
    end


    N = length(dual_vars)
    # storage for the best multiplier found so far
    bestmult = copy(dual_vars)
    # dual problem has the opposite sense to the primal
    dualsense = (JuMP.objective_sense(model) == JuMP.MOI.MIN_SENSE ? JuMP.MOI.MAX_SENSE : JuMP.MOI.MIN_SENSE)
    # gradient
    fdash = zeros(N)

    # The approximate model will be a made from linear hyperplanes
    approx_model = JuMP.Model(mip_solver.optimizer)

    @variables approx_model begin
        θ # objective of the approximate model
        x[i in 1:N] # Lagrangian duals
    end
    JuMP.@objective(approx_model, dualsense, θ)

    if dualsense == MOI.MIN_SENSE
        JuMP.set_lower_bound(θ, obj)
        (best_actual, f_actual, f_approx) = (Inf, Inf, -Inf) # TODO could make first one obj
    else
        JuMP.set_upper_bound(θ, obj)
        (best_actual, f_actual, f_approx) = (-Inf, -Inf, Inf)
    end
    # println(approx_model)

    iteration = 0

    while iteration < mip_solver.max_iter
        iteration += 1
        # Evaluate the real function and a subgradient
        (f_actual, fdash) = _solve_primal(node, dual_vars, slacks)

        # update the model and update best function value so far
        if dualsense == MOI.MIN_SENSE
            # @show dual_vars, fdash
            JuMP.@constraint(approx_model, θ >= f_actual + dot(fdash, x - dual_vars))
            if f_actual <= best_actual
                best_actual = f_actual
                bestmult .= dual_vars
            end
        else
            JuMP.@constraint(approx_model, θ <= f_actual + dot(fdash, x - dual_vars))
            if f_actual >= best_actual
                best_actual = f_actual
                bestmult .= dual_vars
            end
        end
        # println(approx_model)
        # Get a bound from the approximate model
        JuMP.optimize!(approx_model)
        # @show JuMP.value.(approx_model[:x])
        @assert JuMP.termination_status(approx_model) == JuMP.MOI.OPTIMAL
        f_approx = JuMP.objective_value(approx_model)
        # @show f_approx

        # More reliable than checking whether subgradient is zero
        if isapprox(best_actual, f_approx, atol = 1e-6, rtol = 1e-6)
            dual_vars .= bestmult
            if dualsense == JuMP.MOI.MIN_SENSE
                dual_vars .*= -1 # bestmult not the same as getvalue(x), approx_model may have just gotten lucky
            end
            # TODO if done right this may not be necessary
            for (name, state) in node.states
                JuMP.fix(state.in, old_rhs[name], force = true)
            end

            return best_actual
        end

        # next iterate
        dual_vars .= value.(x)
    end
    error("could not solve for Lagrangian duals")

end
