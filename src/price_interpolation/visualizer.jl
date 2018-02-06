#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using Clp

function initialize_minimizer(p::Vector{Float64}, b::Vector{Float64}, A::Array{Float64, 2}, is_minimization::Bool, pmin::Float64, pmax::Float64, lipschitz_constant::Float64, bound::Float64)
    m = Model(solver=ClpSolver())
    @variable(m, -lipschitz_constant <= mu <= lipschitz_constant)
    @variable(m, phi)
    @variable(m, x[1:size(A, 2)])
    if is_minimization
        @constraint(m, [i=1:size(A, 1)], p[i] * mu + phi >= b[i] + dot(A[i,:], x))
        @constraint(m, pmin * mu + phi >= bound)
        @constraint(m, pmax * mu + phi >= bound)
        @objective(m, Min, 0 * mu + phi)
    else
        @constraint(m, [i=1:size(A, 1)], p[i] * mu + phi <= b[i] + dot(A[i,:], x))
        @constraint(m, pmin * mu + phi <= bound)
        @constraint(m, pmax * mu + phi <= bound)
        @objective(m, Max, 0 * mu + phi)
    end
    @constraint(m, rhs, x .== zeros(size(A, 2)))
    m
end

function fix_state!(m, state_idx, val)
    JuMP.setRHS(m[:rhs][state_idx], val)
end
function get_best_value!(yi::Vector{Float64}, m, args, x)
    for (i, v) in zip(args, x)
        fix_state!(m, i, v)
    end
    solve(m)
    push!(yi, JuMP.getobjectivevalue(m))
end

function get_best_value!(yi::Vector{Float64}, m, args, x, price)
    sense = JuMP.getobjectivesense(m)
    JuMP.setobjective(m, sense, price * m[:mu] + m[:phi])
    get_best_value!(yi, m, args, x)
end


function _processvaluefunctiondata(prices, cuts::Vector{Cut}, minprice::Float64, maxprice::Float64, is_minimization::Bool, lipschitz_constant, bound, states::Union{Float64, AbstractVector{Float64}}...)

    A, b = getAb(cuts)
    m = initialize_minimizer(prices, b, A, is_minimization, minprice, maxprice, lipschitz_constant, bound)
    if length(states) != size(A, 2) + 1
        error("Incorrect number of states specified")
    end
    free_states = Int[]
    for (i, state) in enumerate(states[1:(end-1)])
        if isa(state, Float64)
            fix_state!(m, i, state)
        else
            push!(free_states, i)
        end
    end
    sense = JuMP.getobjectivesense(m)
    yi = Float64[]
    if isa(states[end], Real)
        # price is fixed. Check one or two convex state is free
        @assert 1 <= length(free_states) <= 2
        # fix the objective
        JuMP.setobjective(m, sense, states[end] * m[:mu] + m[:phi])
        # build x array
        if length(free_states) == 1
            # one-dimensional
            x = states[free_states[1]]'
        else
            @assert length(free_states) == 2
            # two-dimensional
            x = mesh(states[free_states[1]], states[free_states[2]])
        end
        for i in 1:size(x, 2)
            get_best_value!(yi, m, free_states, x[:, i])
        end
    else
        # price is free. Check 0 or 1 convex state is free
        @assert 0 <= length(free_states) <= 1
        @assert length(states[end]) > 1
        # build x array
        if length(free_states) == 1
            x = mesh(states[free_states[1]], states[end])
            for i in 1:size(x, 2)
                get_best_value!(yi, m, free_states, x[1,i], x[2,i])
            end
        else
            @assert length(free_states) == 0
            x = states[end]'
            for pi in x
                get_best_value!(yi, m, [], [], pi)
            end
        end

    end
    return x, yi
end
