#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using Clp

function initialize_minimizer(p::Vector{Float64}, b::Vector{Float64}, A::Array{Float64, 2}, is_minimization::Bool, pmin::Float64, pmax::Float64)
    m = Model(solver=ClpSolver())
    @variable(m, mu)
    @variable(m, phi)
    @variable(m, x[1:size(A, 2)])
    if is_minimization
        @constraint(m, [i=1:size(A, 1)], p[i] * mu + phi >= b[i] + dot(A[i,:], x))
        @constraint(m, pmin * mu + phi >= -1e9)
        @constraint(m, pmax * mu + phi >= -1e9)
        @objective(m, Min, 0 * mu + phi)
    else
        @constraint(m, [i=1:size(A, 1)], p[i] * mu + phi <= b[i] + dot(A[i,:], x))
        @constraint(m, pmin * mu + phi <= 1e9)
        @constraint(m, pmax * mu + phi <= 1e9)
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


function _processvaluefunctiondata(prices, cuts::Vector{Cut}, minprice::Float64, maxprice::Float64, is_minimization::Bool, states::Union{Float64, AbstractVector{Float64}}...)

    A, b = getAb(cuts)
    m = initialize_minimizer(prices, b, A, is_minimization, minprice, maxprice)
    if length(states) != size(A, 2) + 1
        error("Incorrect number of states specified")
    end
    free_args = Int[]
    for (i, state) in enumerate(states[1:(end-1)])
        if isa(state, Float64)
            fix_state!(m, i, state)
        else
            push!(free_args, i)
        end
    end
    sense = JuMP.getobjectivesense(m)
    yi = Float64[]
    if isa(states[end], Real)
        JuMP.setobjective(m, sense, states[end] * m[:mu] + m[:phi])
        # free states
        if length(free_args) == 1
            x = states[free_args[1]]'
        else
            x = mesh(states[free_args[1]], states[free_args[2]])
        end
        for i in 1:size(x, 2)
            get_best_value!(yi, m, free_args, x[:, i])
        end
    else
        # Must be one other free state
        if length(free_args) == 1
            x = mesh(states[free_args[1]], states[end])
        else
            x = zeros(Float64, (0, length(states[end])))
        end
        for i in 1:size(x, 2)
            get_best_value!(yi, m, free_args, x[1, i], x[2, i])
        end
    end
    return x, yi
end
