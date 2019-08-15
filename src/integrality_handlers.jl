#  Copyright 2017-19, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    AbstractIntegralityHandler

The abstract type for the integrality handlers interface.
"""
abstract type AbstractIntegralityHandler end

"""
    ContinuousRelaxation

The continuous relaxation integrality handler. Duals are obtained in the
backward pass by solving a continuous relaxation for each subproblem.
Integrality constraints are retained in policy simulation.
"""
struct ContinuousRelaxation <: AbstractIntegralityHandler end

"""
    SDDiP

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

update_integrality_handler!(integrality_handler::AbstractIntegralityHandler, ::JuMP.OptimizerFactory, ::Int) = integrality_handler

function update_integrality_handler!(integrality_handler::SDDiP, optimizer::JuMP.OptimizerFactory, num_states::Int)
    integrality_handler.optimizer = optimizer
    integrality_handler.subgradients = Vector{Float64}(undef, num_states)
    integrality_handler.old_rhs = similar(integrality_handler.subgradients)
    integrality_handler.best_mult = similar(integrality_handler.subgradients)
    integrality_handler.slacks = Vector{GenericAffExpr{Float64, VariableRef}}(undef, num_states)
    return integrality_handler
end

const _log2inv = inv(log(2))

function binexpand!(y::Vector{Int}, x::Int)
    if x < 0
        error("Values to be expanded must be nonnegative. Currently x = $x.")
    end
    @inbounds for i in length(y):-1:1
        k = 2^(i - 1)
        if x >= k
            y[i] = 1
            x -= k
        end
    end
    if x > 0
        error("Unable to expand binary. Overflow of $x.")
    end
    return y
end

bitsrequired(x::Int) = floor(Int, log(x) * _log2inv) + 1
bitsrequired(x::Float64, eps::Float64 = 0.1) = floor(Int, log(round(Int, x / eps)) * _log2inv) + 1

initial_value_err = "Cannot perform binary expansion on a negative number." *
    "Initial values of state variables must be nonnegative."
upper_bound_err =  "Cannot perform binary expansion on zero-length vector." *
    "Upper bounds of state variables must be positive."

"""
    binexpand(x::Int, maximum)

Returns a vector of binary coefficients for the binary expansion of `x`.
Length of the result is determined by the number of bits required to represent
`maximum` in binary.
"""
function binexpand(x::Int, maximum) # TODO type maximum- can be float or int
    x < 0 && throw(initial_value_err)
    maximum <= 0 && throw(upper_bound_err)
    y = zeros(Int, bitsrequired(floor(Int, maximum)))
    return binexpand!(y, x)
end

"""
    binexpand(x::Float64, maximum, eps::Float64 = 0.1)

Returns a vector of binary coefficients for the binary expansion of `x`.
Length of the result is determined by the number of bits required to represent
`maximum` in binary.
"""
function binexpand(x::Float64, maximum, eps::Float64 = 0.1)
    @assert eps > 0
    return binexpand(round(Int, x / eps), round(Int, maximum / eps))
end

"""
    bincontract{T}(y::Vector{T})

For vector `y`, evaluates ∑ᵢ 2ⁱ⁻¹yᵢ.
"""
function bincontract(y::Vector{T}) where {T}
    if length(y) > floor(Int, log(typemax(Int)) * _log2inv) + 1
        error("Overflow of input of length $(length(y)).")
    end
    x = zero(T)
    @inbounds for i in eachindex(y)
        x += 2^(i - 1) * y[i]
    end
    return x
end

function bincontract(::Type{Float64}, y::Vector{T}, eps::Float64) where {T}
    @assert eps > 0
    return bincontract(y) * eps
end
