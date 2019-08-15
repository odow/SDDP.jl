#  Copyright 2017-19, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

const _log2inv = inv(log(2))
_bitsrequired(x::Int) = floor(Int, log(x) * _log2inv) + 1
_bitsrequired(x::Float64, eps::Float64 = 0.1) = _bitsrequired(round(Int, x / eps))

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
    y = zeros(Int, _bitsrequired(floor(Int, maximum)))
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

bincontract(y::Vector, eps::Float64) = bincontract(y) * eps
