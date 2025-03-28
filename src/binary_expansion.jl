#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

_bitsrequired(x::T) where {T<:Integer} = floor(T, log(x) * inv(log(2))) + 1

"""
    binexpand(x::Int, maximum::Int)

Returns a vector of binary coefficients for the binary expansion of `x`.
Length of the result is determined by the number of bits required to represent
`maximum` in binary.
"""
function binexpand(x::T, maximum::T) where {T<:Integer}
    if x < 0
        error(
            "Cannot perform binary expansion on a negative number." *
            "Initial values of state variables must be nonnegative.",
        )
    elseif maximum <= 0
        error(
            "Cannot perform binary expansion on zero-length " *
            "vector. Upper bounds of state variables must be positive.",
        )
    end
    y = zeros(T, _bitsrequired(maximum))
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

"""
    binexpand(x::Float64, maximum::Float64, eps::Float64 = 0.1)

Returns a vector of binary coefficients for the binary expansion of `x`.
Length of the result is determined by the number of bits required to represent
`maximum` in binary to precision `eps`.
"""
function binexpand(x::Float64, maximum::Float64, eps::Float64 = 0.1)
    @assert eps > 0
    return binexpand(round(Int, x / eps), round(Int, maximum / eps))
end

"""
    bincontract{T}(y::Vector{T})

For vector `y`, evaluates ∑ᵢ 2ⁱ⁻¹yᵢ.
"""
function bincontract(y::Vector{T}) where {T}
    x = zero(T)
    @inbounds for i in eachindex(y)
        x += 2^(i - 1) * y[i]
    end
    return x
end

"""
    bincontract(y::Vector, eps::Float64)

For vector `y`, evaluates ∑ᵢ 2ⁱ⁻¹yᵢ * `eps`.
"""
bincontract(y::Vector, eps::Float64) = bincontract(y) * eps
