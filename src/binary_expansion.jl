#  Copyright 2017, Oscar Dowson

export binexpand, bincontract, bitsrequired

# Cache these for speed
const log2inv = 1 / log(2)
# dereferencing an array is faster than a 2^x call
const _2i_ = Int[2^(i-1) for i in 1:floor(Int, log(typemax(Int)) * log2inv)+1]
const _2i_L = length(_2i_)

function binexpand!(y::Vector{Int}, x::Int)
    if x < 0
        error("Values to be expanded must be nonnegative. Currently x = $x.")
    end
    @inbounds for i in length(y):-1:1
        k = _2i_[i]
        if x >= k
            y[i] = 1
            x -= k
        end
    end
    if x > 0
        error("Unable to expand binary. Overflow of $x.")
    end
end

function bitsrequired(x::Int)
    floor(Int, log(x) * log2inv) + 1
end
function bitsrequired(x::Float64, eps::Float64=0.1)
    xx = round(Int, x / eps)
    floor(Int, log(xx) * log2inv) + 1
end

"""
    binexpand(x::Int; length::Int=-1, maximum::Real=-1)

Returns an array of 0/1 coefficients for the binary expansion of `x`.
If trailing zeroes are needed, use `length` to specify the length of the output.
"""
function binexpand(x::Int; length::Int=-1, maximum::Real=-1)
    x < 0 && error("Cannot perform binary expansion on a negative number.")
    if maximum != -1
        if length != -1
            warn("Length is being ignored.")
        end
        length = bitsrequired(floor(Int, maximum))
    end
    if length == -1
        y = zeros(Int, bitsrequired(x))
    else
        y = zeros(Int, length)
    end
    binexpand!(y, x)
    y
end

"""
    binexpand(x::Float64, eps::Float64=0.1; length::Int=-1, maximum::Real=-1)

Returns an array of 0/1 coefficients for the binary expansion of `x`.
If trailing zeroes are needed, use `length` to specify the length of the output.
"""
function binexpand(x::Float64, eps::Float64=0.1; length::Int=-1, maximum::Real=-1)
    x < 0 && error("Cannot perform binary expansion on a negative number.")
    if eps <= 0.0
        error("Epsilon tolerance for Float binary expansion must be strictly greater than 0.")
    end
    xx = round(Int, x / eps)
    if maximum != -1
        if length != -1
            warn("Length is being ignored.")
        end
        length = bitsrequired(floor(Int, maximum / eps))
    end
    binexpand(xx, length=length)
end

function bincontract_2i_(y::Vector{T}) where T
    x = zero(T)
    @inbounds for i in 1:length(y)
        x += _2i_[i] * y[i]
    end
    x
end
function bincontract_pow(y::Vector{T}) where T
    x = zero(T)
    @inbounds for i in 1:length(y)
        x += 2^(i-1) * y[i]
    end
    x
end

"""
    bincontract{T}(y::Vector{T})

For vector `y`, evaluates ∑ᵢ 2ⁱ⁻¹yᵢ.
"""
function bincontract(y::Vector{T}) where T
    if length(y) < _2i_L
        bincontract_2i_(y)
    else
        bincontract_pow(y)
    end
end

function bincontract(::Type{Float64}, y::Vector{T}, eps::Float64=0.1) where T
    if eps <= 0
        error("Epsilon tolerance for Float binary contraction must be strictly greater than 0.")
    end
    xx = bincontract(y)
    xx * eps
end
