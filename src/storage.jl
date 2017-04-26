#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

mutable struct CachedVector{T} <: AbstractArray{T, 1}
    data::Vector{T}
    n::Int
end

Base.size{T}(x::CachedVector{T}) = (x.n,)
function Base.getindex{T}(x::CachedVector{T}, i)
    if i > x.n
        throw(BoundsError(x, i))
    end
    x.data[i]
end

function Base.setindex!{T}(x::CachedVector{T}, y::T, i)
    if i > x.n
        throw(BoundsError(x, i))
    end
    x.data[i] = y
end
function Base.push!{T}(x::CachedVector{T}, y::T)
    x.n += 1
    if x.n <= length(x.data)
        x.data[x.n] = y
    else
        push!(x.data, y)
    end
end
CachedVector(T) = CachedVector{T}(T[], 0)

function reset!{T}(x::CachedVector{T})
    x.n = 0
end
