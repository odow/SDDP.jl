#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

function newsolutionstore(X::Vector{Symbol})
    d = Dict(
        :markov         => Int[],
        :scenario       => Int[],
        :obj            => Float64[],
        :stageobjective => Float64[]
    )
    for x in X
        d[x] = Any[]
    end
    d
end

storekey!(s::Symbol, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int) = storekey!(Val{s}, store, markov, scenarioidx, sp, t)

function storekey!(::Type{Val{:markov}}, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    if length(store) < t
        push!(store, markov)
    else
        store[t] = markov
    end
end

function storekey!(::Type{Val{:scenario}}, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    if length(store) < t
        push!(store, scenarioidx)
    else
        store[t] = scenarioidx
    end
end

function storekey!(::Type{Val{:obj}}, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    push!(store, getobjectivevalue(sp))
end

function storekey!(::Type{Val{:stageobjective}}, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    push!(store, getstageobjective(sp))
end

function storekey!(::Type{Val{T}}, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int) where T
    push!(store, getvalue(getvariable(sp, T)))
end

function savesolution!(solutionstore::Dict{Symbol, Any}, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    for (key, store) in solutionstore
        storekey!(key, store, markov, scenarioidx, sp, t)
    end
end
savevaluefunction!(store::Dict{Symbol, Any}, sp::JuMP.Model) = storevaluefunction!(store, valueoracle(sp), sp)
storevaluefunction!(store::Dict{Symbol, Any}, ::DefaultValueFunction{C}, sp::JuMP.Model) where C = nothing

function simulate{C}(m::SDDPModel{DefaultValueFunction{C}},
        variables::Vector{Symbol} = Symbol[];
        scenarios::Vector{Int}    = zeros(Int, length(m.stages)),
        markovstates::Vector{Int} = ones(Int, length(m.stages))
    )
    store = newsolutionstore(variables)
    for t in 1:length(m.stages)
        push!(store[:markov], markovstates[t])
        push!(store[:scenario], scenarios[t])
    end
    obj = forwardpass!(m, Settings(), store)
    store[:objective] = obj
    return store
end

function simulate(m::SDDPModel, N::Int, variables::Vector{Symbol}=Symbol[])
    y = Dict{Symbol, Any}[]
    for i in 1:N
        store = newsolutionstore(variables)
        obj = forwardpass!(m, Settings(), store)
        store[:objective] = obj
        push!(y, store)
    end
    y
end
