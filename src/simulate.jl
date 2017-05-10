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

function savesolution!(solutionstore::Dict{Symbol, Any}, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    for (key, store) in solutionstore
        if key == :markov
            if length(store) < t
                push!(store, markov)
            else
                store[t] = markov
            end
        elseif key == :scenario
            if length(store) < t
                push!(store, scenarioidx)
            else
                store[t] = scenarioidx
            end
        elseif key == :obj
            push!(store, getobjectivevalue(sp))
        elseif key == :stageobjective
            push!(store, getstageobjective(sp))
        else
            push!(store, getvalue(getvariable(sp, key)))
        end
    end
end

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
