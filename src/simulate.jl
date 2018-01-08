#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

function newsolutionstore(X::Vector{Symbol})
    d = Dict(
        :markov         => Int[],
        :noise       => Int[],
        :obj            => Float64[],
        :stageobjective => Float64[]
    )
    for x in X
        d[x] = Any[]
    end
    d
end

storekey!(s::Symbol, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int) = storekey!(Val{s}, store, markov, noiseidx, sp, t)

function storekey!(::Type{Val{:markov}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    if length(store) < t
        push!(store, markov)
    else
        store[t] = markov
    end
end

function storekey!(::Type{Val{:noise}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    if length(store) < t
        push!(store, noiseidx)
    else
        store[t] = noiseidx
    end
end

function storekey!(::Type{Val{:obj}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    push!(store, getobjectivevalue(sp))
end

function storekey!(::Type{Val{:stageobjective}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    push!(store, getstageobjective(sp))
end

function storekey!{T}(::Type{Val{T}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    if JuMPVERSION < v"0.17"
        push!(store, getvalue(getvariable(sp, T)))
    else
        push!(store, getvalue(sp[T]))
    end
end

function savesolution!(solutionstore::Dict{Symbol, Any}, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    for (key, store) in solutionstore
        storekey!(key, store, markov, noiseidx, sp, t)
    end
end

"""
    simulate(m::SDDPPModel,variables::Vector{Symbol};
        noises::Vector{Int}, markovstates::Vector{Int})

# Description

Perform a historical simulation of the current policy in model  `m`.

`noises` is a vector with one element for each stage giving the index of the
(in-sample) stagewise independent noise to sample in each stage. `markovstates`
is a vector with one element for each stage giving the index of the (in-sample)
markov state to sample in each stage.

# Examples

    simulate(m, [:x, :u], noises=[1,2,2], markovstates=[1,1,2])
"""
function simulate(m::SDDPModel,
        variables::Vector{Symbol} = Symbol[];
        noises::AbstractVector{Int}    = zeros(Int, length(m.stages)),
        markovstates::AbstractVector{Int} = ones(Int, length(m.stages))
    )
    store = newsolutionstore(variables)
    for t in 1:length(m.stages)
        push!(store[:markov], markovstates[t])
        push!(store[:noise], noises[t])
    end
    obj = forwardpass!(m, Settings(), store)
    store[:objective] = obj
    return store
end

"""
    results = simulate(m::SDDPPModel, N::Int, variables::Vector{Symbol})

# Description

Perform `N` Monte-Carlo simulations of the current policy in model `m` saving
the values of the variables named in `variables` at every stage.

`results` is a vector containing a dictionary for each simulation. In addition
to the variables specified in the function call, other special keys are:

 - `:stageobjective` - costs incurred during the stage (not future)
 - `:obj`            - objective of the stage including future cost
 - `:markov`         - index of markov state visited
 - `:noise`          - index of noise visited
 - `:objective`      - Total objective of simulation

All values can be accessed as follows

    results[simulation index][key][stage]

with the exception of `:objective` which is just

    results[simulation index][:objective]

# Examples

    results = simulate(m, 10, [:x, :u])
    results[1][:objective] # objective of simulation 1
    mean(r[:objective] for r in results) # mean objective of the simulations
    results[2][:x][3] # value of :x in stage 3 in second simulation
"""
function simulate(m::SDDPModel, N::Int, variables::Vector{Symbol}=Symbol[])
    map(i->randomsimulation(m, variables), 1:N)
end

function randomsimulation(m::SDDPModel, variables::Vector{Symbol})
    store = newsolutionstore(variables)
    obj = forwardpass!(m, Settings(), store)
    store[:objective] = obj
    return store
end
