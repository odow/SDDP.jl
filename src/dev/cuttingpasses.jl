#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################
function sample(x::Vector{Float64})
    r = rand()
    @inbounds for i in 1:length(x)
        r -= x[i]
        if r < eps(Float64)
            return i
        end
    end
    error("x must be a discrete probablity distribution that sums to one.")
end

function samplesubproblem(stage::Stage, last_markov_state::Int)
    newidx = sample(stage.transitionprobabilities[last_markov_state, :])
    return newidx, subproblem(stage, newidx)
end

function samplescenario(sp::JuMP.Model)
    scenarioidx = sample(ext(sp).scenarioprobability)
    return scenairoidx, ext(sp).scenarios[scenarioidx]
end

function newsolutionstore(X::Vector{Symbol})
    d = Dict(
        :markov         = Int[],
        :scenario       = Int[],
        :obj            = Float64[],
        :stageobjective = Float64[]
    )
    for x in X
        d[x] = Any[]
    end
    d
end
savesolution!(solutionstore::Void, markov::Int, scenarioidx::Int, sp::JuMP.Model) = nothing
function savesolution!(solutionstore::Dict{Symbol, Vector}, markov::Int, scenarioidx::Int, sp::JuMP.Model)
    for (key, store) in solutionstore
        if key == :markov
            push!(store, markov)
        elseif key == :scenario
            push!(store, scenarioidx)
        elseif key == :obj
            push!(store, getobjectivevalue(sp))
        elseif key == :stageobjective
            push!(store, getstageobjective(sp)
        else
            push!(store, getvalue(getvariable(sp, key)))
        end
    end
end

function forwardpass!(m::SDDPModel, settings::Settings, solutionstore=nothing)
    last_markov_state = 1
    for (t, stage) in enumerate(stages(m))
        # choose markov state
        (last_markov_state, sp) = samplesubproblem(stage, last_markov_state)
        # choose and set RHS scenario
        (scenarioidx, scenario) = samplescenario(sp)
        setscenario!(sp, scenario)
        # solve subproblem
        @assert solvesubproblem!(ForwardPass, m, sp) == :Optimal
        # store state
        savestates!(stage.state, sp)
        # save solution for simulations (defaults to no-op)
        savesolution!(solutionstore, last_markov_state, scenarioidx, sp)
    end
end

function backwardpass!(m::SDDPModel, settings::Settings)
    # walk backward through the stages
    for t in nstages(m):-1:2
        stage = stage(m, t)
        # solve all stage t problems
        m.storage.idx = 0
        for (i, sp) in enumerate(subproblems(stage))
            solvesubproblem!(BackwardPass, m, sp)
        end
        # add appropriate cuts

    end
end


solvesubproblem!(direction, m::SDDPModel, ex::SubproblemExt, sp::JuMP.Model) = solve(sp)
solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model) = solvesubproblem!(direction, m, ext(sp), sp)

function padstore!(s::BackwardPassStorage, sp::JuMP.Model)
    while s.idx > s.N
        push!(m.storage.objective, 0.0)
        push!(m.storage.scenario, 0)
        push!(m.storage.markovstate, 0)
        push!(m.storage.probability, 0.0)
        push!(m.storage.duals, zeros(nstates(sp)))
        s.N += 1
    end
end
function solvesubproblem!{S,R}(::Type{BackwardPass}, m::SDDPModel, ex::SubproblemExt{S,DefaultValueFunction,R}, sp::JuMP.Model)
    for i in 1:length(ex.scenarioprobaiblity)
        setscenario!(sp, ex.scenarios[i])
        status = solve(sp)
        m.storage.idx += 1
        padstore!(m.storage, sp)
        m.storage.objective[m.storage.idx] = getobjectivevalue(sp)
        m.storage.scenario[m.storage.idx] = i
        m.storage.probability[m.storage.idx] = ex.scenarioprobability[i]
        m.storage.markovstate[m.storage.idx] = ex.markovstate
        saveduals!(m.storage.duals[m.storage.idx], sp)
    end
end

function constructcut(storage)
    # theta <=/>= E[ (y - πᵀx̄) + πᵀx ]
    intercept = 0.0
    coefficients = zeros(nstates(sp))
    for i in 1:storage.n
        intercept += storage.probability[i] * (storage.objective[i] - dot(storage.stage, storage.duals[i]))
        for j in 1:nstates(sp)
            coefficients[j] += storage.probability[i] * storage.duals[i][j]
        end
    end
    Cut(intercept, coefficients)
end
