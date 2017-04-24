#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################
function sample(x::AbstractVector{Float64})
    r = rand()
    for i in 1:length(x)
        @inbounds r -= x[i]
        if r < eps(Float64)
            return i
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end

function samplesubproblem(stage::Stage, last_markov_state::Int)
    newidx = sample(stage.transitionprobabilities[last_markov_state, :])
    return newidx, subproblem(stage, newidx)
end

function samplescenario(sp::JuMP.Model)
    scenarioidx = sample(ext(sp).scenarioprobability)
    return scenarioidx, ext(sp).scenarios[scenarioidx]
end

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
            push!(store, getstageobjective(sp))
        else
            push!(store, getvalue(getvariable(sp, key)))
        end
    end
end

solvesubproblem!(direction, valuefunction, m::SDDPModel, sp::JuMP.Model) = JuMP.solve(sp)
solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model) = solvesubproblem!(direction, vftype(sp), m, sp)
hasscenarios(sp::JuMP.Model) = length(ext(sp).scenarios) > 0
function forwardpass!(m::SDDPModel, settings::Settings, solutionstore=nothing)
    last_markov_state = 1
    scenarioidx = 0
    for (t, stage) in enumerate(stages(m))
        # choose markov state
        (last_markov_state, sp) = samplesubproblem(stage, last_markov_state)
        # choose and set RHS scenario
        if hasscenarios(sp)
            (scenarioidx, scenario) = samplescenario(sp)
            setscenario!(sp, scenario)
        end
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
        s = stage(m, t)
        # solve all stage t problems
        m.storage.idx = 0
        for (i, sp) in enumerate(subproblems(s))
            solvesubproblem!(BackwardPass, m, sp)
        end
        # add appropriate cuts
        for sp in subproblems(stage(m, t-1))
            modifyvaluefunction!(m, sp)
        end
    end
end

function iteration!(m::SDDPModel, settings::Settings)
    forwardpass!(m, settings)
    backwardpass!(m, settings)
end

function solve(m::SDDPModel, settings::Settings=Settings())
    iteration!(m, settings)
end
