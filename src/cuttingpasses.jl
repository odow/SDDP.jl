#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

function samplesubproblem(stage::Stage, last_markov_state::Int)
    newidx = sample(stage.transitionprobabilities[last_markov_state, :])
    return newidx, getsubproblem(stage, newidx)
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

function solvesubproblem!(direction, valuefunction, m::SDDPModel, sp::JuMP.Model)
    @assert JuMP.solve(sp) == :Optimal
end
solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model) = solvesubproblem!(direction, valueoracle(sp), m, sp)
hasscenarios(sp::JuMP.Model) = length(ext(sp).scenarios) > 0
function setstates!(m, sp)
    s = getstage(m, ext(sp).stage-1)
    for (st, v) in zip(states(sp), s.state)
        setvalue!(st, v)
    end
end
function forwardpass!(m::SDDPModel, settings::Settings, solutionstore=nothing)
    last_markov_state = 1
    scenarioidx = 0
    obj = 0.0
    for (t, stage) in enumerate(stages(m))
        # choose markov state
        (last_markov_state, sp) = samplesubproblem(stage, last_markov_state)
        if t > 1
            setstates!(m, sp)
        end

        # choose and set RHS scenario
        if hasscenarios(sp)
            (scenarioidx, scenario) = samplescenario(sp)
            setscenario!(sp, scenario)
        end
        # solve subproblem
        solvesubproblem!(ForwardPass, m, sp)
        # store stage obj
        obj += getstageobjective(sp)
        # store state
        savestates!(stage.state, sp)
        # save solution for simulations (defaults to no-op)
        savesolution!(solutionstore, last_markov_state, scenarioidx, sp)
    end
    return obj
end

function backwardpass!(m::SDDPModel, settings::Settings)
    # walk backward through the stages
    for t in nstages(m):-1:2
        # solve all stage t problems
        reset!(m.storage)
        for sp in subproblems(m, t)
            setstates!(m, sp)
            solvesubproblem!(BackwardPass, m, sp)
        end
        # add appropriate cuts
        for sp in subproblems(m, t-1)
            modifyvaluefunction!(m, sp)
        end
    end

    reset!(m.storage)
    for sp in subproblems(m, 1)
        solvesubproblem!(BackwardPass, m, sp)
    end
    # TODO: improve over just taking mean of first stage subproblems
    bound = mean(m.storage.objective)

    return bound
end

function iteration!(m::SDDPModel, settings::Settings)
    t = time()

    simulation_objective = forwardpass!(m, settings)
    time_forwards = time() - t

    objective_bound = backwardpass!(m, settings)
    time_backwards = time() - time_forwards - t

    return objective_bound, time_backwards, simulation_objective, time_forwards

end
