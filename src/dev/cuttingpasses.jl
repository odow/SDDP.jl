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
function cutiteration(m::SDDPModel)
    # forward pass
    markov_state = 1
    for (t, stage) in enumerate(stages(m))
        sp = subproblem(m, t, markov_state)
        scenario = sample(ext(sp).scenarioprobability)
        setscenario!(sp, ext(sp).scenarios[scenario])

        status = solve(sp)
        @assert status == :Optimal

        savestates!(m.forwardstorage.states[t], sp)

        markov_state = sample(stage.transitionprobabilities)
    end

    # for t in nstages(m)-:-1:2
    #     # m.storage.state
    #     for i in 1:nsubproblems(m, t)
    #         probability = 1.0
    #         sp = subproblem(m, t, i)
    #         solvesubproblem!(m, sp, probability)
    #     end
    #
    #     for i in 1:nsubproblems(m, t-1)
    #         sp = subproblem(m, t, i)
    #         # change probabilities
    #         modifyprobability!(
    #             ex.risk_measure,
    #             m.storage.newprobabilities,
    #             m.storage.oldprobabilities,
    #             sp,
    #             m.storage.state,
    #             m.storage.duals,
    #             m.storage.objective
    #         )
    #
    #         # construct cut
    #         cut = constructcut(m.storage)
    #     end
    # end
end

function solvesubproblem!{S,V<:DefaultValueFunction,R}(ex::SubproblemExt{S, V, R}, m::SDDPModel, sp::JuMP.Model, probability)

    status = solve(sp)
    @assert status == :Optimal

    # store values for risk measure
    preparestore!(m.storage, nstates(sp), n)
    idx = m.storage.idx
    m.storage.objective[idx]   = getobjectivevalue(sp)
    m.storage.probability[idx] = probability
    saveduals!(m.storage.duals[idx], sp)
    m.storage.idx += 1
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

solvesubproblem!(m::SDDPModel, sp::JuMP.Model, probability) = solvesubproblem!(ext(sp), m, sp, probability)

function preparestore!(s::Storage, l::Int, n::Int)
    append!(s.objective, zeros(max(0, length(s.objective) - s.n)))
    append!(s.probability, zeros(max(0, length(s.probability) - s.n)))
    append!(s.duals, [zeros(l) for i in 1:max(0, length(s.probability) - s.n))])
end
function store!{T}(x::Vector{T}, value::T, idx::Int)
    if length(x) > idx
        x[idx] = value
    else
        push!(x, value)
    end
end

# for stage 1 to t-1
#     solve stage i
#     pass new state forward
#     pass new price forward
#     transition to new markov state
#     record state
# end
