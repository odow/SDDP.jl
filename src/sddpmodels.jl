#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

function _nargs(f::Function)
    @assert length(methods(f)) == 1
    return methods(f).mt.max_args-1
end
getel{A, T <: A}(::Type{A}, x::T, t::Int, i::Int) = x
getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int, i::Int) = x[t]
getel{A, T <: A}(::Type{A}, x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
getel{A, T <: A}(::Type{A}, x::T, t::Int) = x
getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int) = x[t]

function SDDPModel(build!::Function;
    sense           = :Min,
    stages          = 1,
    objective_bound = -1e6,
    markov_states   = 1,
    transition      = Array{Float64}(0,0),
    risk_measure    = Expectation(),
    cut_oracle      = DefaultCutOracle(),
    # value_function = DefaultValueFunction{cut_oracle},
    kwargs...)

    includes_markovstate = _nargs(build!) == 3

    m = SDDPModel()

    for t in 1:stages
        # TODO: Transition for stage
        stage = Stage(transition)
        if !includes_markovstate && (getel(Int, markov_states, t) != 1)
            error("""Because you specified a scenario tree in the SDDPModel constructor, you need to use the

                SDDPModel() do sp, stage, markov_state
                    ... model definition ...
                end

            syntax.""")
        end
        for i in 1:getel(Int, markov_states, t)
            mod = Subproblem(
                stage        = t,
                markov_state = i,
                sense        = optimisationsense(sense),
                bound        = float(objective_bound),
                risk_measure = getel(AbstractRiskMeasure, risk_measure, t, i),
                cut_oracle   = getel(AbstractCutOracle, cut_oracle, t, i),
                value_function = DefaultValueFunction
            )
            if includes_markovstate
                build!(mod, t, i)
            else
                build!(mod, t)
            end
            for i in 1:length(ext(mod).scenarios)
                push!(ext(mod).scenarioprobability, 1 / length(ext(mod).scenarios))
            end
            push!(stage.subproblems, mod)
        end

        push!(m.stages, stage)
    end
    m
end
