#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export InterpolatedValueFunction, Noise

immutable Noise{T}
    x::T
    probability::Float64
end

function Noise{T}(x::AbstractVector{T})
    y = Noise{T}[]
    for i in x
        push!(y, Noise(i, 1/length(x)))
    end
    y
end

function Noise{T}(x::AbstractVector{T}, p::AbstractVector{Float64})
    y = Noise{T}[]
    for (xi, pi) in zip(x, p)
        push!(y, Noise(xi, pi))
    end
    y
end

function sample{T}(x::Vector{Noise{T}})
    r = rand()
    for i in 1:length(x)
        @inbounds r -= x[i].probability
        if r < eps(Float64)
            return x[i].x
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end

type InterpolatedValueFunction{C<:AbstractCutOracle, T, T2} <: AbstractValueFunction
    initial_price::T
    location::T
    rib_locations::Vector{T}
    variables::Vector{JuMP.Variable}
    cutoracles::Vector{C}
    noises::Vector{Noise{T2}}
    objective::Function
    dynamics::Function
    A::Array{Float64, 2}
end
InterpolatedValueFunction(;
    cut_oracle = DefaultCutOracle(),
    dynamics = (p,w,t,i)->p,
    initial_price = 0.0,
    rib_locations = [0.0, 1.0],
    noise         = Noise([0.0])) = InterpolatedValueFunction(initial_price,
        initial_price, rib_locations,JuMP.Variable[], typeof(cut_oracle)[], noise, (p)->QuadExpr(p), dynamics, Array{Float64}(0,0)
    )

summarise{C,T,T2}(::Type{InterpolatedValueFunction{C,T,T2}}) = "Interpolated Value Function"

function stageobjective!(vf::InterpolatedValueFunction, sp::JuMP.Model, obj::Function)
    vf.objective = obj
end

getstageobjective(vf::InterpolatedValueFunction, sp::JuMP.Model) = getvalue(vf.objective(vf.location))

function setobjective!(sp::JuMP.Model, price, noise)
    vf = valueoracle(sp)
    p = vf.dynamics(price, noise, ext(sp).stage, ext(sp).markovstate)
    vf.location = p
    # stage objective obj
    stageobj = vf.objective(p)

    # future objective
    future_value = interpolate(vf)
    # set
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), stageobj)
    else
        JuMP.setobjective(sp, getsense(sp), stageobj + future_value)
    end
end
function interpolate{C<:AbstractCutOracle, T,T2}(vf::InterpolatedValueFunction{C, T,T2})
    y = AffExpr(0.0)
    _set = false
    upper_idx = length(vf.rib_locations)
    for i in 2:length(vf.rib_locations)
        if vf.location <= vf.rib_locations[i]
            upper_idx = i
            break
        end
    end
    lower_idx = upper_idx - 1
    lambda = (vf.location - vf.rib_locations[lower_idx]) / (vf.rib_locations[upper_idx] - vf.rib_locations[lower_idx])
    append!(y, vf.variables[lower_idx] * (1- lambda))
    append!(y, vf.variables[upper_idx] * lambda)
    y
end

function init!{C<:AbstractCutOracle, T, T2}(vf::InterpolatedValueFunction{C, T, T2}, m::JuMP.Model, sense, bound)
    for r in vf.rib_locations
        push!(vf.variables, futureobjective!(sense, m, bound))
        push!(vf.cutoracles, C())
    end
    vf
end

function passpriceforward!(m::SDDPModel, sp::JuMP.Model)
    stage = ext(sp).stage
    if stage < length(m.stages)
        # pass price forward
        for sp2 in subproblems(m, stage + 1)
            valueoracle(sp2).location = valueoracle(sp).location
        end
    end
end

function solvesubproblem!(::Type{ForwardPass}, vf::InterpolatedValueFunction, m::SDDPModel, sp::JuMP.Model)
    if ext(sp).stage == 1
        vf.location = vf.initial_price
    end
    p = vf.location
    # learn noise
    w = sample(vf.noises)

    # update price
    setobjective!(sp, p, w)
    passpriceforward!(m, sp)
    @assert JuMP.solve(sp) == :Optimal
end

function backwardpass!{C, T, T2}(m::SDDPModel{InterpolatedValueFunction{C, T, T2}}, settings::Settings)
    for t in (nstages(m)-1):-1:1
        for sp in subproblems(m, t)
            vf = valueoracle(sp)
            if t==1
                vf.location = vf.initial_price
            end
            ex = ext(sp)
            i = 0
            for (rib, theta, cutoracle) in zip(vf.rib_locations, vf.variables, vf.cutoracles)
                i += 1
                reset!(m.storage)
                for sp2 in subproblems(m, t+1)
                    vf2 = valueoracle(sp2)
                    ex2 = ext(sp2)
                    markov_prob = getstage(m, ex2.stage).transitionprobabilities[ex.markovstate, ex2.markovstate]
                    for noise in vf2.noises
                        setobjective!(sp2, rib, noise.x)
                        if hasscenarios(sp2)
                            for (scenario, probability) in zip(ex2.scenarios, ex2.scenarioprobability)
                                setscenario!(sp2, scenario)
                                @assert JuMP.solve(sp2) == :Optimal
                                # save
                                push!(m.storage.objective, getobjectivevalue(sp2))
                                push!(m.storage.probability, markov_prob * noise.probability * probability)
                                push!(m.storage.modifiedprobability, 0.0)
                                push!(m.storage.markov, ex2.markovstate)
                                push!(m.storage.duals, zeros(nstates(sp2)))
                                saveduals!(m.storage.duals[end], sp2)
                            end
                        else
                            @assert JuMP.solve(sp2) == :Optimal
                            # save
                            push!(m.storage.objective, getobjectivevalue(sp2))
                            push!(m.storage.probability, markov_prob * noise.probability)
                            push!(m.storage.modifiedprobability, 0.0)
                            push!(m.storage.markov, ex2.markovstate)
                            push!(m.storage.duals, zeros(nstates(sp2)))
                            saveduals!(m.storage.duals[end], sp2)
                        end

                    end
                end
                # add cut
                I = 1:length(m.storage.objective)
                modifyprobability!(ex.riskmeasure,
                    view(m.storage.modifiedprobability.data, I),
                    m.storage.probability.data[I],
                    sp,
                    m.storage.state,
                    m.storage.duals.data[I],
                    m.storage.objective.data[I]
                )
                cut = constructcut(m, sp)
                if settings.cut_output_file != ""
                    writecut!(settings.cut_output_file, cut, ex.stage, ex.markovstate, i)
                end
                storecut!(cutoracle, m, sp, cut)
                storecut!(m, sp, cut, i)
                addcut!(vf, sp, theta, cut)
            end
        end
    end
    reset!(m.storage)
    for sp in subproblems(m, 1)
        vf = valueoracle(sp)
        ex = ext(sp)
        for noise in vf.noises
            setobjective!(sp, vf.initial_price, noise.x)
            if hasscenarios(sp)
                for (scenario, probability) in zip(ex.scenarios, ex.scenarioprobability)
                    setscenario!(sp, scenario)
                    @assert JuMP.solve(sp) == :Optimal
                    push!(m.storage.objective, getobjectivevalue(sp))
                end
            else
                @assert JuMP.solve(sp) == :Optimal
                push!(m.storage.objective, getobjectivevalue(sp))
            end
        end
    end
    mean(m.storage.objective)
end

function writecut!(filename::String, cut::Cut, stage::Int, markovstate::Int, rib::Int)
    open(filename, "a") do file
        write(file, "$(stage), $(markovstate), $(rib), $(cut.intercept)")
        for pi in cut.coefficients
            write(file, ",$(pi)")
        end
        write(file, "\n")
    end
end

function addcut!(vf::InterpolatedValueFunction, sp::JuMP.Model, theta::JuMP.Variable, cut::Cut)
    ex = ext(sp)
    affexpr = AffExpr(cut.intercept)
    for i in 1:nstates(sp)
        append!(affexpr, cut.coefficients[i] * ex.states[i].variable)
    end
    _addcut!(ex.sense, sp, theta, affexpr)
end

function storekey!(::Type{Val{:price}}, store, markov::Int, scenarioidx::Int, sp::JuMP.Model, t::Int)
    push!(store, valueoracle(sp).location)
end

function rebuildsubproblem!{C<:AbstractCutOracle,T,T2}(vf::InterpolatedValueFunction{C,T,T2}, m::SDDPModel, sp::JuMP.Model)
    n = n_args(m.build!)
    ex = ext(sp)
    for i in 1:nstates(sp)
        pop!(ex.states)
    end
    for i in 1:length(ex.scenarios)
        pop!(ex.scenarios)
    end
    sp = Model(solver = m.lpsolver)

    empty!(vf.variables)
    for r in vf.rib_locations
        push!(vf.variables, futureobjective!(getsense(m.sense), sp, ex.problembound))
    end

    sp.ext[:SDDP] = ex
    if n == 2
        m.build!(sp, ex.stage)
    elseif n == 3
        m.build!(sp, ex.stage, ex.markovstate)
    end

    # re-add cuts
    for i in 1:length(vf.variables)
        for cut in validcuts(vf.cutoracles[i])
            addcut!(vf, sp, vf.variables[i], cut)
        end
    end
    m.stages[ex.stage].subproblems[ex.markovstate] = sp
end
rebuildsubproblem!{T,T2}(vf::InterpolatedValueFunction{DefaultCutOracle,T,T2}, m::SDDPModel, sp::JuMP.Model) = nothing
