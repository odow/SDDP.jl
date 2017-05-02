#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export InterpolatedValueFunction, Noise

struct Noise{T}
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

function sample{T}(x::Vector{Noise{T}})
    r = rand()
    for i in 1:length(x)
        @inbounds r -= x[i].probability
        if r < eps(Float64)
            return i
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end

mutable struct InterpolatedValueFunction{C<:AbstractCutOracle, T} <: AbstractValueFunction
    location::Float64
    rib_locations::Vector{Float64}
    variables::Vector{JuMP.Variable}
    cutoracles::Vector{C}
    noises::Vector{Noise{T}}
    objective::Function
    dynamics::Function
end
InterpolatedValueFunction(;
    cut_oracle = DefaultCutOracle(),
    dynamics = (p,w,t,i)->p,
    initial_price = 0.0,
    rib_locations = [0.0, 1.0],
    noise         = Noise([0.0])) = InterpolatedValueFunction(
        initial_price, rib_locations,JuMP.Variable[], typeof(cut_oracle)[], noise, (p)->QuadExpr(p), dynamics
    )

function stageobjective!{C<:AbstractCutOracle, T}(vf::InterpolatedValueFunction{C, T}, sp::JuMP.Model, obj::Function)
    vf.objective = obj
end

getstageobjective{C<:AbstractCutOracle, T}(vf::InterpolatedValueFunction{C, T}, sp::JuMP.Model) = getvalue(vf.objective(vf.location))

function setobjective!(sp::JuMP.Model, price::Float64, noise)
    vf = valueoracle(sp)
    p = vf.dynamics(price, noise)
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
function interpolate{C<:AbstractCutOracle, T}(vf::InterpolatedValueFunction{C, T})
    y = AffExpr(0.0)
    _set = false
    for i in 2:length(vf.rib_locations)
        if vf.location <= vf.rib_locations[i]
            lambda = (vf.location - vf.rib_locations[i-1]) / (vf.rib_locations[i] - vf.rib_locations[i-1])
            append!(y, vf.variables[i-1] * lambda)
            append!(y, vf.variables[i] * (1 - lambda))
            _set = true
        end
    end
    if !_set
        error("Price must lie inside ribs")
    end
    y
end


init!{C<:AbstractCutOracle, T}(vf::InterpolatedValueFunction{C, T}, m::JuMP.Model, sense, bound, cutmanager) = vf

function passpriceforward!(m::SDDPModel, sp::JuMP.Model)
    stage = ext(sp).stage
    if stage < length(m.stages)
        # pass price forward
        for sp2 in subproblems(m, stage + 1)
            valueoracle(sp2).location = valueoracle(sp).location
        end
    end
end


function solvesubproblem!{C<:AbstractCutOracle, T}(::Type{ForwardPass}, vf::InterpolatedValueFunction{C, T}, m::SDDPModel, sp::JuMP.Model)
    p = vf.location
    # learn noise
    w = sample(vf.noises)
    # update price
    setobjective!(sp, p, w)
    passpriceforward!(m, sp)
end

function backwardpass!{C, T}(m::SDDPModel{InterpolatedValueFunction{C, T}}, settings::Settings)
    for t in (nstages(m)-1):-1:1
        for sp in subproblems(m, t)
            vf = valueoracle(sp)
            ex = ext(sp)

            for (rib, theta) in zip(vf.rib_locations, vf.variables)
                reset!(m.storage)
                for sp2 in subproblems(m, t+1)
                    vf2 = valueoracle(sp2)
                    ex2 = ext(sp2)
                    markov_prob = getstage(m, ex.stage2).transitionprobabilities[ex.markovstate, ex2.markovstate]
                    for noise in zip(vf2.noises)
                        setobjective!(sp2, rib, noise.x)
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
                storecut!(vf.cutmanager, m, sp, cut)
                addcut!(vf, sp, cut)
            end
        end
    end
    # reset!(m.storage)
    # for sp in subproblems(m, 1)
    #     solvesubproblem!(BackwardPass, m, sp)
    # end
end

function addcut!{C<:AbstractCutOracle, T}(vf::InterpolatedValueFunction{C, T}, sp::JuMP.Model, theta::JuMP.Variable, cut::Cut)
    ex = ext(sp)
    affexpr = AffExpr(cut.intercept)
    for i in 1:nstates(sp)
        append!(affexpr, cut.coefficients[i] * ex.states[i].variable)
    end
    _addcut!(ex.sense, sp, theta, affexpr)
end








# function modifyvaluefunction!{C<:AbstractCutOracle}(vf::InterpolatedValueFunction{C}, m::SDDPModel, sp::JuMP.Model)
# end
#
# function rebuildsubproblem!{C<:AbstractCutOracle}(vf::InterpolatedValueFunction{C}, m::SDDPModel, sp::JuMP.Model)
# end
