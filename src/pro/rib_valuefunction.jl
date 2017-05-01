#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

struct Noise{T}
    x::T
    probability::Float64
end
function sample{T}(x::Vector{Noise{T}})
    r = rand()
    for i in 1:length(x)
        @inbounds r -= x[i].probablity
        if r < eps(Float64)
            return i
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end

mutable struct InterpolatedValueFunction{C} <: AbstractValueFunction
    location::Float64
    rib_locations::Vector{Float64}
    variables::Vector{JuMP.Variable}
    cutoracles::Vector{C}
    noises::Vector{Noise}
    objective::Function
    dynamics::Function
end

function interpolate(vf::InterpolatedValueFunction)
    y = AffExpr(0.0)
    _set = false
    for i in 2:length(vf.rib_locations)
        if vf.location <= vf.rib_locations[i]
            lambda = (vf.location - vf.rib_locations[i-1]) / (vf.rib_locations[i] - vf.rib_locations[i-1])
            y += vf.variables[i-1] * lambda
            y += vf.variables[i] * (1 - lambda)
            _set = true
        end
    end
    if !_set
        error("Price must lie inside ribs")
    end
    y
end

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

function stageobjective!(vf::InterpolatedValueFunction, sp::JuMP.Model, obj...)
    error("You've chosen the InterpolatedValueFunction. You should use priceprocess!() instead of stageobjective!.")
end

function getstageobjective(vf::InterpolatedValueFunction, sp::JuMP.Model)
    getvalue(vf.objective(vf.ribs.location))
end

function init!(vf::InterpolatedValueFunction, m::JuMP.Model, sense, bound, cutmanager)
end

function passpriceforward!(m::SDDPModel, sp::JuMP.Model)
    stage = ext(sp).stage
    if stage < length(m.stages)
        # pass price forward
        for sp2 in subproblems(m, stage + 1)
            valueoracle(sp2).location = vf.location
        end
    end
end


function solvesubproblem!(::Type{ForwardPass}, vf::InterpolatedValueFunction, m::SDDPModel, sp::JuMP.Model)
    p = vf.location
    # learn noise
    w = sample(vf.noises)
    # update price
    setobjective!(sp, p, w)
    passpriceforward!(m, sp)
end


function solvesubproblem!(::Type{BackwardPass}, vf::InterpolatedValueFunction, m::SDDPModel, sp::JuMP.Model)
end
# function modifyvaluefunction!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model)
# end
# function rebuildsubproblem!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model)
# end
