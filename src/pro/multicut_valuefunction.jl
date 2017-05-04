#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export MultiCutValueFunction

struct MultiCutValueFunction{C<:AbstractCutManager} <: AbstractValueFunction
    n::Int
    stageobjective::QuadExpr
    variables::Vector{JuMP.Variable}
    cutmanager::Vector{C}
end
MultiCutValueFunction(nscenarios::Int, cutoracle=DefaultCutOracle()) = MultiCutValueFunction(
    nscenarios, QuadExpr(0.0), JuMP.Variable[], typeof(cutoracle)[]
)

function stageobjective!(vf::MultiCutValueFunction, sp::JuMP.Model, obj)
    append!(vf.stageobjective, QuadExpr(obj))
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), obj)
    else
        JuMP.setobjective(sp, getsense(sp), obj + vf.theta)
    end
end

getstageobjective(vf::MultiCutValueFunction, sp::JuMP.Model) = getvalue(vf.stageobjective)

# init!(::Type{MultiCutValueFunction}, m::JuMP.Model, sense, bound, cutmanager)

# modifyvaluefunction!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model)

# rebuildsubproblem!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model) = nothing

function solvesubproblem!(::Type{BackwardPass}, vf::MultiCutValueFunction, m::SDDPModel, sp::JuMP.Model)
    ex = ext(sp)
    if hasscenarios(sp)
        for i in 1:length(ex.scenarioprobability)
            setscenario!(sp, ex.scenarios[i])
            @assert JuMP.solve(sp) == :Optimal
            push!(m.storage.objective, getobjectivevalue(sp))
            push!(m.storage.scenario, i)
            push!(m.storage.probability, ex.scenarioprobability[i])
            push!(m.storage.modifiedprobability, ex.scenarioprobability[i])
            push!(m.storage.markov, ex.markovstate)
            push!(m.storage.duals, zeros(nstates(sp)))
            saveduals!(m.storage.duals[end], sp)
        end
    else
        @assert JuMP.solve(sp) == :Optimal
        push!(m.storage.objective, getobjectivevalue(sp))
        push!(m.storage.scenario, 0)
        push!(m.storage.probability, 1.0)
        push!(m.storage.modifiedprobability, 1.0)
        push!(m.storage.markov, ex.markovstate)
        push!(m.storage.duals, zeros(nstates(sp)))
        saveduals!(m.storage.duals[end], sp)
    end
end
