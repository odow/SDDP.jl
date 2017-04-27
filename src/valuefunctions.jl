#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

stageobjective!(vf::AbstractValueFunction, sp::JuMP.Model, obj) = error("You need this method")
getstageobjective(vf::AbstractValueFunction, sp::JuMP.Model) = error("You need this method")
init!(vf::AbstractValueFunction, m::JuMP.Model, sense, bound, cutmanager) = error("You need this method")
modifyvaluefunction!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model) = error("You need this method")

stageobjective!(sp::JuMP.Model, obj::AffExpr) = stageobjective!(valueoracle(sp), sp, obj)
stageobjective!(sp::JuMP.Model, obj::JuMP.Variable) = stageobjective!(sp, AffExpr(obj))
getstageobjective(sp::JuMP.Model) = getstageobjective(valueoracle(sp), sp)
modifyvaluefunction!(m::SDDPModel, sp::JuMP.Model) = modifyvaluefunction!(valueoracle(sp), m, sp)


struct DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
    cutmanager::C
    stageobjective::QuadExpr
    theta::JuMP.Variable
end
function DefaultValueFunction(m::JuMP.Model, sense, bound, cutmanager=DefaultCutOracle())
    DefaultValueFunction(cutmanager, QuadExpr(0), futureobjective!(sense, m, bound))
end

function init!(::Type{DefaultValueFunction}, m::JuMP.Model, sense, bound, cutmanager)
    DefaultValueFunction(m::JuMP.Model, sense, bound, cutmanager)
end

function stageobjective!{C<:AbstractCutOracle}(vf::DefaultValueFunction{C}, sp::JuMP.Model, obj::AffExpr)
    append!(vf.stageobjective, QuadExpr(obj))
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), obj)
    else
        JuMP.setobjective(sp, getsense(sp), obj + vf.theta)
    end
end
getstageobjective{C<:AbstractCutOracle}(vf::DefaultValueFunction{C}, sp::JuMP.Model) = getvalue(vf.stageobjective)

# optional method that can be overloaded
function solvesubproblem!{C<:AbstractCutOracle}(::Type{BackwardPass}, vf::DefaultValueFunction{C}, m::SDDPModel, sp::JuMP.Model)
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

function constructcut(m::SDDPModel, sp::JuMP.Model)
    storage = m.storage
    # theta <=/>= E[ (y - πᵀx̄) + πᵀx ]
    intercept = 0.0
    coefficients = zeros(nstates(sp))
    for i in 1:length(storage.objective)
        intercept += storage.modifiedprobability[i] * (storage.objective[i] - dot(getstage(m, ext(sp).stage).state, storage.duals[i]))
        for j in 1:nstates(sp)
            coefficients[j] += storage.modifiedprobability[i] * storage.duals[i][j]
        end
    end
    Cut(intercept, coefficients)
end

_addcut!(::Type{Min}, sp, theta, affexpr) = @constraint(sp, theta >= affexpr)
_addcut!(::Type{Max}, sp, theta, affexpr) = @constraint(sp, theta <= affexpr)
# valuefunction(sp::JuMP.Model) = ext(sp).valuefunction
function modifyvaluefunction!{C<:AbstractCutOracle}(vf::DefaultValueFunction{C}, m::SDDPModel, sp::JuMP.Model)
    ex = ext(sp)
    I = 1:length(m.storage.objective)
    for i in I
        m.storage.probability[i] *= getstage(m, ex.stage).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
    modifyprobability!(ex.riskmeasure, m.storage.modifiedprobability.data[I], m.storage.probability.data[I], m.storage.objective.data[I])
    cut = constructcut(m, sp)
    storecut!(vf.cutmanager, m, sp, cut)

    affexpr = AffExpr(cut.intercept)
    for i in 1:nstates(sp)
        affexpr += cut.coefficients[i] * ex.states[i].variable
    end
    c = _addcut!(ex.sense, sp, vf.theta, affexpr)

    for i in I
        m.storage.probability[i] /= getstage(m, ex.stage).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
end
