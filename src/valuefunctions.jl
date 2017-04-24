#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

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

function stageobjective!{C<:AbstractCutOracle}(::Type{DefaultValueFunction{C}}, sp::JuMP.Model, obj::AffExpr)
    append!(ext(sp).valueoracle.stageobjective, QuadExpr(obj))
    JuMP.setobjective(sp, getsense(sp), obj + ext(sp).valueoracle.theta)
end
getstageobjective{C<:AbstractCutOracle}(::Type{DefaultValueFunction{C}}, sp::JuMP.Model) = getvalue(ext(sp).valueoracle.stageobjective)

function padstore!(s::BackwardPassStorage, sp::JuMP.Model)
    while s.idx > s.N
        push!(s.objective, 0.0)
        push!(s.scenario, 0)
        push!(s.markovstate, 0)
        push!(s.probability, 0.0)
        push!(s.duals, zeros(nstates(sp)))
        s.N += 1
    end
end

# optional method that can be overloaded
function solvesubproblem!{C<:AbstractCutOracle}(::Type{BackwardPass}, ::Type{DefaultValueFunction{C}}, m::SDDPModel, sp::JuMP.Model)
    ex = ext(sp)
    if hasscenarios(sp)
        for i in 1:length(ex.scenarioprobability)
            setscenario!(sp, ex.scenarios[i])
            status = JuMP.solve(sp)
            m.storage.idx += 1
            padstore!(m.storage, sp)
            m.storage.objective[m.storage.idx] = getobjectivevalue(sp)
            m.storage.scenario[m.storage.idx] = i
            m.storage.probability[m.storage.idx] = ex.scenarioprobability[i]
            m.storage.markovstate[m.storage.idx] = ex.markovstate
            saveduals!(m.storage.duals[m.storage.idx], sp)
        end
    else
        status = JuMP.solve(sp)
        m.storage.idx += 1
        padstore!(m.storage, sp)
        m.storage.objective[m.storage.idx] = getobjectivevalue(sp)
        m.storage.scenario[m.storage.idx] = 0
        m.storage.probability[m.storage.idx] = 1.0
        m.storage.markovstate[m.storage.idx] = ex.markovstate
        saveduals!(m.storage.duals[m.storage.idx], sp)
    end
end

function constructcut(storage)
    # theta <=/>= E[ (y - πᵀx̄) + πᵀx ]
    intercept = 0.0
    coefficients = zeros(nstates(sp))
    for i in 1:storage.n
        intercept += storage.newprobability[i] * (storage.objective[i] - dot(storage.stage, storage.duals[i]))
        for j in 1:nstates(sp)
            coefficients[j] += storage.newprobability[i] * storage.duals[i][j]
        end
    end
    Cut(intercept, coefficients)
end

_addcut!(::Type{Min}, sp, theta, affexpr) = @constraint(sp, theta >= affexpr)
_addcut!(::Type{Max}, sp, theta, affexpr) = @constraint(sp, theta <= affexpr)
# valuefunction(sp::JuMP.Model) = ext(sp).valuefunction
function modifyvaluefunction!{C<:AbstractCutOracle}(::Type{DefaultValueFunction{C}}, m::SDDPModel, sp::JuMP.Model)
    ex = ext(sp)
    I = 1:m.storage.idx
    for i in I
        m.storage.probability[i] *= stage(m, ex.stage).transitionprobabilities[ex.stage, m.storage.markovstate[i]]
    end
    modifyprobability!(ex.riskmeasure, view(m.storage.newprobability, I), view(m.storage.probability, I), view(m.storage.objective, I))
    cut = constructcut(m.storage)

    storecut!(ex.valuefunction.cutmanager, m, ex.stage, ex.markovstate, cut)

    ex = ext(sp)
    affexpr = cut.intercept
    for i in 1:nstates(sp)
        affexpr += cut.coefficients[i] * ex.states[i].variable
    end
    _addcut!(ex.sense, sp, ex.valuefunction.cutmanager, affexpr)

    for i in 1:m.storage.idx
        m.storage.probability[i] /= stage.transitionprobabilities[ex.stage, m.storage.markovstate[i]]
    end
end

stageobjective!{T<:AbstractValueFunction}(::Type{T}, sp::JuMP.Model, obj) = error("You need this method")
getstageobjective{T<:AbstractValueFunction}(::Type{T}, sp::JuMP.Model) = error("You need this method")
init!{T<:AbstractValueFunction}(::Type{T}, m::JuMP.Model, sense, bound, cutmanager) = error("You need this method")
modifyvaluefunction!{T<:AbstractValueFunction}(::Type{T}, m::SDDPModel, sp::JuMP.Model) = error("You need this method")

stageobjective!(sp::JuMP.Model, obj::AffExpr) = stageobjective!(vftype(sp), sp, obj)
stageobjective!(sp::JuMP.Model, obj::JuMP.Variable) = stageobjective!(sp, AffExpr(obj))
getstageobjective(sp::JuMP.Model) = getstageobjective(ext(sp), sp)
modifyvaluefunction!(m::SDDPModel, sp::JuMP.Model) = modifyvaluefunction!(vftype(sp), m, sp)
