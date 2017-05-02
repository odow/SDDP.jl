#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

mutable struct DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
    cutmanager::C
    stageobjective::QuadExpr
    theta::JuMP.Variable
end
DefaultValueFunction() = DefaultValueFunction

init!(::Type{DefaultValueFunction}, m::JuMP.Model, sense, bound, cutmanager) = DefaultValueFunction(
    cutmanager,
    QuadExpr(0.0),
    futureobjective!(sense, m, bound)
)

function stageobjective!(vf::DefaultValueFunction, sp::JuMP.Model, obj)
    append!(vf.stageobjective, QuadExpr(obj))
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), obj)
    else
        JuMP.setobjective(sp, getsense(sp), obj + vf.theta)
    end
end

getstageobjective(vf::DefaultValueFunction, sp::JuMP.Model) = getvalue(vf.stageobjective)

function modifyvaluefunction!(vf::DefaultValueFunction, m::SDDPModel, sp::JuMP.Model)
    ex = ext(sp)
    I = 1:length(m.storage.objective)
    for i in I
        m.storage.probability[i] *= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
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
    for i in I
        m.storage.probability[i] /= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
end

function addcut!(vf::DefaultValueFunction, sp, cut::Cut)
    ex = ext(sp)
    affexpr = AffExpr(cut.intercept)
    for i in 1:nstates(sp)
        append!(affexpr, cut.coefficients[i] * ex.states[i].variable)
    end
    _addcut!(ex.sense, sp, vf.theta, affexpr)
end

function solvesubproblem!(::Type{BackwardPass}, vf::DefaultValueFunction, m::SDDPModel, sp::JuMP.Model)
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

function rebuildsubproblem!{C<:AbstractCutOracle}(vf::DefaultValueFunction{C}, m::SDDPModel, sp::JuMP.Model)
    n = n_args(m.build!)
    ex = ext(sp)
    for i in 1:nstates(sp)
        pop!(ex.states)
    end
    for i in 1:length(ex.scenarios)
        pop!(ex.scenarios)
    end
    sp = Model(solver = m.lpsolver)

    vf.stageobjective = QuadExpr(0.0)
    vf.theta = futureobjective!(ex.sense, sp, ex.problembound)

    sp.ext[:SDDP] = ex
    if n == 2
        m.build!(sp, ex.stage)
    elseif n == 3
        m.build!(sp, ex.stage, ex.markovstate)
    end
    for cut in validcuts(vf.cutmanager)
        addcut!(vf, sp, cut)
    end
    m.stages[ex.stage].subproblems[ex.markovstate] = sp
end
rebuildsubproblem!(vf::DefaultValueFunction{DefaultCutOracle}, m::SDDPModel, sp::JuMP.Model) = nothing
