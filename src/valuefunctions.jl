#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

stageobjective!(sp::JuMP.Model, obj...) = stageobjective!(valueoracle(sp), sp, obj...)
getstageobjective(sp::JuMP.Model) = getstageobjective(valueoracle(sp), sp)

type DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
    cutmanager::C
    stageobjective::QuadExpr
    theta::JuMP.Variable
end
DefaultValueFunction(cutoracle=DefaultCutOracle()) = DefaultValueFunction(cutoracle, QuadExpr(0.0), JuMP.Variable(JuMP.Model(), 0))

summarise{C}(::Type{DefaultValueFunction{C}}) = "Default"

function init!{C}(vf::DefaultValueFunction{C}, m::JuMP.Model, sense, bound)
    vf.theta = futureobjective!(sense, m, bound)
    vf
end

function stageobjective!(vf::DefaultValueFunction, sp::JuMP.Model, obj)
    append!(vf.stageobjective, QuadExpr(obj))
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), obj)
    else
        JuMP.setobjective(sp, getsense(sp), obj + vf.theta)
    end
    JuMP.setobjectivesense(sp, getsense(sp))
end

getstageobjective(vf::DefaultValueFunction, sp::JuMP.Model) = getvalue(vf.stageobjective)

function writecut!(filename::String, cut::Cut, stage::Int, markovstate::Int)
    open(filename, "a") do file
        write(file, "$(stage), $(markovstate), $(cut.intercept)")
        for pi in cut.coefficients
            write(file, ",$(pi)")
        end
        write(file, "\n")
    end
end

function loadcuts!{C}(m::SDDPModel{DefaultValueFunction{C}}, filename::String)
    open(filename, "r") do file
        while true
            line      = readline(file)
            line == nothing || line == "" && break
            items     = split(line, ",")
            stage     = parse(Int, items[1])
            ms        = parse(Int, items[2])
            intercept = parse(Float64, items[3])
            coefficients = [parse(Float64, i) for i in items[4:end]]
            cut = Cut(intercept, coefficients)
            sp = getsubproblem(m, stage, ms)
            vf = valueoracle(sp)
            # Add cut to the cut manager
            storecut!(vf.cutmanager, m, sp, cut)
            # Add cut as a constraint to the subproblem
            addcut!(vf, sp, cut)
        end
    end
end

function modifyvaluefunction!{V<:DefaultValueFunction}(m::SDDPModel{V}, settings::Settings, sp::JuMP.Model)
    ex = ext(sp)
    vf = valueoracle(sp)
    I = 1:length(m.storage.objective)
    for i in I
        m.storage.probability[i] *= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
    modifyprobability!(ex.riskmeasure,
        view(m.storage.modifiedprobability.data, I),
        m.storage.probability.data[I],
        m.storage.objective.data[I],
        m,
        sp
    )
    cut = constructcut(m, sp)

    if settings.cut_output_file != ""
        writecut!(settings.cut_output_file, cut, ex.stage, ex.markovstate)
    end

    storecut!(vf.cutmanager, m, sp, cut)
    addcut!(vf, sp, cut)
    storecut!(m, sp, cut)
    for i in I
        m.storage.probability[i] /= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
end
storecut!(m::SDDPModel, sp::JuMP.Model, cut::Cut, args...) = nothing

function addcut!(vf::DefaultValueFunction, sp::JuMP.Model, cut::Cut)
    ex = ext(sp)
    affexpr = AffExpr(cut.intercept)
    for i in 1:nstates(sp)
        append!(affexpr, cut.coefficients[i] * ex.states[i].variable)
    end
    _addcut!(ex.sense, sp, vf.theta, affexpr)
end

function solvesubproblem!{V<:DefaultValueFunction}(::Type{BackwardPass}, m::SDDPModel{V}, sp::JuMP.Model)
    ex = ext(sp)
    if hasnoises(sp)
        for i in 1:length(ex.noiseprobability)
            setnoise!(sp, ex.noises[i])
            JuMPsolve(BackwardPass, m, sp)
            push!(m.storage.objective, getobjectivevalue(sp))
            push!(m.storage.noise, i)
            push!(m.storage.probability, ex.noiseprobability[i])
            push!(m.storage.modifiedprobability, ex.noiseprobability[i])
            push!(m.storage.markov, ex.markovstate)
            push!(m.storage.duals, zeros(nstates(sp)))
            saveduals!(m.storage.duals[end], sp)
        end
    else
        JuMPsolve(BackwardPass, m, sp)
        push!(m.storage.objective, getobjectivevalue(sp))
        push!(m.storage.noise, 0)
        push!(m.storage.probability, 1.0)
        push!(m.storage.modifiedprobability, 1.0)
        push!(m.storage.markov, ex.markovstate)
        push!(m.storage.duals, zeros(nstates(sp)))
        saveduals!(m.storage.duals[end], sp)
    end
end

function rebuildsubproblem!{C<:AbstractCutOracle}(m::SDDPModel{DefaultValueFunction{C}}, sp::JuMP.Model)
    vf = valueoracle(sp)
    n = n_args(m.build!)
    ex = ext(sp)
    for i in 1:nstates(sp)
        pop!(ex.states)
    end
    for i in 1:length(ex.noises)
        pop!(ex.noises)
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

rebuildsubproblem!(m::SDDPModel{DefaultValueFunction{DefaultCutOracle}}, sp::JuMP.Model) = nothing
