#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

stageobjective!(vf::AbstractValueFunction, sp::JuMP.Model, obj...) = error("You need this method")
getstageobjective(vf::AbstractValueFunction, sp::JuMP.Model) = error("You need this method")
init!(vf::AbstractValueFunction, m::JuMP.Model, sense, bound) = error("You need this method")
modifyvaluefunction!(vf::AbstractValueFunction, m::SDDPModel, settings::Settings, sp::JuMP.Model) = error("You need this method")
rebuildsubproblem!(vf::AbstractValueFunction, m::SDDPModel, sp::JuMP.Model) = nothing
summarise{T<:AbstractValueFunction}(::Type{T}) = "$T"

stageobjective!(sp::JuMP.Model, obj...) = stageobjective!(valueoracle(sp), sp, obj...)
getstageobjective(sp::JuMP.Model) = getstageobjective(valueoracle(sp), sp)
modifyvaluefunction!(m::SDDPModel, settings::Settings, sp::JuMP.Model) = modifyvaluefunction!(valueoracle(sp), m, settings, sp)
rebuildsubproblem!(m::SDDPModel, sp::JuMP.Model) = rebuildsubproblem!(valueoracle(sp), m, sp)

# mutable struct DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
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

function readcuts!{C}(m::SDDPModel{DefaultValueFunction{C}}, filename::String)
    open(filename, "r") do file
        while true
            line      = readline(f)
            items     = split(line, ",")
            stage     = parse(Int, items[1])
            ms        = parse(Int, items[2])
            intercept = parse(Float64, items[3])
            coefficients = [parse(Float64, i) for i in items[4:end]]
            cut = Cut(intercept, coefficients)
            sp = getsubproblem(m, stage, ms)
            storecut!(valueoracle(sp).cutmanager, m, sp, cut)
        end
    end
    for (t, stage) in enumerate(stages(m))
        t == length(stages(m)) && continue
        for sp in subproblems(stage)
            rebuildsubproblem!(m, sp)
        end
    end
end

function modifyvaluefunction!(vf::DefaultValueFunction, m::SDDPModel, settings::Settings, sp::JuMP.Model)
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
    if hasnoises(sp)
        for i in 1:length(ex.noiseprobability)
            setnoise!(sp, ex.noises[i])
            if sp.solvehook == nothing
                @assert JuMP.solve(sp) == :Optimal
            else
                @assert JuMP.solve(sp, require_duals = true) == :Optimal
            end
            push!(m.storage.objective, getobjectivevalue(sp))
            push!(m.storage.noise, i)
            push!(m.storage.probability, ex.noiseprobability[i])
            push!(m.storage.modifiedprobability, ex.noiseprobability[i])
            push!(m.storage.markov, ex.markovstate)
            push!(m.storage.duals, zeros(nstates(sp)))
            saveduals!(m.storage.duals[end], sp)
        end
    else
        if sp.solvehook == nothing
            @assert JuMP.solve(sp) == :Optimal
        else
            @assert JuMP.solve(sp, require_duals = true) == :Optimal
        end
        push!(m.storage.objective, getobjectivevalue(sp))
        push!(m.storage.noise, 0)
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
rebuildsubproblem!(vf::DefaultValueFunction{DefaultCutOracle}, m::SDDPModel, sp::JuMP.Model) = nothing
