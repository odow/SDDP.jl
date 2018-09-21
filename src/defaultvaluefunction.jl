#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

setstageobjective!(sp::JuMP.Model, obj...) = setstageobjective!(valueoracle(sp), sp, obj...)
getstageobjective(sp::JuMP.Model) = getstageobjective(valueoracle(sp), sp)


mutable struct DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
    cutmanager::C
    theta::JuMP.Variable
end

"""
    DefaultValueFunction(cutoracle=DefaultCutOracle())

The default value function.
"""
DefaultValueFunction(cutoracle=DefaultCutOracle()) = DefaultValueFunction(cutoracle, JuMP.Variable(JuMP.Model(), 0))

cutoracle(vf::DefaultValueFunction) = vf.cutmanager

"""
    summarise(::Type{V})

Return a short string describing the value function
"""
function summarise end

summarise(::Type{DefaultValueFunction{C}}) where {C} = "Default"

"""
    initialise_value_function(vf::V, sp::JuMP.Model, sense::OptimisationSense, bound::Real)

Initialize a value function
"""
function initializevaluefunction end

function initializevaluefunction(vf::DefaultValueFunction{C}, m::JuMP.Model, sense::OptimisationSense, bound::Real) where C
    vf.theta = futureobjective!(sense, m, bound)
    vf
end


"""
    getstageobjective(v::V, sp::JuMP.Model)

Return a Float64 of the stage objective
"""
function getstageobjective end

function getstageobjective(vf::DefaultValueFunction, sp::JuMP.Model)
    if ext(sp).finalstage
        return JuMP.getobjectivevalue(sp)
    else
        return JuMP.getobjectivevalue(sp) - JuMP.getvalue(vf.theta)
    end
end

"""
    setstageobjective!(v::V, sp::JuMP.Model, obj)

Set the full objective
"""
function setstageobjective! end

function setstageobjective!(vf::DefaultValueFunction, sp::JuMP.Model, obj)
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), obj)
    else
        JuMP.setobjective(sp, getsense(sp), obj + vf.theta)
    end
    JuMP.setobjectivesense(sp, getsense(sp))
end

"""
    modifyvaluefunction!(m::SDDPModel{V}, settings::Settings, sp::JuMP.Model, writecuts::Bool=true)
"""
function modifyvaluefunction! end

function modifyvaluefunction!(m::SDDPModel{V}, settings::Settings, sp::JuMP.Model) where V<:DefaultValueFunction
    ex = ext(sp)
    vf = valueoracle(sp)
    I = 1:length(m.storage.objective)
    # TODO: improve this to reduce the extra memory usage
    current_transition = copy(m.storage.probability.data[I])
    for i in I
        m.storage.probability[i] *= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
    @timeit TIMER "risk measure" begin
        modifyprobability!(ex.riskmeasure,
            view(m.storage.modifiedprobability.data, I),
            m.storage.probability.data[I],
            m.storage.objective.data[I],
            m,
            sp
        )
    end
    cut = constructcut(m, sp)

    if !settings.is_asynchronous && isopen(settings.cut_output_file)
        # only write the cut to file if it is open and
        # we are in serial mode. Async cut writing happens on
        # master thread
        @timeit TIMER "cut to file" begin
            writecut!(settings.cut_output_file, cut, ex.stage, ex.markovstate)
        end
    end

    # add the cut to value function and JuMP model
    addcut!(m, sp, cut)

    # if we are solving asynchronously, need to save it
    # for passing
    if settings.is_asynchronous
        storeasynccut!(m, sp, cut)
    end

    for i in I
        m.storage.probability[i] = current_transition[i]
    end
end

function addcut!(m::SDDPModel{V}, sp::JuMP.Model, cut::Cut) where V<:DefaultValueFunction
    vf = valueoracle(sp)
    # store cut in oracle
    store_cut(cutoracle(vf), m, sp, cut)
    addcuttoJuMPModel!(sp, cut)
end

function addcuttoJuMPModel!(sp::JuMP.Model, cut::Cut)
    # add constraint to JuMP model
    ex = ext(sp)
    vf = valueoracle(sp)
    affexpr = cuttoaffexpr(sp, cut)
    addcutconstraint!(ex.sense, sp, vf.theta, affexpr)
end

function solvesubproblem!(::Type{BackwardPass}, m::SDDPModel, sp::JuMP.Model, incoming_probability::Float64=1.0)
    ex = ext(sp)
    if hasnoises(sp)
        for i in 1:length(ex.noiseprobability)
            setnoise!(sp, ex.noises[i])
            JuMPsolve(BackwardPass, m, sp)
            push!(m.storage.objective, getobjectivevalue(sp))
            push!(m.storage.noise, i)
            push!(m.storage.probability, incoming_probability * ex.noiseprobability[i])
            push!(m.storage.modifiedprobability, incoming_probability * ex.noiseprobability[i])
            push!(m.storage.markov, ex.markovstate)
            push!(m.storage.duals, zeros(nstates(sp)))
            saveduals!(m.storage.duals[end], sp)
        end
    else
        JuMPsolve(BackwardPass, m, sp)
        push!(m.storage.objective, getobjectivevalue(sp))
        push!(m.storage.noise, 0)
        push!(m.storage.probability, incoming_probability)
        push!(m.storage.modifiedprobability, incoming_probability)
        push!(m.storage.markov, ex.markovstate)
        push!(m.storage.duals, zeros(nstates(sp)))
        saveduals!(m.storage.duals[end], sp)
    end
end

"""
    rebuildsubproblem!(m::SDDPModel, sp::JuMP.Model)

"""
function rebuildsubproblem! end

function rebuildsubproblem!(m::SDDPModel{DefaultValueFunction{C}}, sp::JuMP.Model) where C<:AbstractCutOracle
    vf = valueoracle(sp)
    n = n_args(m.build!)
    ex = ext(sp)
    for i in 1:nstates(sp)
        pop!(ex.states)
    end
    for i in 1:length(ex.noises)
        pop!(ex.noises)
    end
    sp2 = Model(solver = sp.solver)

    vf.theta = futureobjective!(ex.sense, sp2, ex.problembound)

    sp2.ext[:SDDP] = ex
    if n == 2
        m.build!(sp2, ex.stage)
    elseif n == 3
        m.build!(sp2, ex.stage, ex.markovstate)
    end
    for cut in valid_cuts(vf.cutmanager)
        addcuttoJuMPModel!(sp2, cut)
    end
    m.stages[ex.stage].subproblems[ex.markovstate] = sp2
end

"""
    loadcuts!(m::SDDPModel, filename::String)

Load cuts from the file created using the `cut_output_file` argument
in `solve`.

### Example

    m = SDDPModel() do ... end
    status = solve(m; cut_output_file="path/to/m.cuts")`
    m2 = SDDPModel() do ... end
    loadcuts!(m2, "path/to/m.cuts")

"""
function loadcuts! end

function loadcuts!(m::SDDPModel{DefaultValueFunction{C}}, filename::String) where C
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
            addcut!(m, sp, cut)
        end
    end
end

#=
    To enable asynchronous solutions, it is necessary to implement
        asynccutstoragetype(::Type{V})::Type{T}
        asynccutstorage(m::SDDPModel{V}, sp::JuMP.Model, cut::Cut)::T
        addasynccut!(m::SDDPModel{V}, cut::T)
=#
"""
    asynccutstoragetype(::Type{V})

Used to pass cuts between processors
"""
function asynccutstoragetype end

asynccutstoragetype(::Type{DefaultValueFunction{C}}) where {C} = Tuple{Int, Int, Cut}

"""
    addasynccut!(m::SDDPModel{V}, cut::T)

Where T is asyncutstoragetype(V)
"""
function addasynccut end

function addasynccut!(m::SDDPModel{DefaultValueFunction{C}}, cut::Tuple{Int, Int, Cut}) where C
    sp = getsubproblem(m, cut[1], cut[2])
    addcut!(m, sp, cut[3])
end
