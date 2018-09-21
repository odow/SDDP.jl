#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const KW_SYM = (VERSION < v"0.6-") ? (:kw) : (:(=))
const START  = (VERSION < v"0.6-") ? (:start) : (esc(:start))

ext(m::JuMP.Model) = m.ext[:SDDP]::SubproblemExt
isext(m::JuMP.Model) = isa(m.ext[:SDDP], SubproblemExt)
valueoracle(sp::JuMP.Model) = ext(sp).valueoracle
cutoracle(sp::JuMP.Model) = cutoracle(valueoracle(sp))

"""
    getbound(m)

# Description

Get the lower (if minimizing), or upper (if maximizing) bound of the solved SDDP
model `m`.
"""
function getbound(m::SDDPModel)
    if length(m.log) > 0
        return m.log[end].bound
    else
        error("Model not solved.")
    end
end

getsense(::Max) = :Max
getsense(::Min) = :Min
getsense(m::JuMP.Model) = getsense(ext(m).sense)
function worstcase(s)
    if s==:Min
        return Inf
    elseif s == :Max
        return -Inf
    end
end
function optimisationsense(s::Symbol)
    if s==:Min
        return Min()
    elseif s==:Max
        return Max()
    else
        error("Unknown optimisation sense $s. Must be :Max or :Min")
    end
end

futureobjective!(::Max, m::JuMP.Model, bound) = @variable(m, upperbound = bound)
futureobjective!(::Min, m::JuMP.Model, bound) = @variable(m, lowerbound = bound)

stages(m::SDDPModel) = m.stages
getstage(m::SDDPModel, t) = stages(m)[t]

subproblems(s::Stage) = s.subproblems
getsubproblem(s::Stage, i::Int) = subproblems(s)[i]
subproblems(m::SDDPModel, t) = subproblems(getstage(m, t))

"""
    getsubproblem(m::SDDPModel, t::Int, i::Int)

Get the subproblem in stage `t` and Markov state `i` from the SDDPModel `m`.
"""
getsubproblem(m::SDDPModel, t, i) = getsubproblem(getstage(m, t), i)

nstages(m::SDDPModel) = length(stages(m))

function n_args(f::Function)
    @assert length(methods(f)) == 1
    return methods(f).mt.max_args-1
end

getel(::Type{A}, x::T, t::Int, i::Int=1) where {A, T <: A} = x
getel(::Type{A}, x::Vector{T}, t::Int, i::Int=1) where {A, T <: A} = x[t]
getel(::Type{A}, x::Vector{Vector{T}}, t::Int, i::Int) where {A, T <: A} = x[t][i]
# Issue #64
getel(T, f::Function, t::Int, i::Int=1) = f(t, i)::T

# Construct a confidence interval
function confidenceinterval(x, conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    err = tstar * std(x)/sqrt(length(x))
    mu = mean(x)
    return mu - err, mu + err
end

function sample(x::AbstractVector{Float64})
    r = rand()
    for i in 1:length(x)
        @inbounds r -= x[i]
        if r < eps(Float64)
            return i
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end

contains(x, l, u) = x >= l && x <= u

applicable(iteration, frequency) = frequency > 0 && mod(iteration, frequency) == 0

function constructcut(m::SDDPModel, sp::JuMP.Model)
    # theta <=/>= E[ (y - πᵀx̄) + πᵀx ]
    intercept = 0.0
    coefficients = zeros(nstates(sp))
    for i in 1:length(m.storage.objective)
        intercept += m.storage.modifiedprobability[i] * (m.storage.objective[i] - dot(getstage(m, ext(sp).stage).state, m.storage.duals[i]))
        for j in 1:nstates(sp)
            coefficients[j] += m.storage.modifiedprobability[i] * m.storage.duals[i][j]
        end
    end
    Cut(intercept, coefficients)
end

addcutconstraint!(sense::Min, sp, theta, affexpr) = @constraint(sp, theta >= affexpr)
addcutconstraint!(sense::Max, sp, theta, affexpr) = @constraint(sp, theta <= affexpr)


Base.size(x::CachedVector{T}) where {T} = (x.n,)
function Base.getindex(x::CachedVector{T}, i) where T
    if i > x.n
        throw(BoundsError(x, i))
    end
    x.data[i]
end

function Base.setindex!(x::CachedVector{T}, y::T, i) where T
    if i > x.n
        throw(BoundsError(x, i))
    end
    x.data[i] = y
end
function Base.push!(x::CachedVector{T}, y::T) where T
    x.n += 1
    if x.n <= length(x.data)
        x.data[x.n] = y
    else
        push!(x.data, y)
    end
end
CachedVector(T) = CachedVector{T}(T[], 0)

function reset!(x::CachedVector{T}) where T
    x.n = 0
end


function samplesubproblem(stage::Stage, last_markov_state::Int)
    newidx = sample(stage.transitionprobabilities[last_markov_state, :])
    return newidx, getsubproblem(stage, newidx)
end
samplesubproblem(stage::Stage, last_markov_state::Int, solutionstore::Void) = samplesubproblem(stage, last_markov_state)

function samplesubproblem(stage::Stage, last_markov_state, solutionstore::Dict{Symbol, Any})
    if length(solutionstore[:noise]) >= stage.t
        idx = solutionstore[:markov][stage.t]
        return idx, getsubproblem(stage, idx)
    else
        return samplesubproblem(stage, last_markov_state)
    end
end

savesolution!(solutionstore::Void, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int) = nothing


function solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model, solutionstore=nothing)
    JuMPsolve(direction, m, sp)
end
hasnoises(sp::JuMP.Model) = length(ext(sp).noises) > 0

function presolve!(direction, m, sp)
    nothing
end
function postsolve!(direction, m, sp)
    nothing
end

function JuMPsolve(direction::Type{T}, m::SDDPModel, sp::JuMP.Model) where T<:IterationDirection
    @timeit TIMER "JuMP.solve" begin
        # @assert JuMP.solve(sp) == :Optimal
        presolve!(direction, m, sp)
        status = jumpsolve(sp)
        if status != :Optimal
            filepath = joinpath(pwd(), "infeasible_subproblem.lp")
            JuMP.writeLP(sp, filepath; genericnames=false)
            error("""Model in stage $(ext(sp).stage) and markov state $(ext(sp).markovstate)
            was not solved to Optimality. I wrote the offending MPS file to
            $(filepath).

            This is most commonly caused by numerical issues with the solver.
            Consider reformulating the model or try different solver parameters.
            """)
        end
        postsolve!(direction, m, sp)
    end
end



"""
    SDDP.savemodel!(filename::String, m::SDDPModel)

Save the SDDPModel `m` to the location `filename`. Can be loaded at a later date
with `m = SDDP.loadmodel(filename)`.

Note: this function relies in the internal Julia `Base.serialize`function. It
should not be relied on to save an load models between versions of Julia (i.e
between v0.5 and v0.6). For a longer-term solution, see `SDDP.loadcuts!` for help.
"""
function savemodel!(filename::String, m::SDDPModel)
    for stage in stages(m)
        for sp in subproblems(stage)
            sp.internalModelLoaded = false
        end
    end
    m.ext[:serializer_version] = Base.Serializer.ser_version
    save!(filename, m)
end

function save!(filename::String, x)
    open(filename, "w") do io
        serialize(io, x)
    end
end
function load(filename::String)
    io = open(filename, "r")
    x = deserialize(io)
    close(io)
    x
end

"""
    loadmodel(filename::String)

Load a model from the location `filename` that was saved using `SDDP.savemodel!`.

Note: this function relies in the internal Julia `Base.serialize`function. It
should not be relied on to save an load models between versions of Julia (i.e
between v0.5 and v0.6). For a longer-term solution, see `SDDP.loadcuts!` for help.
"""
function loadmodel(filename::String)
    # We're going to need this to be merged into the copy of Julia that the user
    # is running
    # https://github.com/JuliaLang/julia/pull/21799
    m = load(filename)
    if Base.Serializer.ser_version != m.ext[:serializer_version]
        error("The Base serializer has changed. You should use the JLD package
        or similar to save models long-term.")
    end
    #==
        Deepcopying (or serializing) an ObjectIdDict causes the hash table to
        get out of whack. We should just rehash it once we've finished.
    ==#
    for stage in stages(m)
        for sp in subproblems(stage)
            Base.rehash!(sp.varData)
        end
    end
    m
end

function dominates(sense, trial, incumbent)
    if sense == :Min
        return trial > incumbent
    elseif sense == :Max
        return trial < incumbent
    end
    error("Sense $sense not recognised")
end

function cuttoaffexpr(sp::Model, cut::Cut)
    x = AffExpr(cut.intercept)
    for (idx, coef) in enumerate(cut.coefficients)
        append!(x, coef * states(sp)[idx].variable)
    end
    return x
end

writecut!(io, cut::Cut, stage::Int, markovstate::Int) = writecut!(io, stage, markovstate, cut)

function writecut!(io::IO, stage::Int, markovstate::Int, cut::Cut)
    write(io, string(stage), ",", string(markovstate), ",", string(cut.intercept))
    for pi in cut.coefficients
        write(io, ",", string(pi))
    end
    write(io, "\n")
end
writecut!(io, cut::Tuple) = writecut!(io, cut...)

"""
    writecuts!(filename::String, m::SDDPModel; onlyvalid=false)

Writes all cuts from model m to `filename`.

If `onlyvalid` is true, write the cuts returned from `valid_cuts`, else write the
cuts returned from `all_cuts`.
"""
function writecuts!(filename::String, m::SDDPModel; onlyvalid=false)
    open(filename, "w") do io
        writecuts!(io, m, onlyvalid=onlyvalid)
    end
end
function writecuts!(io::IO, m::SDDPModel; onlyvalid=false)
    for (t, stage) in enumerate(stages(m))
        for (i, sp) in enumerate(subproblems(stage))
            cut_or = cutoracle(sp)
            cuts = onlyvalid ? valid_cuts(cut_or) : all_cuts(cut_or)
            for cut in cuts
                writecut!(io, t, i, cut)
            end
        end
    end
end
