#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const KW_SYM = (VERSION < v"0.6-")?(:kw):(:(=))
const START  = (VERSION < v"0.6-")?(:start):(esc(:start))

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

getsense(::Type{Max}) = :Max
getsense(::Type{Min}) = :Min
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
        return Min
    elseif s==:Max
        return Max
    else
        error("Unknown optimisation sense $s. Must be :Max or :Min")
    end
end

futureobjective!(::Type{Max}, m::JuMP.Model, bound) = @variable(m, upperbound = bound)
futureobjective!(::Type{Min}, m::JuMP.Model, bound) = @variable(m, lowerbound = bound)

stages(m::SDDPModel) = m.stages
getstage(m::SDDPModel, t) = stages(m)[t]

subproblems(s::Stage) = s.subproblems
getsubproblem(s::Stage, i::Int) = subproblems(s)[i]
subproblems(m::SDDPModel, t) = subproblems(getstage(m, t))
getsubproblem(m::SDDPModel, t, i) = getsubproblem(getstage(m, t), i)

nstages(m::SDDPModel) = length(stages(m))

function n_args(f::Function)
    @assert length(methods(f)) == 1
    return methods(f).mt.max_args-1
end
# No triangular dispatch in v0.5
if VERSION >= v"0.6"
    getel{A, T <: A}(::Type{A}, x::T, t::Int, i::Int=1) = x
    getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int, i::Int=1) = x[t]
    getel{A, T <: A}(::Type{A}, x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
    # Issue #64
    getel(T, f::Function, t::Int, i::Int=1) = f(t, i)::T
else
    # fake triangular dispatch
    singleelement(x, t::Int, i::Int=1) = x
    vectorelement(x, t::Int, i::Int=1) = x[t]
    jaggedelement(x, t::Int, i::Int)   = x[t][i]

    typecompare{T1, T2}(::Type{T1}, t2::T2) = T2 <: T1?(:single):(:false)
    typecompare{T1, T2}(::Type{T1}, t2::Vector{T2}) = T2 <: T1?(:vector):(:false)
    typecompare{T1, T2}(::Type{T1}, t2::Vector{Vector{T2}}) = T2 <: T1?(:jagged):(:false)
    function getel{T1, T2}(typ::Type{T1}, x::T2, t::Int, i::Int=1)
        sym = typecompare(typ, x)
        if sym == :false && isa(x, Function)
            return x(t, i)::T1
        elseif sym == :single
            return singleelement(x, t, i)::T1
        elseif sym == :vector
            return vectorelement(x, t, i)::T1
        elseif sym == :jagged
            return jaggedelement(x, t, i)::T1
        else
            error("Type mismatch: unable to get an object of type $(T1) from $(x).")
        end
    end
end

# Construct a confidence interval
function confidenceinterval(x, conf_level=0.95)
    tstar = quantile(TDist(length(x)-1), 1 - (1 - conf_level)/2)
    err = tstar * std(x)/sqrt(length(x))
    mu = mean(x)
    return mu - err, mu + err
end

function rtol(x, y)
    if abs(y) < 1e-6
        return x - y
    else
        (x - y) / abs(y)
    end
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

_addcut!(::Type{Min}, sp, theta, affexpr) = @constraint(sp, theta >= affexpr)
_addcut!(::Type{Max}, sp, theta, affexpr) = @constraint(sp, theta <= affexpr)


Base.size{T}(x::CachedVector{T}) = (x.n,)
function Base.getindex{T}(x::CachedVector{T}, i)
    if i > x.n
        throw(BoundsError(x, i))
    end
    x.data[i]
end

function Base.setindex!{T}(x::CachedVector{T}, y::T, i)
    if i > x.n
        throw(BoundsError(x, i))
    end
    x.data[i] = y
end
function Base.push!{T}(x::CachedVector{T}, y::T)
    x.n += 1
    if x.n <= length(x.data)
        x.data[x.n] = y
    else
        push!(x.data, y)
    end
end
CachedVector(T) = CachedVector{T}(T[], 0)

function reset!{T}(x::CachedVector{T})
    x.n = 0
end


function samplesubproblem(stage::Stage, last_markov_state::Int)
    newidx = sample(stage.transitionprobabilities[last_markov_state, :])
    return newidx, getsubproblem(stage, newidx)
end

savesolution!(solutionstore::Void, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int) = nothing


function solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model, solutionstore=nothing)
    JuMPsolve(direction, m, sp)
end
hasnoises(sp::JuMP.Model) = length(ext(sp).noises) > 0

function JuMPsolve{T<:IterationDirection}(::Type{T}, ::SDDPModel, sp::JuMP.Model)
    @timeit TIMER "JuMP.solve" begin
        # @assert JuMP.solve(sp) == :Optimal
        @assert jumpsolve(sp) == :Optimal
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
