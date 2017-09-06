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

function getbound(m::SDDPModel)
    if length(m.log) > 0
        return m.log[end].bound
    else
        error("Model not solved.")
    end
end

getsense(::Type{Max}) = :Max
getsense(::Type{Min}) = :Min
getsense(s::Symbol) = (s==:Min)?Min:Max

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
nsubproblems(m::SDDPModel, t::Int) = length(subproblems(m, t))

function n_args(f::Function)
    @assert length(methods(f)) == 1
    return methods(f).mt.max_args-1
end
# No triangular dispatch in v0.5
# getel{A, T <: A}(::Type{A}, x::T, t::Int, i::Int) = x
# getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int, i::Int) = x[t]
# getel{A, T <: A}(::Type{A}, x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
# getel{A, T <: A}(::Type{A}, x::T, t::Int) = x
# getel{A, T <: A}(::Type{A}, x::Vector{T}, t::Int) = x[t]
getel{T}(::Type{T}, x::T, t::Int, i::Int) = x
getel{T}(::Type{T}, x::Vector{T}, t::Int, i::Int) = x[t]
getel{T}(::Type{T}, x::Vector{Vector{T}}, t::Int, i::Int) = x[t][i]
getel{T}(::Type{T}, x::T, t::Int) = x
getel{T}(::Type{T}, x::Vector{T}, t::Int) = x[t]

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


function solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model)
    JuMPsolve(direction, m, sp)
end
hasnoises(sp::JuMP.Model) = length(ext(sp).noises) > 0

function JuMPsolve{T<:IterationDirection}(::Type{T}, ::SDDPModel, sp::JuMP.Model)
    @assert JuMP.solve(sp) == :Optimal
end

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
    Deserialize a serialized model. Not guaranteed to work or remain supported.
    The Base serializer is subject to change at any point. Use the JLD package
    if you want to save objects long-term.
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
