#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export DynamicPriceInterpolation, StaticPriceInterpolation,
    DiscreteDistribution, observation, probability

include("discrete_distribution.jl")

#=
    SDP hybrid of Gjelsvik et al.
=#

mutable struct StaticPriceInterpolation{C<:AbstractCutOracle, T, T2} <: AbstractValueFunction
    initial_price::T
    location::T
    rib_locations::Vector{T}
    variables::Vector{JuMP.Variable}
    cutoracles::Vector{C}
    noises::DiscreteDistribution{T2}
    objective::Function
    dynamics::Function
end

#=
    Baucke's method
=#
include("dynamic_price_interpolation_oracle.jl")
mutable struct DynamicPriceInterpolation{C<:DynamicCutOracle, T, T2} <: AbstractValueFunction
    initial_price::T
    location::T
    minprice::T
    maxprice::T
    noises::DiscreteDistribution{T2}
    objective::Function
    dynamics::Function
    mu::Vector{JuMP.Variable}
    oracle::C
    lipschitz_constant::Float64
    bound::Float64
end

# methods used by both types

const PriceInterpolationMethods = Union{DynamicPriceInterpolation, StaticPriceInterpolation}

function getstageobjective(vf::PriceInterpolationMethods, sp::JuMP.Model)
    if ext(sp).finalstage
        return getobjectivevalue(sp)
    else
        return getobjectivevalue(sp) - getvalue(interpolate(vf))
    end
end

function setstageobjective!(vf::PriceInterpolationMethods, sp::JuMP.Model, obj::Function)
    vf.objective = obj
end

function writeasynccut!{T}(io, cut::Tuple{Int, Int, T, Cut})
    writecut!(io, cut...)#[1], cut[2], cut[3], cut[4])
end

# ==============================================================================

samplepricenoise{T}(stage::Int, noises::DiscreteDistribution{T}, solutionstore::Void) = sample(noises)

function samplepricenoise{T}(stage::Int, noises::DiscreteDistribution{T}, solutionstore::Dict{Symbol, Any})
    if haskey(solutionstore, :pricenoise)
        @assert length(solutionstore[:pricenoise])>=stage
        return observation(getnoise(noises, solutionstore[:pricenoise][stage]))
    else
        return sample(noises)
    end
end

function solvesubproblem!{V<:PriceInterpolationMethods}(::Type{ForwardPass}, m::SDDPModel{V}, sp::JuMP.Model, solutionstore)
    vf = valueoracle(sp)
    if ext(sp).stage == 1
        vf.location = vf.initial_price
    end
    p = vf.location
    # learn noise
    w = samplepricenoise(ext(sp).stage, vf.noises, solutionstore)

    # update price
    setobjective!(sp, p, w)
    passpriceforward!(m, sp)
    JuMPsolve(ForwardPass, m, sp)
end

function setobjective!(sp::JuMP.Model, price, noise)
    vf = valueoracle(sp)
    p = vf.dynamics(price, noise, ext(sp).stage, ext(sp).markovstate)
    vf.location = p
    # stage objective obj
    stageobj = vf.objective(p)

    # future objective
    future_value = interpolate(vf)
    # set
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), stageobj)
    else
        JuMP.setobjective(sp, getsense(sp), stageobj + future_value)
    end
end

function passpriceforward!{V<:PriceInterpolationMethods}(m::SDDPModel{V}, sp::JuMP.Model)
    stage = ext(sp).stage
    if stage < length(m.stages)
        # pass price forward
        for sp2 in subproblems(m, stage + 1)
            valueoracle(sp2).location = valueoracle(sp).location
        end
    end
end

function storekey!(::Type{Val{:price}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    push!(store, valueoracle(sp).location)
end

function storekey!(::Type{Val{:pricenoise}}, store, markov::Int, noiseidx::Int, sp::JuMP.Model, t::Int)
    nothing
end

function solvepricenoises!(m::SDDPModel, sp::JuMP.Model, last_markov_state, price)
    vf = valueoracle(sp)
    ex = ext(sp)
    markov_prob = getstage(m, ex.stage).transitionprobabilities[last_markov_state, ex.markovstate]
    for price_noise in vf.noises
        setobjective!(sp, price, observation(price_noise))
        solvesubproblem!(BackwardPass, m, sp, markov_prob * probability(price_noise))
    end
end

function constructcut{V<:PriceInterpolationMethods}(m::SDDPModel{V}, sp::JuMP.Model, ex, t, price)
    reset!(m.storage)
    for sp2 in subproblems(m, t+1)
        setstates!(m, sp2)
        solvepricenoises!(m, sp2, ex.markovstate, price)
    end
    # add cut
    I = 1:length(m.storage.objective)
    modifyprobability!(ex.riskmeasure,
        view(m.storage.modifiedprobability.data, I),
        m.storage.probability.data[I],
        m.storage.objective.data[I],
        m,
        sp
    )
    cut = constructcut(m, sp)
end

function calculatefirststagebound(m::SDDPModel)
    reset!(m.storage)
    for sp in subproblems(m, 1)
        vf = valueoracle(sp)
        solvepricenoises!(m, sp, 1, vf.initial_price)
    end
    dot(m.storage.objective, m.storage.probability)
end

function backwardpass!{V<:PriceInterpolationMethods}(m::SDDPModel{V}, settings::Settings)
    for t in (nstages(m)-1):-1:1
        for sp in subproblems(m, t)
            updatevaluefunction!(m, settings, t, sp)
        end
    end
    calculatefirststagebound(m)
end

function writecut!(filename::String, stage::Int, markovstate::Int, price, cut::Cut)
    open(filename, "a") do file
        writecut!(file, stage, markovstate, price, cut)
    end
end

function writecut!(io::IO, stage::Int, markovstate::Int, price, cut::Cut)
    write(io, string(stage), ",", string(price), ",", string(markovstate), ",", string(cut.intercept))
    for pi in cut.coefficients
        write(io, ",", string(pi))
    end
    write(io, "\n")
end

function simulate{V<:PriceInterpolationMethods}(m::SDDPModel{V},
        variables::Vector{Symbol}         = Symbol[];
        noises::AbstractVector{Int}       = zeros(Int, length(m.stages)),
        markovstates::AbstractVector{Int} = ones(Int, length(m.stages)),
        pricenoises::AbstractVector{Int}  = zeros(Int, length(m.stages))
    )
    store = newsolutionstore(variables)
    store[:pricenoise] = Int[]
    for t in 1:length(m.stages)
        push!(store[:markov], markovstates[t])
        push!(store[:noise], noises[t])
        push!(store[:pricenoise], pricenoises[t])
    end
    obj = forwardpass!(m, Settings(), store)
    store[:objective] = obj
    return store
end

# ==============================================================================
#   loadcuts!
function parsesinglepriceline(line)
    # stage, price, markovstate, intercept, coefficients...
    items = split(line, ",")
    stage = parse(Int, items[1])
    price = parse(Float64, items[2])
    markov_state = parse(Int, items[3])
    intercept = parse(Float64, items[4])
    coefficients = [parse(Float64, i) for i in items[5:end]]
    cut = Cut(intercept, coefficients)
    return stage, markov_state, price, cut
end
function parsemultipriceline(pricebegin, line)
    # stage, (price1, price2, ...), markovstate, intercept, coefficients...
    stage = parse(Int, split(line[1:pricebegin], ",")[1])
    priceend = findfirst(line, ')')
    price_items = split(line[pricebegin+1:priceend-1], ",")
    price = tuple([parse(Float64, p) for p in price_items]...)
    items = split(line[priceend+2:end], ",")
    ms = parse(Int, items[1])
    intercept = parse(Float64, items[2])
    coefficients = [parse(Float64, i) for i in items[3:end]]
    cut = Cut(intercept, coefficients)
    return stage, markov_state, price, cut
end

function loadcuts!{V<:PriceInterpolationMethods}(m::SDDPModel{V}, filename::String)
    open(filename, "r") do file
        while true
            line      = readline(file)
            line == nothing || line == "" && break
            pricebegin = findfirst(line, '(')
            if pricebegin == 0
                addasynccut!(m, parsesinglepriceline(line))
            else
                addasynccut!(m, parsemultipriceline(pricebegin, line))
            end
        end
    end
end

include("static_price_interpolation.jl")
include("dynamic_price_interpolation.jl")
