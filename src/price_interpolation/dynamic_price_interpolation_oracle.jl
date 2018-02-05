"""
    Since these oracles are slightly different...
"""
abstract type DynamicCutOracle <: AbstractCutOracle end

struct DefaultDynamicOracle{T} <: DynamicCutOracle
    cuts::Vector{Tuple{Cut, T}}
end
DefaultDynamicOracle(T) = DefaultDynamicOracle(Tuple{Cut, T}[])

function storecut!(oracle::DefaultDynamicOracle{T}, m::SDDPModel, sp::JuMP.Model, cut::Cut, price::T) where T
    push!(oracle.cuts, (cut, price))
end

function validcuts(oracle::DefaultDynamicOracle)
    oracle.cuts
end

"""
    NanniciniOracle(ρ::Int)

# Description

As proposed in

    Nannicini, G., Traversi, E., and Calvo, R. (2017). A Benders Squared (B2)
    Framework for Infinite-Horizon Stochastic Linear Programs. Optimization
    Online. http://www.optimization-online.org/DB_HTML/2017/06/6101.html.

This oracle permanently removes a cut is it is not used in the last `ρ` iterations.
"""
mutable struct NanniciniOracle{T} <: DynamicCutOracle
    ρ::Int
    cutsinmodel::Int
    cuts::Vector{Tuple{Cut, T}}
    iterations_since_last_active::Vector{Int}
end
NanniciniOracle(T,ρ::Int=typemax(Int)) = NanniciniOracle(ρ,0,Tuple{Cut, T}[], Int[])

function storecut!(oracle::NanniciniOracle, m::SDDPModel, sp::JuMP.Model, cut::Cut, price)
    push!(oracle.cuts, (cut, price))
    push!(oracle.iterations_since_last_active, 0)
    oracle.cutsinmodel += 1
end

function validcuts(oracle::NanniciniOracle)
    # sort in order of least used, to most recently used
    p = sortperm(oracle.iterations_since_last_active, rev=true)
    permute!(oracle.cuts, p)
    permute!(oracle.iterations_since_last_active, p)
    # find cut off
    idx = findfirst(x->x < oracle.ρ, oracle.iterations_since_last_active)
    if idx == 0
        error("No cuts in model used in the last $(oracle.ρ) iterations.")
    end
    oracle.cutsinmodel = length(oracle.cuts) - idx + 1
    oracle.cuts[idx:end]
end
