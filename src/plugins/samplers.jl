include("../user_interface.jl")
using StatsBase

abstract type AbstractBackwardPassSampler end

struct CompleteSampler <: AbstractBackwardPassSampler end

struct MonteCarloSampler <: AbstractBackwardPassSampler end

sample_backward_noise_terms(::CompleteSampler, node) = node.noise_terms

function sample_backward_noise_terms(::MonteCarloSampler, node)
    sampled = Noise[]
    terms = [noise.term for noise in node.noise_terms]

    sampled_t = sample(terms,10)
    u=unique(sampled_t)
    d=Dict([(i,count(x->x==i,sampled_t)) for i in u])

    for i in 1:length(u)
        push!(sampled,Noise(u[i],1/d[u[i]]))
    end
    return sampled
end
