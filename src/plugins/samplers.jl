import StatsBase


"""
     CompleteSampler()

 Backward sampler that returns all noises of the cosrresponding node.

 """
struct CompleteSampler <: AbstractBackwardPassSampler end
sample_backward_noise_terms(::CompleteSampler, node) = node.noise_terms


"""
     MonteCarloSampler(number_of_samples ::Integer)

 Backward sampler that returns `number_of_samples` noises sampled
 with replacement from noises on the corresponding node.

 """
struct MonteCarloSampler <: AbstractBackwardPassSampler
    number_of_samples ::Integer
end


function sample_backward_noise_terms(mcs::MonteCarloSampler, node)
    sampled = Noise[]
    terms = [noise.term for noise in node.noise_terms]

    n_samples = min(length(node.noise_terms),mcs.number_of_samples)
    sampled_t = StatsBase.sample(terms,n_samples)
    u=unique(sampled_t)
    d=Dict([(i,count(x->x==i,sampled_t)) for i in u])

    for i in 1:length(u)
        push!(sampled,Noise(u[i],d[u[i]]/n_samples))
    end
    return sampled
end
