import StatsBase


"""
     CompleteSampler()

 Backward sampler that returns all noises of the corresponding node.
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
    weights = [noise.probability for noise in node.noise_terms]

    sampled_t = StatsBase.sample(terms,StatsBase.Weights(weights),mcs.number_of_samples)
    u=unique(sampled_t)
    d=Dict([(i,count(x->x==i,sampled_t)) for i in u])

    for i in 1:length(u)
        push!(sampled,Noise(u[i],d[u[i]]/mcs.number_of_samples))
    end
    return sampled
end
