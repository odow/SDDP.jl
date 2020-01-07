#  Copyright 2017-20, Oscar Dowson and contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
     CompleteSampler()

Backward sampler that returns all noises of the corresponding node.
"""
struct CompleteSampler <: AbstractBackwardSamplingScheme end

sample_backward_noise_terms(::CompleteSampler, node) = node.noise_terms


"""
     MonteCarloSampler(number_of_samples::Int)

Backward sampler that returns `number_of_samples` noises sampled with
replacement from noises on the corresponding node.
"""
struct MonteCarloSampler <: AbstractBackwardSamplingScheme
    number_of_samples::Int
end

function sample_backward_noise_terms(sampler::MonteCarloSampler, node::Node)
    prob = 1 / sampler.number_of_samples
    return [Noise(sample_noise(node.noise_terms), prob) for _ = 1:sampler.number_of_samples]
end
