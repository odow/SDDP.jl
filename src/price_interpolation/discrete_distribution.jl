#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

export DiscreteDistribution, observation, probability

"""
    A single realization of a noise in a DiscreteDistribution
"""
struct NoiseRealization{T}
    observation::T
    probability::Float64
end

"""
    observation(n::Noise)

Return the value that is observed when the noise `n` is sampled.
"""
observation(n::NoiseRealization) = n.observation

"""
    probability(n::Noise)

Return the probability of sampling the noise `n`.
"""
probability(n::NoiseRealization) = n.probability

struct DiscreteDistribution{T}
    noises::Vector{NoiseRealization{T}}
end

"""
    DiscreteDistribution{T}(observations::AbstractVector{T}, probabilities::AbstractVector{Float64})

Create a finite discrete distribution of `observations` supported with probability
`probabilities`.
"""
function DiscreteDistribution{T}(observations::AbstractVector{T}, probabilities::AbstractVector{Float64})
    if !isapprox(sum(probabilities), 1.0, atol=1e-6)
        error("Finite discrete distribution must sum to 1.0")
    end
    y = NoiseRealization{T}[]
    for (xi, pi) in zip(observations, probabilities)
        push!(y, NoiseRealization(xi, pi))
    end
    DiscreteDistribution(y)
end

"""
    DiscreteDistribution{T}(observations::AbstractVector{T})

Create a finite discrete distribution of `observations` supported with
uniform probability.
"""
function DiscreteDistribution{T}(observations::AbstractVector{T})
    DiscreteDistribution(observations, fill(1.0 / length(observations), length(observations)))
end

Base.start(d::DiscreteDistribution) = 1
Base.next(d::DiscreteDistribution, state) = (d.noises[state], state+1)
Base.done(d::DiscreteDistribution, state) = (state > length(d.noises))
Base.getindex(d::DiscreteDistribution, idx::Int) = d.noises[idx]
Base.length(d::DiscreteDistribution) = length(d.noises)

function sample{T}(d::DiscreteDistribution{T})
    r = rand()
    for noise in d
        @inbounds r -= probability(noise)
        if r < eps(Float64)
            return observation(noise)
        end
    end
    error("x must be a discrete probablity distribution that sums to one. sum= $(sum(x))")
end
