#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.


# ========================= Monte Carlo Sampling Scheme ====================== #

"""
    InSampleMonteCarlo(;
        terminate_on_cycle = false,
        max_depth = 0,
        terminate_on_dummy_leaf = true
    )

A Monte Carlo sampling scheme using the in-sample data from the policy graph
definition.

If `terminate_on_cycle`, terminate the forward pass once a cycle is detected.
If `max_depth > 0`, return once `max_depth` nodes have been sampled.
If `terminate_on_dummy_leaf`, terminate the forward pass with 1 - probability of
sampling a child node.

Note that if `terminate_on_cycle = false` and `terminate_on_dummy_leaf = false`
then `max_depth` must be set > 0.
"""
struct InSampleMonteCarlo <: AbstractSamplingScheme
    terminate_on_cycle::Bool
    terminate_on_dummy_leaf::Bool
    max_depth::Int
    function InSampleMonteCarlo(;
        terminate_on_cycle::Bool = false,
        terminate_on_dummy_leaf::Bool = true,
        max_depth::Int = 0)
        if !terminate_on_cycle && !terminate_on_dummy_leaf && max_depth == 0
            error("terminate_on_cycle and terminate_on_dummy_leaf cannot both" *
                  " be false when max_depth=0.")
        end
        return new(terminate_on_cycle, terminate_on_dummy_leaf, max_depth)
    end
end

# A helper utility for sampling a Noise using Monte Carlo.
function sample_noise(::InSampleMonteCarlo, noise_terms::Vector{<:Noise})
    if length(noise_terms) == 0
        return nothing
    end
    cumulative_probability = sum(noise.probability for noise in noise_terms)
    if cumulative_probability > 1.0 + 1e-6
        error("Cumulative probability cannot be greater than 1.0.")
    end
    rnd = rand() * cumulative_probability
    for noise in noise_terms
        rnd -= noise.probability
        if rnd <= 0.0
            return noise.term
        end
    end
    error("Internal Kokako error: unable to sample noise from $(noise_terms)" *
          " using Kokako.InSampleMonteCarlo().")
end

"""
    sample_scenario(graph, ::InSampleMonteCarlo; kwargs...)

Sample a scenario using the InSampleMonteCarlo sampler.
"""
function sample_scenario(graph::PolicyGraph{T},
                         sampling_scheme::InSampleMonteCarlo) where T
    # Storage for our scenario. Each tuple is (node_index, noise.term).
    scenario_path = Tuple{T, Any}[]
    # We only use visited_nodes if terminate_on_cycle=true. Just initialize
    # anyway.
    visited_nodes = Set{T}()
    # Begin by sampling a node from the children of the root node.
    node_index = sample_noise(sampling_scheme, graph.root_children)::T
    while true
        node = graph[node_index]
        noise = sample_noise(sampling_scheme, node.noise_terms)
        push!(scenario_path, (node_index, noise))
        # Termination conditions:
        if length(node.children) == 0
            # 1. Our node has no children, i.e., we are at a leaf node.
            return scenario_path, false
        elseif sampling_scheme.terminate_on_cycle && node_index in visited_nodes
            # 2. terminate_on_cycle = true and we have detected a cycle.
            return scenario_path, true
        elseif 0 < sampling_scheme.max_depth <= length(scenario_path)
            # 3. max_depth > 0 and we have explored max_depth number of nodes.
            return scenario_path, false
        elseif sampling_scheme.terminate_on_dummy_leaf &&
                rand() < 1 - sum(child.probability for child in node.children)
            # 4. we sample a "dummy" leaf node in the next step due to the
            # probability of the child nodes summing to less than one.
            return scenario_path, false
        end
        # We only need to store a list of visited nodes if we want to terminate
        # due to the presence of a cycle.
        if sampling_scheme.terminate_on_cycle
            push!(visited_nodes, node_index)
        end
        # Sample a new node to transition to.
        node_index = sample_noise(sampling_scheme, node.children)::T
    end
    # Throw an error because we should never end up here.
    error("Internal Kokako error: something went wrong sampling a scenario.")
end

# ========================= Historical Sampling Scheme ======================= #

struct Historical{NodeIndex, NoiseTerm} <: AbstractSamplingScheme
    scenarios::Vector{Noise{Tuple{NodeIndex, NoiseTerm}}}
end

"""
    Historical(scenarios::Vector{Vector{Tuple{T, S}}},
               probability::Vector{Float64})

A sampling scheme that samples a scenario from the vector of scenarios
`scenarios` according to `probability`. If probability omitted, defaults to
uniform probability.

# Example

    Historical(
        [
            [(1, 0.5), (2, 1.0), (3, 0.5)],
            [(1, 0.5), (2, 0.0), (3, 1.0)],
            [(1, 1.0), (2, 0.0), (3, 0.0)]
        ],
        [0.2, 0.5, 0.3]
    )
"""
function Historical(scenarios::Vector{Vector{Tuple{NodeIndex, NoiseTerm}}},
                    probability::Vector{Float64} =
                        fill(1.0 / length(scenarios), length(scenarios))
                    ) where {NodeIndex, NoiseTerm}
    if sum(probability) != 1.0
        error("Probability of historical scenarios must sum to 1. Currently: " *
              "$(sum(probability)).")
    end
    output = Noise{Vector{Tuple{NodeIndex, NoiseTerm}}}[]
    for (scenario, prob) in zip(scenarios, probability)
        push!(output, Noise(scenario, prob))
    end
    return Historical(output)
end

function sample_scenario(graph::PolicyGraph{T},
                         sampling_scheme::Historical{T, NoiseTerm};
                         # Ignore the other kwargs because the user is giving
                         # us the full scenario.
                         kwargs...) where {T, NoiseTerm}
    return sample_noise(InSampleMonteCarlo(), sampling_scheme.scenarios), false
end
