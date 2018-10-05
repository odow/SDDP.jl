"""
    AbstractSamplingScheme

The abstract type for the sampling-scheme interface.

You need to define the following methods:
 - Kokako.sample_scenario
"""
abstract type AbstractSamplingScheme end

"""
    sample_scenario(graph::PolicyGraph{T}, ::AbstractSamplingScheme;
                    terminate_on_cycle::Bool = true,
                    max_depth::Int = 0) where T

Sample a scenario from the policy graph `graph` based on the sampling scheme.

If `terminate_on_cycle`, terminate the forward pass once a cycle is detected.

*Important note:* If a cycle is detected, the scenario should include the node
that forms the cycle. For example, if there is a single node in the graph that
loops onto itself, the returned scenario should contain two elements.

If `max_depth > 0`, return once `max_depth` nodes have been sampled, otherwise,
return once a leaf node is reached. Note that if `terminate_on_cycle=false` then
`max_depth` must be set > 0.

The scenario is a list of tuples (type `Vector{Tuple{T, Any}}`) where the first
component of each tuple is the index of the node, and the second component is
the stagewise-independent noise term observed in that node.
"""
function sample_scenario(graph::PolicyGraph{T},
                         sampling_scheme::AbstractSamplingScheme;
                         terminate_on_cycle::Bool = true,
                         max_depth::Int = 0) where T
    error("You need to overload the function Kokako.sample_scenario for the " *
          "sampling scheme (sampling_scheme).")
end

# ========================= Monte Carlo Sampling Scheme ========================

"""
    InSampleMonteCarlo

A Monte Carlo sampling scheme using the in-sample data from the policy graph
definition.
"""
struct InSampleMonteCarlo <: AbstractSamplingScheme end

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

function sample_scenario(graph::PolicyGraph{T},
                         sampling_scheme::InSampleMonteCarlo;
                         terminate_on_cycle::Bool = true,
                         max_depth::Int = 0) where T
    if !terminate_on_cycle && max_depth == 0
        error("If terminate_on_cycle=false, then max_depth must be >0.")
    end
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
        #   1. Our node has no children, i.e., we are at a leaf node.
        #   2. terminate_on_cycle = true and we have detected a cycle.
        #   3. max_depth > 0 and we have explored max_depth number of nodes.
        if length(node.children) == 0
            return scenario_path
        elseif terminate_on_cycle && node_index in visited_nodes
            return scenario_path
        elseif max_depth > 0 && length(scenario_path) > max_depth
            return scenario_path
        end
        # We only need to store a list of visited nodes if we want to terminate
        # due to the presence of a cycle.
        if terminate_on_cycle
            push!(visited_nodes, node_index)
        end
        # Sample a new node to transition to.
        node_index = sample_noise(sampling_scheme, node.children)::T
    end
    # Throw an error because we should never end up here.
    error("Internal Kokako error: something went wrong sampling a scenario.")
end

# ========================= Historical Sampling Scheme ========================

struct Historical{NodeIndex, NoiseTerm} <: AbstractSamplingScheme
    scenarios::Vector{Noise{Tuple{NodeIndex, NoiseTerm}}}
end

"""
    Historical(scenarios::Vector{Vector{Tuple{T, S}}}, probability::Vector{Float64})

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
                         terminate_on_cycle::Bool = true,
                         max_depth::Int = 0) where {T, NoiseTerm}
    sample_noise(InSampleMonteCarlo(), sampling_scheme.scenarios)
end
