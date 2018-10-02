abstract type AbstractSamplingScheme end

"""
    sample_scenario(graph::PolicyGraph{T}, ::AbstractSamplingScheme;
                    max_cycles::Int = 1) where T

Sample a scenario from the policy graph `graph`. If the graph is cyclic, return
once `max_cycles` number of cycles have been detected.

The scenario is a list of tuples (type `Vector{Tuple{T, Any}}`) where the first
component of each tuple is the index of the node, and the second component is
the stagewise-independent noise term observed in that node.
"""
function sample_scenario(graph::PolicyGraph{T},
                         sampling_scheme::AbstractSamplingScheme;
                         max_cycles::Int = 1) where T
    error("You need to overload the function Kokako.sample_scenario for the " *
          "sampling scheme (sampling_scheme).")
end

# ========================= Monte Carlo Sampling Scheme ========================

struct MonteCarlo end

# A helper utility for sampling a Noise using Monte Carlo.
function sample_noise(::MonteCarlo, noise_terms::Vector{<:Noise})
    if length(noise_terms) == 0
        return nothing
    end
    cumulative_probability = sum(noise.probability for noise in noise_terms)
    if cumulative_probability > 1.0
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
          " using Kokako.MonteCarlo().")
end

function sample_scenario(graph::PolicyGraph{T}, sampling_scheme::MonteCarlo;
                         max_cycles::Int = 1) where T
    scenario_path = Tuple{T, Any}[]
    visited_nodes = Set{T}()
    node_index = sample_noise(sampling_scheme, graph.root_children)::T
    while true
        node = graph[node_index]
        noise = sample_noise(sampling_scheme, node.noise_terms)
        push!(scenario_path, (node_index, noise))
        if length(node.children) == 0
            break  # We are at a leaf node.
        elseif node_index in visited_nodes
            # We are at the start of a cycle. Note that we still sample a noise
            # from this node for completeness.
            max_cycles -= 1
            if max_cycles <= 0
                break
            end
        end
        push!(visited_nodes, node_index)
        node_index = sample_noise(sampling_scheme, node.children)::T
    end
    return scenario_path
end
