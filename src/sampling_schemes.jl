abstract type AbstractSamplingScheme end

"""
    sample_scenario(graph::PolicyGraph{T}, ::AbstractSamplingScheme;
                    terminate_on_cycle::Bool = true,
                    max_depth::Int = 0) where T

Sample a scenario from the policy graph `graph` based on the sampling scheme.

If `terminate_on_cycle`, terminate the forward pass once a cycle is detected.

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

struct MonteCarlo end

# A helper utility for sampling a Noise using Monte Carlo.
function sample_noise(::MonteCarlo, noise_terms::Vector{<:Noise})
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
          " using Kokako.MonteCarlo().")
end

function sample_scenario(graph::PolicyGraph{T}, sampling_scheme::MonteCarlo;
                         terminate_on_cycle::Bool = true,
                         max_depth::Int = 0) where T
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
            if terminate_on_cycle
                break
            elseif max_depth == 0
                error("If terminate_on_cycle=false, then max_depth must be >0.")
            elseif length(scenario_path) >= max_depth
                break
            end
        end
        push!(visited_nodes, node_index)
        node_index = sample_noise(sampling_scheme, node.children)::T
    end
    return scenario_path
end
