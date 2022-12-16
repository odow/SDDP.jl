#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# ========================= Monte Carlo Sampling Scheme ====================== #

struct InSampleMonteCarlo <: AbstractSamplingScheme
    max_depth::Int
    terminate_on_cycle::Bool
    terminate_on_dummy_leaf::Bool
    rollout_limit::Function
    initial_node::Any
end

"""
    InSampleMonteCarlo(;
        max_depth::Int = 0,
        terminate_on_cycle::Function = false,
        terminate_on_dummy_leaf::Function = true,
        rollout_limit::Function = (i::Int) -> typemax(Int),
        initial_node::Any = nothing,
    )

A Monte Carlo sampling scheme using the in-sample data from the policy graph
definition.

If `terminate_on_cycle`, terminate the forward pass once a cycle is detected.
If `max_depth > 0`, return once `max_depth` nodes have been sampled.
If `terminate_on_dummy_leaf`, terminate the forward pass with 1 - probability of
sampling a child node.

Note that if `terminate_on_cycle = false` and `terminate_on_dummy_leaf = false`
then `max_depth` must be set > 0.

Control which node the trajectories start from using `initial_node`. If it is
left as `nothing`, the root node is used as the starting node.

You can use `rollout_limit` to set iteration specific depth limits. For example:

    InSampleMonteCarlo(rollout_limit = i -> 2 * i)
"""
function InSampleMonteCarlo(;
    max_depth::Int = 0,
    terminate_on_cycle::Bool = false,
    terminate_on_dummy_leaf::Bool = true,
    rollout_limit::Function = i -> typemax(Int),
    initial_node::Any = nothing,
)
    if !terminate_on_cycle && !terminate_on_dummy_leaf && max_depth == 0
        error(
            "terminate_on_cycle and terminate_on_dummy_leaf cannot both be " *
            "false when max_depth=0.",
        )
    end
    new_rollout = let i = 0
        () -> (i += 1; rollout_limit(i))
    end
    return InSampleMonteCarlo(
        max_depth,
        terminate_on_cycle,
        terminate_on_dummy_leaf,
        new_rollout,
        initial_node,
    )
end

# ==================== OutOfSampleMonteCarlo Sampling Scheme ================= #

struct OutOfSampleMonteCarlo{T} <: AbstractSamplingScheme
    noise_terms::Dict{T,Vector{Noise}}
    root_children::Vector{Noise{T}}
    children::Dict{T,Vector{Noise{T}}}
    terminate_on_cycle::Bool
    terminate_on_dummy_leaf::Bool
    max_depth::Int
    rollout_limit::Function
    initial_node::Union{Nothing,T}
end

"""
    OutOfSampleMonteCarlo(
        f::Function,
        graph::PolicyGraph;
        use_insample_transition::Bool = false,
        max_depth::Int = 0,
        terminate_on_cycle::Bool = false,
        terminate_on_dummy_leaf::Bool = true,
        rollout_limit::Function = i -> typemax(Int),
        initial_node = nothing,
    )

Create a Monte Carlo sampler using out-of-sample probabilities and/or supports
for the stagewise-independent noise terms, and out-of-sample probabilities for
the node-transition matrix.

`f` is a function that takes the name of a node and returns a tuple containing
a vector of new [`SDDP.Noise`](@ref) terms for the children of that node, and
a vector of new [`SDDP.Noise`](@ref) terms for the stagewise-independent
noise.

If `f` is called with the name of the root node (e.g., `0` in a linear policy
graph, `(0, 1)` in a Markovian Policy Graph), then return a vector of
[`SDDP.Noise`](@ref) for the children of the root node.

If `use_insample_transition`, the in-sample transition probabilities will be
used. Therefore, `f` should only return a vector of the stagewise-independent
noise terms, and `f` will not be called for the root node.

If `terminate_on_cycle`, terminate the forward pass once a cycle is detected.
If `max_depth > 0`, return once `max_depth` nodes have been sampled.
If `terminate_on_dummy_leaf`, terminate the forward pass with 1 - probability of
sampling a child node.

Note that if `terminate_on_cycle = false` and `terminate_on_dummy_leaf = false`
then `max_depth` must be set > 0.

Control which node the trajectories start from using `initial_node`. If it is
left as `nothing`, the root node is used as the starting node.

You can use `rollout_limit` to set iteration specific depth limits. For example:

```julia
OutOfSampleMonteCarlo(rollout_limit = i -> 2 * i)
```

## Examples

Given linear policy graph `graph` with `T` stages:
```julia
sampler = OutOfSampleMonteCarlo(graph) do node
    if node == 0
        return [SDDP.Noise(1, 1.0)]
    else
        noise_terms = [SDDP.Noise(node, 0.3), SDDP.Noise(node + 1, 0.7)]
        children = node < T ? [SDDP.Noise(node + 1, 0.9)] : SDDP.Noise{Int}[]
        return children, noise_terms
    end
end
```

Given linear policy graph `graph` with `T` stages:
```julia
sampler = OutOfSampleMonteCarlo(graph, use_insample_transition=true) do node
    return [SDDP.Noise(node, 0.3), SDDP.Noise(node + 1, 0.7)]
end
```
"""
function OutOfSampleMonteCarlo(
    f::Function,
    graph::PolicyGraph{T};
    use_insample_transition::Bool = false,
    max_depth::Int = 0,
    terminate_on_cycle::Bool = false,
    terminate_on_dummy_leaf::Bool = true,
    rollout_limit::Function = i -> typemax(Int),
    initial_node::Union{Nothing,T} = nothing,
) where {T}
    if !terminate_on_cycle && !terminate_on_dummy_leaf && max_depth == 0
        error(
            "terminate_on_cycle and terminate_on_dummy_leaf cannot both be " *
            "false when max_depth=0.",
        )
    end
    noise_terms = Dict{T,Vector{Noise}}()
    children = Dict{T,Vector{Noise{T}}}()
    root_children = if use_insample_transition
        graph.root_children
    else
        f(graph.root_node)::Vector{Noise{T}}
    end
    for key in keys(graph.nodes)
        if use_insample_transition
            child = graph.nodes[key].children
            noise = f(key)
        else
            child, noise = f(key)
        end
        noise_terms[key] = convert(Vector{Noise}, noise)
        children[key] = child
    end
    new_rollout = let i = 0
        () -> (i += 1; rollout_limit(i))
    end
    return OutOfSampleMonteCarlo{T}(
        noise_terms,
        root_children,
        children,
        terminate_on_cycle,
        terminate_on_dummy_leaf,
        max_depth,
        new_rollout,
        initial_node,
    )
end

function get_noise_terms(
    sampling_scheme::InSampleMonteCarlo,
    node::Node{T},
    node_index::T,
) where {T}
    return node.noise_terms
end

function get_noise_terms(
    sampling_scheme::OutOfSampleMonteCarlo{T},
    node::Node{T},
    node_index::T,
) where {T}
    return sampling_scheme.noise_terms[node_index]
end

function get_children(
    sampling_scheme::InSampleMonteCarlo,
    node::Node{T},
    node_index::T,
) where {T}
    return node.children
end

function get_children(
    sampling_scheme::OutOfSampleMonteCarlo{T},
    node::Node{T},
    node_index::T,
) where {T}
    return sampling_scheme.children[node_index]
end

function get_root_children(
    sampling_scheme::InSampleMonteCarlo,
    graph::PolicyGraph{T},
) where {T}
    return graph.root_children
end

function get_root_children(
    sampling_scheme::OutOfSampleMonteCarlo{T},
    graph::PolicyGraph{T},
) where {T}
    return sampling_scheme.root_children
end

function sample_noise(noise_terms::Vector{<:Noise})
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
    return error(
        "Internal SDDP error: unable to sample noise from $(noise_terms)",
    )
end

function sample_scenario(
    graph::PolicyGraph{T},
    sampling_scheme::Union{InSampleMonteCarlo,OutOfSampleMonteCarlo{T}},
) where {T}
    max_depth = min(sampling_scheme.max_depth, sampling_scheme.rollout_limit())
    # Storage for our scenario. Each tuple is (node_index, noise.term).
    scenario_path = Tuple{T,Any}[]
    # We only use visited_nodes if terminate_on_cycle=true. Just initialize
    # anyway.
    visited_nodes = Set{T}()
    # Begin by sampling a node from the children of the root node.
    node_index = something(
        sampling_scheme.initial_node,
        sample_noise(get_root_children(sampling_scheme, graph)),
    )::T
    while true
        node = graph[node_index]
        noise_terms = get_noise_terms(sampling_scheme, node, node_index)
        children = get_children(sampling_scheme, node, node_index)
        noise = sample_noise(noise_terms)
        push!(scenario_path, (node_index, noise))
        # Termination conditions:
        if length(children) == 0
            # 1. Our node has no children, i.e., we are at a leaf node.
            return scenario_path, false
        elseif sampling_scheme.terminate_on_cycle && node_index in visited_nodes
            # 2. terminate_on_cycle = true and we have detected a cycle.
            return scenario_path, true
        elseif 0 < sampling_scheme.max_depth <= length(scenario_path)
            # 3. max_depth > 0 and we have explored max_depth number of nodes.
            return scenario_path, false
        elseif sampling_scheme.terminate_on_dummy_leaf &&
               rand() < 1 - sum(child.probability for child in children)
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
        node_index = sample_noise(children)::T
    end
    # Throw an error because we should never end up here.
    return error(
        "Internal SDDP error: something went wrong sampling a scenario.",
    )
end

# ========================= Historical Sampling Scheme ======================= #

mutable struct Historical{T,S} <: AbstractSamplingScheme
    scenarios::Vector{Noise{Vector{Tuple{T,S}}}}
    sequential::Bool
    counter::Int
end

function Base.show(io::IO, h::Historical)
    print(
        io,
        "A Historical sampler with $(length(h.scenarios)) scenarios sampled ",
        h.sequential ? "sequentially." : "probabilistically.",
    )
    return
end

"""
    Historical(
        scenarios::Vector{Vector{Tuple{T,S}}},
        probability::Vector{Float64},
    ) where {T,S}

A sampling scheme that samples a scenario from the vector of scenarios
`scenarios` according to `probability`.

## Examples

```julia
Historical(
    [
        [(1, 0.5), (2, 1.0), (3, 0.5)],
        [(1, 0.5), (2, 0.0), (3, 1.0)],
        [(1, 1.0), (2, 0.0), (3, 0.0)]
    ],
    [0.2, 0.5, 0.3],
)
```
"""
function Historical(
    scenarios::Vector{Vector{Tuple{T,S}}},
    probability::Vector{Float64},
) where {T,S}
    if !(sum(probability) ≈ 1.0)
        error(
            "Probability of historical scenarios must sum to 1. Currently: " *
            "$(sum(probability)).",
        )
    end
    output = [Noise(s, p) for (s, p) in zip(scenarios, probability)]
    return Historical(output, false, 0)
end

"""
    Historical(scenarios::Vector{Vector{Tuple{T,S}}}) where {T,S}

A deterministic sampling scheme that iterates through the vector of provided
`scenarios`.

## Examples

```julia
Historical([
    [(1, 0.5), (2, 1.0), (3, 0.5)],
    [(1, 0.5), (2, 0.0), (3, 1.0)],
    [(1, 1.0), (2, 0.0), (3, 0.0)],
])
```
"""
function Historical(scenarios::Vector{Vector{Tuple{T,S}}}) where {T,S}
    return Historical(Noise.(scenarios, NaN), true, 0)
end

"""
    Historical(scenario::Vector{Tuple{T,S}}) where {T,S}

A deterministic sampling scheme that always samples `scenario`.

## Examples

```julia
Historical([(1, 0.5), (2, 1.5), (3, 0.75)])
```
"""
Historical(scenario::Vector{Tuple{T,S}}) where {T,S} = Historical([scenario])

function sample_scenario(
    graph::PolicyGraph{T},
    sampling_scheme::Historical{T,NoiseTerm};
    # Ignore the other kwargs because the user is giving
    # us the full scenario.
    kwargs...,
) where {T,NoiseTerm}
    if sampling_scheme.sequential
        sampling_scheme.counter += 1
        if sampling_scheme.counter > length(sampling_scheme.scenarios)
            sampling_scheme.counter = 1
        end
        return sampling_scheme.scenarios[sampling_scheme.counter].term, false
    end
    return sample_noise(sampling_scheme.scenarios), false
end

"""
    PSRSamplingScheme(N::Int; sampling_scheme = InSampleMonteCarlo())

A sampling scheme with `N` scenarios, similar to how PSR does it.
"""
mutable struct PSRSamplingScheme{A} <: AbstractSamplingScheme
    N::Int
    sampling_scheme::A
    scenarios::Vector{Any}
    counter::Int

    function PSRSamplingScheme(
        N::Int;
        sampling_scheme::AbstractSamplingScheme = InSampleMonteCarlo(),
    )
        return new{typeof(sampling_scheme)}(N, sampling_scheme, Any[], 0)
    end
end

function Base.show(io::IO, h::PSRSamplingScheme)
    print(io, "A sampler with $(length(h.scenarios)) scenarios like PSR does.")
    return
end

function sample_scenario(
    graph::PolicyGraph{T},
    s::PSRSamplingScheme{A};
    kwargs...,
) where {T,A}
    s.counter += 1
    if s.counter > s.N
        s.counter = 1
    end
    if s.counter > length(s.scenarios)
        push!(s.scenarios, sample_scenario(graph, s.sampling_scheme; kwargs...))
    end
    return s.scenarios[s.counter]
end
