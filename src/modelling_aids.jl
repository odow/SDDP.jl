#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
    MarkovianPolicyGraph(
        simulator::Function;
        stages::Int,
        simulations::Int,
        markov_states::Union{Int, Vector{Int}})

Create a discretized Markov chain representation with `stages` stages of the
stochastic process defined by the simulator function `simulator`.

In stage `t`, the process is discretized into `markov_states[t]` (unless
`markov_states::Int`).

`simulator` is a function that takes two arguments:
    1) the stage `t::Int`
    2) the previous state `x`
When `t=1`, `x=nothing`.

## Examples

    MarkovianGraph(stages = 4, simulations = 100, markov_states = 3) do t, x
        return t == 1 ? 0.0 : x + rand()
    end
"""
function MarkovianGraph(simulator::Function; stages::Int, n_simulations::Int,
                        n_markov_states::Union{Int, Vector{Int}})
    if n_markov_states isa Int
        n_markov_states = fill(n_markov_states, stages)
    end
    # Type of Markov state
    T = typeof(simulator(1, nothing))
    # Simulate the stochastic process.
    simulations = zeros(T, (stages, n_simulations))
    for i in 1:n_simulations
        x = nothing
        for t in 1:stages
            simulations[t, i] = x = simulator(t, x)
        end
    end
    # Fit the supports of the Markov process and the transition matrices.
    supports = _fit_markov_supports(simulations, n_markov_states)
    transition_matrices = _fit_probabilities(simulations, supports)
    # Turn this into a SDDP.Graph object.
    node_type, root_node = Tuple{Int, Float64}, (0, zero(T))
    nodes, edges = node_type[], Tuple{Pair{node_type, node_type}, Float64}[]
    for (t, transition) in enumerate(transition_matrices)
        for i in 1:size(transition, 2)
            push!(nodes, (t, supports[t][i]))
        end
        for m in 1:size(transition, 1)
            support = t == 1 ? zero(T) : supports[t - 1][m]
            for m′ in 1:size(transition, 2)
                support′ = supports[t][m′]
                transition[m, m′] == 0.0 && continue
                push!(edges, ((t - 1, support) => (t, support′), transition[m, m′]))
            end
        end
    end
    return Graph(root_node, nodes, edges)
end

function _fit_markov_supports(
        simulations::Array{T, 2}, n_markov_states::Vector{Int}) where {T}
    supports = Vector{T}[]
    for t in 1:size(simulations, 1)
        random_k_means = [_k_means(simulations[t, :], n_markov_states[t])
            for _ in 1:20]
        _, idx = findmin(map(x -> x[2], random_k_means))
        push!(supports, random_k_means[idx][1])
    end
    return supports
end

function _k_means(samples::Vector{T}, K::Int) where {T}
    centroids = rand(samples, K)
    converged = false
    assignments = zeros(Int, length(samples))
    assignments′ = copy(assignments)
    counts = zeros(Int, length(centroids))
    l₂penalty = 0.0
    for iteration in 1:100
        # Assignment step: assign each observation to the nearest centroid.
        _assign_to_centroid(samples, centroids, assignments)
        # Update step: move centroid to mean of each assigned observation.
        centroids .= zero(T)
        counts .= 0
        for (sample, assignment) in zip(samples, assignments)
            centroids[assignment] += sample
            counts[assignment] += 1
        end
        centroids ./= counts
        if assignments == assignments′
            break
        end
        copyto!(assignments′, assignments)
    end
    l₂penalty = 0.0
    for (sample, assignment) in zip(samples, assignments)
        l₂penalty += (sample - centroids[assignment])^2
    end
    sort!(centroids)
    return centroids, l₂penalty / length(samples)
end

function _assign_to_centroid(samples, centroids, assignments)
    for (i, sample) in enumerate(samples)
        best_assignment = 1
        best_distance = Inf
        for (j, centroid) in enumerate(centroids)
            distance = (centroid - sample)^2
            if distance < best_distance
                best_assignment, best_distance = j, distance
            end
        end
        assignments[i] = best_assignment
    end
end

function _fit_probabilities(simulations, supports)
    T, N = size(simulations)
    # Assign simulated points to their nearest support.
    assignments = zeros(Int, size(simulations))
    for t in 1:T
        _assign_to_centroid(
            view(simulations, t, :), supports[t], view(assignments, t, :))
    end
    # Initialize the storage of transition matrices.
    transition_matrices = Array{Float64, 2}[]
    push!(transition_matrices, zeros(1, length(supports[1])))
    for t in 2:T
        push!(transition_matrices,
            zeros(length(supports[t-1]), length(supports[t])))
    end
    # Sum number of points transitioning from Markov state m to Markov state m′.
    for i in 1:N
        for t in 2:T
        transition_matrices[t][1, assignments[t, i]] += 1
            m = assignments[t - 1, i]
            m′ = assignments[t, i]
            transition_matrices[t][m, m′] += 1
        end
    end
    # Normalize the counts into a proportion.
    for t in 1:T
        transition_matrices[t] ./= N
    end
    return transition_matrices
end
