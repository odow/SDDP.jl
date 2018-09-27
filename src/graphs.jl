struct Graph{T}
    root_node::T
    nodes::Vector{T}
    edges::Dict{Tuple{T, T}, Float64}
end

function Graph(
    root_node::T, nodes::Vector{T}, edges::Vector{Pair{Tuple{T, T}, Float64}}) where T
    edge_dict = Dict{Tuple{T, T}, Float64}()
    for edge in edges
        edge_dict[edge.first] = edge.second
    end
    return Graph(root_node, nodes, edges)
end

function Graph(root_node::T) where T
    return Graph(root_node, T[], Dict{Tuple{T, T}, Float64}())
end

"""
    add_node(graph::Graph{T}, node::T) where T

Add a node to the graph `graph`.
"""
function add_node(graph::Graph{T}, node::T) where T
    if node in graph.nodes || node == graph.root_node
        error("Node $(node) already exists!")
    end
    push!(graph.nodes, node)
    return
end

"""
    add_node(graph::Graph{T}, node::T) where T

Add an edge to the graph `graph`.
"""
function add_edge(graph::Graph{T}, edge::Pair{Tuple{T, T}, Float64}) where T
    if haskey(graph.edges, edge.first)
        error("Edge $(edge.first) already exists!")
    elseif !(edge.first[1] == graph.root_node || edge.first[1] in graph.nodes)
        error("Node $(edge.first[1]) does not exist.")
    elseif !(edge.first[2] in graph.nodes)
        error("Node $(edge.first[2]) does not exist.")
    elseif edge.first[2] == graph.root_node
        error("Cannot have an edge entering the root node.")
    else
        graph.edges[edge.first] = edge.second
    end
    return
end

"""
    LinearGraph(stages::Int)
"""
function LinearGraph(stages::Int)
    node_type = Int
    root_node = 0
    nodes = collect(1:stages)
    edges = Dict{Tuple{node_type, node_type}, Float64}()
    for t in 1:stages
        edges[(t - 1, t)] = 1.0
    end
    return Graph(root_node, nodes, edges)
end

"""
    MarkovianGraph(transition_matrices::Vector{Matrix{Float64}})
"""
function MarkovianGraph(transition_matrices::Vector{Matrix{Float64}})
    if size(transition_matrices[1], 1) != 1
        error("Expected the first transition matrix to be of size (1, N). It " *
              "is of size $(size(transition_matrices[1])).")
    end
    node_type = Tuple{Int, Int}
    root_node = (0, 1)
    nodes = node_type[]
    edges = Dict{Tuple{node_type, node_type}, Float64}()
    for (stage, transition) in enumerate(transition_matrices)
        if !all(transition .>= 0.0)
            error("Entries in the transition matrix must be non-negative.")
        end
        if !all(0.0 .<= sum(transition; dims=2) .<= 1.0)
            error("Rows in the transition matrix must sum to between [0.0, 1.0].")
        end
        if stage > 1
            if size(transition_matrices[stage-1], 2) != size(transition, 1)
                error("Transition matrix for stage $(stage) is the wrong size.")
            end
        end
        for markov_state in 1:size(transition, 2)
            push!(nodes, (stage, markov_state))
        end
        for markov_state in 1:size(transition, 2)
            for last_markov_state in 1:size(transition, 1)
                probability = transition[last_markov_state, markov_state]
                if 0.0 < probability <= 1.0
                    edges[(stage - 1, last_markov_state), (stage, markov_state)] = probability
                end
            end
        end
    end
    return Graph(root_node, nodes, edges)
end

"""
    MarkovianGraph(; stages::Int,
                   transition_matrix::Matrix{Float64},
                   root_node_transition::Vector{Float64})

Construct a Markovian graph object.
"""
function MarkovianGraph(; stages::Int = 1,
                        transition_matrix::Matrix{Float64}=[1.0],
                        root_node_transition::Vector{Float64}=[1.0])
    return MarkovianGraph(
        vcat([reshape(root_node_transition, 1, length(root_node_transition))],
             [transition_matrix for stage in 1:(stages - 1)])
    )
end
