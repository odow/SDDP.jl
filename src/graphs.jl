struct Graph{T}
    # The root node of the policy graph.
    root_node::T
    # nodes[x] returns a vector of the children of node x and their
    # probabilities.
    nodes::Dict{T, Vector{Tuple{T, Float64}}}
end

function validate_graph(graph)
    for (node, children) in nodes
        probability = sum(child[2] for child in children)
        if !(0.0 <= probability <= 1.0)
            error("Probability on edges leaving node $(node) sum to " *
                  "$(probability), but this must be in [0.0, 1.0]")
        end
    end
end

function Graph(root_node::T) where T
    return Graph(root_node, Dict{T, Vector{Tuple{T, Float64}}}(
        root_node => Tuple{T, Float64}[]))
end

"""
    add_node(graph::Graph{T}, node::T) where T

Add a node to the graph `graph`.
"""
function add_node(graph::Graph{T}, node::T) where T
    if haskey(graph.nodes, node) || node == graph.root_node
        error("Node $(node) already exists!")
    end
    graph.nodes[node] = Tuple{T, Float64}[]
    return
end
function add_node(graph::Graph{T}, node) where T
    error("Unable to add node $(node). Nodes must be of type $(T).")
end

"""
    add_node(graph::Graph{T}, node::T) where T

Add an edge to the graph `graph`.
"""
function add_edge(graph::Graph{T}, edge::Pair{T, T}, probability::Float64) where T
    (parent, child) = edge
    if !(parent == graph.root_node || haskey(graph.nodes, parent))
        error("Node $(parent) does not exist.")
    elseif !haskey(graph.nodes, child)
        error("Node $(child) does not exist.")
    elseif child == graph.root_node
        error("Cannot have an edge entering the root node.")
    else
        push!(graph.nodes[parent], (child, probability))
    end
    return
end

function Graph(root_node::T, nodes::Vector{T},
               edges::Vector{Tuple{Pair{T, T}, Float64}}) where T
    graph = Graph(root_node)
    add_node.(Ref(graph), nodes)
    for (edge, probability) in edges
        add_edge(graph, edge, probability)
    end
    return graph
end

"""
    LinearGraph(stages::Int)
"""
function LinearGraph(stages::Int)
    edges = Tuple{Pair{Int, Int}, Float64}[]
    for t in 1:stages
        push!(edges, (t - 1 => t, 1.0))
    end
    return Graph(0, collect(1:stages), edges)
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
    edges = Tuple{Pair{node_type, node_type}, Float64}[]
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
                    push!(edges, (
                        (stage - 1, last_markov_state) => (stage, markov_state),
                        probability
                    ))
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
