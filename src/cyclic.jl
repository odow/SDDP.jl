# This source code is licensed under the Creative Commons ShareAlike license 3.
# For more details, see:
# https://en.wikipedia.org/wiki/Wikipedia:Text_of_Creative_Commons_Attribution-ShareAlike_3.0_Unported_License

"""
    is_cyclic(G::PolicyGraph{T}) where {T}

Return `true` or `false` if the graph `G` contains a cycle.

We implement Tarjan's strongly connected components algorithm to detect cycles
in a directed graph in O(|V| + |E|) time. See this Wiki for details
https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
The notation here follows the pseudocode in the Wikipedia article, rather than
the typical JuMP style guide.

Since we're only checking for cyclic graphs, we can stop as soon as one is
found. A cyclic graph has a stongly connected component with at least two
components, or it has a node with connects to itself. That means we don't need
to store the set of all strongly connected components.
"""
function is_cyclic(G::PolicyGraph{T}) where {T}
    index_counter = 0
    S = T[]
    low_link = Dict{T,Int}()
    index = Dict{T,Int}()
    on_stack = Dict{T,Bool}()
    function strong_connect(v)
        index[v] = index_counter
        low_link[v] = index_counter
        index_counter += 1
        push!(S, v)
        on_stack[v] = true
        for child in G[v].children
            w = child.term
            if v == w
                # Cycle detected: Type I: a node that loops to itself.
                return true
            end
            if !haskey(index, w)
                if strong_connect(w)
                    # A cycle was detected further down the tree. Propogate it
                    # upwards.
                    return true
                end
                low_link[v] = min(low_link[v], low_link[w])
            elseif on_stack[w]
                low_link[v] = min(low_link[v], index[w])
            end
        end
        if low_link[v] == index[v]
            scc = T[]
            w = G.root_node
            while v != w
                w = pop!(S)
                on_stack[w] = false
                push!(scc, w)
            end
            if length(scc) > 1
                # Cycle detected: Type II: a strongly connected component with
                # more than one element.
                return true
            end
        end
        return false  # No cycle detected.
    end
    for v in keys(G.nodes)
        if !haskey(index, v)
            if strong_connect(v)
                # Cycle detected!
                return true
            end
        end
    end
    return false
end
