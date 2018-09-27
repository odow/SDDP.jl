struct Noise
    # The noise term.
    term  # TODO(odow): make this a concrete type?
    # The probability of sampling the noise term.
    probability::Float64
end

struct State
    # The incoming state variable in the subproblem.
    incoming::JuMP.VariableRef
    # The outgoing state variable in the subproblem.
    outgoing::JuMP.VariableRef
    # The "fishing" constraint so we can query the dual.
    fishing_dual  # TODO(odow): type this appropriately
end

struct Child
    # The child node.
    node::JuMP.Model
    # The probability of transitioning to the child node from the parent.
    probability::Float64
end

mutable struct NodeExtension
    # A vector of the child nodes.
    children::Vector{Child}
    # A vector of the discrete stagewise-independent noise terms.
    noise_terms::Vector{Noise}
    # A function parameterize(model::JuMP.Model, noise::T) that modifies the
    # JuMP model based on the observation of the noise.
    parameterize::Function  # TODO(odow): make this a concrete type?
    # A list of the state variables in the model.
    states::Vector{State}
end

struct PolicyGraph{T}
    # The value of the initial state variables.
    root_state::Vector{Float64}
    # Children of the root node. child => probability.
    root_children::Vector{Pair{T, Float64}}
    # All nodes in the graph.
    nodes::Dict{T, JuMP.Model}
    PolicyGraph(T) = new{T}(Float64[], Pair{T, Float64}[], Dict{T, JuMP.Model}())
end

function construct_subproblem(optimizer_factory, direct_mode::Bool)
    subproblem = if direct_mode
        instance = optimizer_factory.constructor(
            optimizer_factory.args...; optimizer_factory.kwargs...)
        JuMP.direct_model(instance)
    else
        JuMP.Model(optimizer_factory)
    end
    ext = NodeExtension(Child[], Noise[], (sp, index)->(), State[])
    subproblem.ext[:kokako] = ext
    return subproblem
end

function construct_subproblem(optimizer_factory::Nothing, direct_mode::Bool)
    if direct_mode
        error("You must specify an optimizer in the form:\n" *
              "    with_optimizer(Module.Opimizer, args...) if " *
              "direct_mode=true.")
    end
    subproblem = JuMP.Model()
    ext = NodeExtension(Child[], Noise[], (sp, index)->(), State[])
    subproblem.ext[:kokako] = ext
    return subproblem
end
extension(subproblem::JuMP.Model) = subproblem.ext[:kokako]::NodeExtension

function PolicyGraph(builder::Function, graph::Graph{T};
                     optimizer = nothing,
                     direct_mode = true,
                     sense = :Min) where T
    policy_graph = PolicyGraph(T)
    # Initialize nodes.
    for (node_index, children) in graph.nodes
        if node_index == graph.root_node
            continue
        end
        subproblem = construct_subproblem(optimizer, direct_mode)
        builder(subproblem, node_index)
        policy_graph.nodes[node_index] = subproblem
    end
    # Loop back through and add the arcs/children.
    for (node_index, children) in graph.nodes
        if node_index == graph.root_node
            continue
        end
        ext = extension(policy_graph.nodes[node_index])
        for (child, probability) in children
            push!(ext.children,
                Child(policy_graph.nodes[child], probability))
        end
    end
    # Add root nodes
    for (child, probability) in graph.nodes[graph.root_node]
        push!(policy_graph.root_children, child => probability)
    end
    return policy_graph
end
