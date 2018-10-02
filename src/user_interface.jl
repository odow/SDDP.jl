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
function add_edge(graph::Graph{T}, edge::Pair{T, T},
                  probability::Float64) where T
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

struct Noise{T}
    # The noise term.
    term::T
    # The probability of sampling the noise term.
    probability::Float64
end

struct State
    # The incoming state variable in the subproblem.
    incoming::JuMP.VariableRef
    # The outgoing state variable in the subproblem.
    outgoing::JuMP.VariableRef
end

mutable struct Node{T}
    # The index of the node in the policy graph.
    index::T
    # The JuMP subproblem.
    subproblem::JuMP.Model
    # A vector of the child nodes.
    children::Vector{Noise{T}}
    # A vector of the discrete stagewise-independent noise terms.
    noise_terms::Vector{Noise}
    # A function parameterize(model::JuMP.Model, noise) that modifies the JuMP
    # model based on the observation of the noise.
    parameterize::Function  # TODO(odow): make this a concrete type?
    # A list of the state variables in the model.
    states::Dict{Symbol, State}
    # Stage objective
    stage_objective  # TODO(odow): make this a concrete type?
    # The optimization sense. Must be :Min or :Max.
    optimization_sense::Symbol
    # Bellman function
    bellman_function  # TODO(odow): make this a concrete type?
end

struct PolicyGraph{T}
    # Children of the root node. child => probability.
    root_children::Vector{Noise{T}}
    # Starting value of the state variables.
    initial_root_state::Dict{Symbol, Float64}
    # All nodes in the graph.
    nodes::Dict{T, Node{T}}
    PolicyGraph(T) = new{T}(Noise{T}[], Dict{Symbol, Float64}(), Dict{T, Node{T}}())
end

# So we can query nodes in the graph as graph[node].
function Base.getindex(graph::PolicyGraph{T}, index::T) where T
    return graph.nodes[index]
end

function get_subproblem(graph::PolicyGraph{T}, index::T) where T
    return graph[index].subproblem::JuMP.Model
end

# Work around different JuMP modes (Automatic / Manual / Direct).
function construct_subproblem(optimizer_factory, direct_mode::Bool)
    subproblem = if direct_mode
        instance = optimizer_factory.constructor(
            optimizer_factory.args...; optimizer_factory.kwargs...)
        JuMP.direct_model(instance)
    else
        JuMP.Model(optimizer_factory)
    end
    return subproblem
end

# Work around different JuMP modes (Automatic / Manual / Direct).
function construct_subproblem(optimizer_factory::Nothing, direct_mode::Bool)
    if direct_mode
        error("You must specify an optimizer in the form:\n" *
              "    with_optimizer(Module.Opimizer, args...) if " *
              "direct_mode=true.")
    end
    return JuMP.Model()
end

"""
    PolicyGraph(builder::Function, graph::Graph{T};
                bellman_function = AverageCut,
                optimizer = nothing,
                direct_mode = true) where T

Construct a a policy graph based on the graph structure of `graph`. (See `Graph`
for details.)

# Example

    function builder(subproblem::JuMP.Model, index)
        # ... subproblem definition ...
    end
    model = PolicyGraph(builder, graph;
                        bellman_function = AverageCut,
                        optimizer = with_optimizer(GLPK.Optimizer),
                        direct_mode = false)

Or, using the Julia `do ... end` syntax:

    model = PolicyGraph(graph;
                        bellman_function = AverageCut,
                        optimizer = with_optimizer(GLPK.Optimizer),
                        direct_mode = true) do subproblem, index
        # ... subproblem definitions ...
    end
"""
function PolicyGraph(builder::Function, graph::Graph{T};
                     bellman_function = AverageCut(),
                     optimizer = nothing,
                     direct_mode = true) where T
    policy_graph = PolicyGraph(T)
    # Initialize nodes.
    for (node_index, children) in graph.nodes
        if node_index == graph.root_node
            continue
        end
        subproblem = construct_subproblem(optimizer, direct_mode)
        node = Node(
            node_index,
            subproblem,
            Noise{T}[],
            Noise[],
            (ω) -> nothing,
            Dict{Symbol, State}(),
            nothing,
            :Min,
            # Delay initializing the bellman function until later so that it can
            # use information about the children and number of
            # stagewise-independent noise realizations.
            nothing

        )
        subproblem.ext[:kokako_policy_graph] = policy_graph
        policy_graph.nodes[node_index] = subproblem.ext[:kokako_node] = node
        builder(subproblem, node_index)
        # Add a dummy noise here so that all nodes have at least one noise term.
        if length(node.noise_terms) == 0
            push!(node.noise_terms, Noise(nothing, 1.0))
        end
    end
    # Loop back through and add the arcs/children.
    for (node_index, children) in graph.nodes
        if node_index == graph.root_node
            continue
        end
        node = policy_graph.nodes[node_index]
        for (child, probability) in children
            push!(node.children, Noise(child, probability))
        end
        # Intialize the bellman function. (See note in creation of Node above.)
        node.bellman_function = initialize_bellman_function(
            bellman_function, policy_graph, node)
    end
    # Add root nodes
    for (child, probability) in graph.nodes[graph.root_node]
        push!(policy_graph.root_children, Noise(child, probability))
    end
    return policy_graph
end

function get_node(subproblem::JuMP.Model)
    return subproblem.ext[:kokako_node]::Node
end
function get_policy_graph(subproblem::JuMP.Model)
    return subproblem.ext[:kokako_policy_graph]::PolicyGraph
end

# Internal function: at the state variable to the subproblem. Requires
# JuMP.name(incoming) and JuMP.start_value(incoming) to be set.
function add_state_variable(subproblem::JuMP.Model,
                            incoming::JuMP.VariableRef,
                            outgoing::JuMP.VariableRef)
    node = get_node(subproblem)
    name = Symbol(JuMP.name(incoming))
    if haskey(node.states, name)
        error("The state $(name) already exists.")
    end
    node.states[name] = State(incoming, outgoing)
    graph = get_policy_graph(subproblem)
    graph.initial_root_state[name] = 0.0  # JuMP.start_value(incoming)
    return
end

# Internal function: an instance of add_state_variable for base-Julia Arrays.
function add_state_variable(subproblem::JuMP.Model,
                            incoming::Array{JuMP.VariableRef},
                            outgoing::Array{JuMP.VariableRef})
    for (incoming_variable, outgoing_variable) in zip(incoming, outgoing)
        add_state_variable(subproblem, incoming_variable, outgoing_variable)
    end
end

# Internal function: an instance of add_state_variable for JuMPArrays.
function add_state_variable(subproblem::JuMP.Model,
                            incoming::JuMPArray{JuMP.VariableRef},
                            outgoing::JuMPArray{JuMP.VariableRef})
    for (incoming_variable, outgoing_variable) in zip(incoming, outgoing)
        add_state_variable(subproblem, incoming_variable, outgoing_variable)
    end
end

# Internal function: an instance of add_state_variable for base-Julia Dicts.
function add_state_variable(subproblem::JuMP.Model,
                            incoming::Dict{T, JuMP.VariableRef},
                            outgoing::Dict{T, JuMP.VariableRef}) where T
    for (incoming_key, incoming_variable) in incoming
        outgoing_variable = outgoing[incoming_key]
        add_state_variable(subproblem, incoming_variable, outgoing_variable)
    end
end

# Internal function: overload for copy since copy(::Symbol) isn't defined.
_copy(x::Symbol) = x
_copy(x::Expr) = copy(x)

# Internal function: given inputs like `1 <= x[i=1:3]` and `0 <= x[i=1:3] <= 1`,
# return the `x[i=1:3]` component.
function get_outgoing_name(outgoing::Expr)
    if outgoing.head == :comparison && length(outgoing.args) == 5
        return _copy(outgoing.args[3])
    elseif (outgoing.head == :call && length(outgoing.args) == 3 &&
            outgoing.args[1] in (:(<=), :(>=), :(==)))
        return _copy(outgoing.args[2])
    else
        error("@state input error: $(outgoing) not valid syntax.")
    end
end
get_outgoing_name(outgoing::Symbol) = outgoing

"""
    @state(subproblem, outgoing_state, incoming_state)

Define a new state variable in the subproblem `subproblem`.

# Examples

One-dimensional state variables:

    @state(subproblem, x′, x == 0.5)

State variables in a JuMPArray:

    @state(subproblem, x′[1:3, [:A, :B]] >= 1, x == 0)

Julia arrays with upper and lower bounds, and the value at the root node depends
upon the index:

    @state(subproblem, 0 <= x′[i=1:3] <= 1, x == i)
"""

macro state(subproblem, outgoing, incoming)
    # subproblem = esc(subproblem)
    if !isa(incoming, Expr) || incoming.head != :call || incoming.args[1] != :(==)
        error("The third argument to @state must be of the form " *
              "[name]==[value at root node].")
    end
    comparison_symbol, in_name, root_value = incoming.args
    # At the moment, we have something like
    #     outgoing = :(0 <= x′[i=1:3] <= 1)
    #     in_name = :(x).
    # We want to create the incoming expression :(x[i=1:3]). To do so, we're
    # going to copy the middle portion of the outgoing expression, and replace
    # the x′ with the incoming name.
    incoming_expression = get_outgoing_name(outgoing)
    if isa(incoming_expression, Expr)
        incoming_expression.args[1] = in_name
    else
        incoming_expression = in_name
    end
    @show incoming_expression, outgoing
    macro bar(args...)
        @show args
        model = esc(args[1])
        state_in = gensym()
        code = quote end
        push!(code.args, :(@variable $model $(Meta.quot(args[2]))))
        @show code
        # push!(code.args, Expr(
            # :(=),
            # $(state_in),
            # Expr(:macrocall, Symbol("@variable"),
                # esc(model)
            #
            # )
        # )
        return code
    end
    return quote
        let
            state_in = @variable($subproblem, $incoming_expression)
            state_out = @variable($subproblem, $outgoing)
            add_state_variable($subproblem, state_in, state_out)
            nothing
        end
    end
end


"""
    parameterize(modify::Function,
                 subproblem::JuMP.Model,
                 realizations::Vector{T},
                 probability::Vector{Float64} = fill(1.0 / length(realizations))
                     ) where T

Add a parameterization function `modify` to `subproblem`. The `modify` function
takes one argument and modifies `subproblem` based on the realization of the
noise sampled from `realizations` with corresponding probabilities
`probability`.

In order to conduct an out-of-sample simulation, `modify` should accept
arguments that are not in realizations (but still of type T).

# Example

    Kokako.parameterize(subproblem, [1, 2, 3], [0.4, 0.3, 0.3]) do ω
        JuMP.set_upper_bound(x, ω)
    end
"""
function parameterize(modify::Function,
                      subproblem::JuMP.Model,
                      realizations::AbstractVector{T},
                      probability::AbstractVector{Float64} = fill(1.0 / length(realizations), length(realizations))
                          ) where T
    node = get_node(subproblem)
    if length(node.noise_terms) != 0
        error("Duplicate calls to Kokako.parameterize detected. Only " *
              "a subproblem at most one time.")
    end
    for (realization, prob) in zip(realizations, probability)
        push!(node.noise_terms, Noise(realization, prob))
    end
    node.parameterize = modify
    return
end

"""
    set_stage_objective(subproblem::JuMP.Model, sense::Symbol,
                        stage_objective)

Set the stage-objective of `subproblem` to `stage_objective` and the
optimization sense to `sense` (which must be `:Min` or `:Max`).

# Example

    Kokako.set_stage_objective(subproblem, :Min, 2x + 1)
"""
function set_stage_objective(subproblem::JuMP.Model, sense::Symbol,
                             stage_objective)
    if !(sense == :Min || sense == :Max)
        error("The optimization sense must be :Min or :Max. It is $(sense).")
    end
    node = get_node(subproblem)
    node.stage_objective = stage_objective
    node.optimization_sense = sense
    return
end

macro stageobjective(subproblem, sense, expr)
    sense = Expr(:quote, sense)
    code = quote
        set_stage_objective(
            $(esc(subproblem)),
            $(esc(sense)),
            $(Expr(:macrocall,
                Symbol("@expression"),
                :LineNumber,
                esc(subproblem),
                esc(expr)
            ))
        )
    end
    return code
end
