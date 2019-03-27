#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

struct Graph{T}
    # The root node of the policy graph.
    root_node::T
    # nodes[x] returns a vector of the children of node x and their
    # probabilities.
    nodes::Dict{T, Vector{Tuple{T, Float64}}}
end

# Helper utilities to sort the nodes for printing. This helps linear and
# Markovian policy graphs where the nodes might be stored in an unusual ordering
# in the dictionary.
sort_nodes(nodes::Vector{Int}) = sort!(nodes)
sort_nodes(nodes::Vector{Tuple{Int, Int}}) = sort!(nodes)
sort_nodes(nodes) = nodes

function Base.show(io::IO, graph::Graph)
    println(io, "Root")
    println(io, " ", graph.root_node)
    println(io, "Nodes")
    nodes = sort_nodes(collect(keys(graph.nodes)))
    println.(Ref(io), " ", filter(n -> n != graph.root_node, nodes))
    println(io, "Arcs")
    for node in nodes
        for (child, probability) in graph.nodes[node]
            println(io, " ", node, " => ", child, " w.p. ", probability)
        end
    end
end

"""
    Graph(root_node::T) where T

Create an empty graph struture with the root node `root_node`.
"""
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
            error("Rows in the transition matrix must sum to between 0.0 and " *
                  "1.0.")
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
        vcat([Base.reshape(root_node_transition, 1, length(root_node_transition))],
             [transition_matrix for stage in 1:(stages - 1)])
    )
end

struct Noise{T}
    # The noise term.
    term::T
    # The probability of sampling the noise term.
    probability::Float64
end

struct State{T}
    # The incoming state variable.
    in::T
    # The outgoing state variable.
    out::T
end

mutable struct ObjectiveState{N}
    update::Function
    initial_value::NTuple{N, Float64}
    state::NTuple{N, Float64}
    lower_bound::NTuple{N, Float64}
    upper_bound::NTuple{N, Float64}
    μ::NTuple{N, JuMP.VariableRef}
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
    states::Dict{Symbol, State{JuMP.VariableRef}}
    # Stage objective
    stage_objective  # TODO(odow): make this a concrete type?
    stage_objective_set::Bool
    # Bellman function
    bellman_function  # TODO(odow): make this a concrete type?
    # For dynamic interpolation of objective states.
    objective_state::Union{Nothing, ObjectiveState}
end

mutable struct PolicyGraph{T}
    # Must be MOI.MIN_SENSE or MOI.MAX_SENSE
    objective_sense::MOI.OptimizationSense
    # Children of the root node. child => probability.
    root_children::Vector{Noise{T}}
    # Starting value of the state variables.
    initial_root_state::Dict{Symbol, Float64}
    # All nodes in the graph.
    nodes::Dict{T, Node{T}}
    # Storage for the most recent training results.
    most_recent_training_results
    function PolicyGraph(T, sense::Symbol)
        optimization_sense = if sense == :Min
            MOI.MIN_SENSE
        elseif sense == :Max
            MOI.MAX_SENSE
        else
            error("The optimization sense must be :Min or :Max. It is $(sense).")
        end
        new{T}(optimization_sense, Noise{T}[], Dict{Symbol, Float64}(), Dict{T, Node{T}}(), nothing)
    end
end

function Base.show(io::IO, graph::PolicyGraph)
    println(io, "A policy graph with $(length(graph.nodes)) nodes.")
    println(io, " Node indices: ",
        join(sort_nodes(collect(keys(graph.nodes))), ", "))
end

# So we can query nodes in the graph as graph[node].
function Base.getindex(graph::PolicyGraph{T}, index::T) where T
    return graph.nodes[index]
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
    LinearPolicyGraph(builder::Function; stages::Int, kwargs...)

Create a linear policy graph with `stages` number of stages.

See [`PolicyGraph`](@ref) for the other keyword arguments.
"""
function LinearPolicyGraph(builder::Function; stages::Int, kwargs...)
    if stages < 1
        error("You must create a LinearPolicyGraph with `stages >= 1`.")
    end
    return PolicyGraph(builder, LinearGraph(stages); kwargs...)
end

"""
    MarkovianPolicyGraph(builder::Function;
        transition_matrices::Vector{Array{Float64, 2}}, kwargs...)

Create a Markovian policy graph based on the transition matrices given in
`transition_matrices`.

See [`PolicyGraph`](@ref) for the other keyword arguments.
"""
function MarkovianPolicyGraph(builder::Function;
        transition_matrices::Vector{Array{Float64, 2}}, kwargs...)
    return PolicyGraph(builder, MarkovianGraph(transition_matrices); kwargs...)
end

"""
    PolicyGraph(builder::Function, graph::Graph{T};
                bellman_function = BellmanFunction,
                optimizer = nothing,
                direct_mode = true) where T

Construct a a policy graph based on the graph structure of `graph`. (See `Graph`
for details.)

# Example

    function builder(subproblem::JuMP.Model, index)
        # ... subproblem definition ...
    end
    model = PolicyGraph(builder, graph;
                        bellman_function = BellmanFunction,
                        optimizer = with_optimizer(GLPK.Optimizer),
                        direct_mode = false)

Or, using the Julia `do ... end` syntax:

    model = PolicyGraph(graph;
                        bellman_function = BellmanFunction,
                        optimizer = with_optimizer(GLPK.Optimizer),
                        direct_mode = true) do subproblem, index
        # ... subproblem definitions ...
    end
"""
function PolicyGraph(builder::Function, graph::Graph{T};
                     sense = :Min,
                     bellman_function = nothing,
                     lower_bound = -Inf,
                     upper_bound = Inf,
                     optimizer = nothing,
                     direct_mode = true) where T
    policy_graph = PolicyGraph(T, sense)

    # Create a Bellman function if one is not given.
    if bellman_function === nothing
        if lower_bound === -Inf && upper_bound === Inf
            error("You must specify a bound on the objective value, through " *
                  "`lower_bound` if minimizing, or `upper_bound` if maximizing.")
        else
            bellman_function = BellmanFunction(
                lower_bound = lower_bound, upper_bound = upper_bound)
        end
    end

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
            Dict{Symbol, State{JuMP.VariableRef}}(),
            nothing,
            false,
            # Delay initializing the bellman function until later so that it can
            # use information about the children and number of
            # stagewise-independent noise realizations.
            nothing,
            # Likewise for the objective states.
            nothing

        )
        subproblem.ext[:kokako_policy_graph] = policy_graph
        policy_graph.nodes[node_index] = subproblem.ext[:kokako_node] = node
        JuMP.set_objective_sense(subproblem, policy_graph.objective_sense)
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

# Internal functino: helper to get the node given a subproblem.
function get_node(subproblem::JuMP.Model)
    return subproblem.ext[:kokako_node]::Node
end

# Internal functino: helper to get the policy graph given a subproblem.
function get_policy_graph(subproblem::JuMP.Model)
    return subproblem.ext[:kokako_policy_graph]::PolicyGraph
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

    SDDP.parameterize(subproblem, [1, 2, 3], [0.4, 0.3, 0.3]) do ω
        JuMP.set_upper_bound(x, ω)
    end
"""
function parameterize(modify::Function,
                      subproblem::JuMP.Model,
                      realizations::AbstractVector{T},
                      probability::AbstractVector{Float64} =
                          fill(1.0 / length(realizations), length(realizations))
                          ) where T
    node = get_node(subproblem)
    if length(node.noise_terms) != 0
        error("Duplicate calls to SDDP.parameterize detected.")
    end
    for (realization, prob) in zip(realizations, probability)
        push!(node.noise_terms, Noise(realization, prob))
    end
    node.parameterize = modify
    return
end

"""
    set_stage_objective(subproblem::JuMP.Model, stage_objective)

Set the stage-objective of `subproblem` to `stage_objective`.

# Example

    SDDP.set_stage_objective(subproblem, 2x + 1)
"""
function set_stage_objective(subproblem::JuMP.Model, stage_objective)
    node = get_node(subproblem)
    node.stage_objective = stage_objective
    node.stage_objective_set = false
    return
end

"""
    @stageobjective(subproblem, expr)

Set the stage-objective of `subproblem` to `expr`.

### Example

    @stageobjective(subproblem, 2x + y)
"""
macro stageobjective(subproblem, expr)
    code = quote
        set_stage_objective(
            $(esc(subproblem)),
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

# ============================================================================ #
#
#   Code to implement a JuMP variable extension.
#
#   Usage:
#   julia> @variable(subproblem, 0 <= x[i=1:2] <= i,
#              SDDP.State, initial_value = i)
#
#   julia> x
#   2-element Array{State{VariableRef},1}:
#     State(x[1]_in,x[1]_out)
#     State(x[2]_in,x[2]_out)
#
#   julia> x[1].in
#   x[1]_in
#
#   julia> typeof(x[1].in)
#   VariableRef
#
#   julia> x[2].out
#   x[2]_out
#
#   Assuming subproblem has been solved, and there exists a primal solution
#   julia> x_values = JuMP.value.(x)
#   2-element Array{State{Float64},1}:
#     State(0.0,1.0)
#     State(1.2,3.0)
#
#   julia> x_values[1].out
#   1.0
# ============================================================================ #

struct StateInfo
    in::JuMP.VariableInfo
    out::JuMP.VariableInfo
    initial_value::Float64
end

function JuMP.build_variable(
        _error::Function, info::JuMP.VariableInfo, ::Type{State};
        initial_value = NaN,
        kwargs...)
    if isnan(initial_value)
        _error("When creating a state variable, you must set the " *
               "`initial_value` keyword to the value of the state variable at" *
               " the root node.")
    end
    return StateInfo(
        JuMP.VariableInfo(
            false, NaN,  # lower bound
            false, NaN,  # upper bound
            false, NaN,  # fixed value
            false, NaN,  # start value
            false, false # binary and integer
        ),
        info,
        initial_value
    )
end

function JuMP.add_variable(
        subproblem::JuMP.Model, state_info::StateInfo, name::String)
    state = State(
        JuMP.add_variable(
            subproblem, JuMP.ScalarVariable(state_info.in), name * "_in"),
        JuMP.add_variable(
            subproblem, JuMP.ScalarVariable(state_info.out), name * "_out")
    )
    node = get_node(subproblem)
    sym_name = Symbol(name)
    @assert !haskey(node.states, sym_name)  # JuMP prevents duplicate names.
    node.states[sym_name] = state
    graph = get_policy_graph(subproblem)
    graph.initial_root_state[sym_name] = state_info.initial_value
    return state
end

JuMP.variable_type(model::JuMP.Model, ::Type{State}) = State

function JuMP.value(state::State{JuMP.VariableRef})
    return State(JuMP.value(state.in), JuMP.value(state.out))
end

Broadcast.broadcastable(x::State{T}) where {T} = Ref{State{T}}(x)

# ==============================================================================

"""
    add_objective_state(update::Function, subproblem::JuMP.Model; kwargs...)

Add an objective state variable to `subproblem`.

Required `kwargs` are:

 - `initial_value`: The initial value of the objective state variable at the
    root node.
 - `lipschitz`: The lipschitz constant of the objective state variable.

Setting a tight value for the lipschitz constant can significantly improve the
speed of convergence.

Optional `kwargs` are:

 - `lower_bound`: A valid lower bound for the objective state variable. Can be
    `-Inf`.
 - `upper_bound`: A valid upper bound for the objective state variable. Can be
    `+Inf`.

Setting tight values for these optional variables can significantly improve the
speed of convergence.

If the objective state is `N`-dimensional, each keyword argument must be an
`NTuple{N, Float64}`. For example, `initial_value = (0.0, 1.0)`.
"""
function add_objective_state(
        update::Function, subproblem::JuMP.Model;
        initial_value, lipschitz, lower_bound = -Inf, upper_bound = Inf)
    return add_objective_state(update, subproblem, initial_value, lower_bound,
        upper_bound, lipschitz)
end

# Internal function: add_objective_state with positional Float64 arguments.
function add_objective_state(
        update::Function, subproblem::JuMP.Model,
        initial_value::Float64, lower_bound::Float64,
        upper_bound::Float64, lipschitz::Float64)
    return add_objective_state(update, subproblem, (initial_value,),
        (lower_bound,), (upper_bound,), (lipschitz,))
end

# Internal function: add_objective_state with positional NTuple arguments.
function add_objective_state(
        update::Function, subproblem::JuMP.Model,
        initial_value::NTuple{N, Float64}, lower_bound::NTuple{N, Float64},
        upper_bound::NTuple{N, Float64}, lipschitz::NTuple{N, Float64}) where {N}
    node = get_node(subproblem)
    if node.objective_state !== nothing
        error("add_objective_state can only be called once.")
    end
    μ = @variable(
        subproblem, [i = 1:N],
        lower_bound = -lipschitz[i], upper_bound = lipschitz[i])
    node.objective_state = ObjectiveState(
        update, initial_value, initial_value, lower_bound, upper_bound,
        tuple(μ...))
    return
end

"""
    objective_state(subproblem::JuMP.Model)
"""
function objective_state(subproblem::JuMP.Model)
    objective_state = get_node(subproblem).objective_state
    if objective_state !== nothing
        if length(objective_state.state) == 1
            return objective_state.state[1]
        else
            return objective_state.state
        end
    else
        error("No objective state defined.")
    end
end

# Internal function: calculate <y, μ>.
function get_objective_state_component(node::Node)
    objective_state_component = JuMP.AffExpr(0.0)
    objective_state = node.objective_state
    if objective_state !== nothing
        for (y, μ) in zip(objective_state.state, objective_state.μ)
            objective_state_component += y * μ
        end
    end
    return objective_state_component
end
