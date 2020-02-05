# Create a general policy graph

```@meta
DocTestSetup = quote
    using SDDP, GLPK
end
```

SDDP.jl uses the concept of a _policy graph_ to formulate multistage stochastic
programming problems. We _highly_ recommend that you read the following paper
before continuing with this tutorial.

 - Dowson, O. (2018). The policy graph decomposition of multistage stochastic
   optimization problems. Optimization Online. [link](http://www.optimization-online.org/DB_HTML/2018/11/6914.html)

## Creating a [`SDDP.Graph`](@ref)

### Linear graphs

Linear policy graphs can be created using the [`SDDP.LinearGraph`](@ref)
function.

```jldoctest linear_graph
julia> graph = SDDP.LinearGraph(3)
Root
 0
Nodes
 1
 2
 3
Arcs
 0 => 1 w.p. 1.0
 1 => 2 w.p. 1.0
 2 => 3 w.p. 1.0
```

We can add nodes to a graph using [`SDDP.add_node`](@ref) and edges using
[`SDDP.add_edge`](@ref).

```jldoctest linear_graph
julia> SDDP.add_node(graph, 4)

julia> SDDP.add_edge(graph, 3 => 4, 1.0)

julia> SDDP.add_edge(graph, 4 => 1, 0.9)

julia> graph
Root
 0
Nodes
 1
 2
 3
 4
Arcs
 0 => 1 w.p. 1.0
 1 => 2 w.p. 1.0
 2 => 3 w.p. 1.0
 3 => 4 w.p. 1.0
 4 => 1 w.p. 0.9
```

Look! We just made a cyclic graph! SDDP.jl can solve infinite horizon problems.
The probability on the arc that completes a cycle should be interpreted as a
discount factor.

### Markovian policy graphs

Markovian policy graphs can be created using the [`SDDP.MarkovianGraph`](@ref)
function.

```jldoctest
julia> SDDP.MarkovianGraph(Matrix{Float64}[[1.0]', [0.4 0.6]])
Root
 (0, 1)
Nodes
 (1, 1)
 (2, 1)
 (2, 2)
Arcs
 (0, 1) => (1, 1) w.p. 1.0
 (1, 1) => (2, 1) w.p. 0.4
 (1, 1) => (2, 2) w.p. 0.6
```

### General graphs

Arbitrarily complicated graphs can be constructed using [`SDDP.Graph`](@ref),
[`SDDP.add_node`](@ref) and [`SDDP.add_edge`](@ref). For example

```jldoctest
julia> graph = SDDP.Graph(:root_node)
Root
 root_node
Nodes
Arcs

julia> SDDP.add_node(graph, :decision_node)

julia> SDDP.add_edge(graph, :root_node => :decision_node, 1.0)

julia> SDDP.add_edge(graph, :decision_node => :decision_node, 0.9)

julia> graph
Root
 root_node
Nodes
 decision_node
Arcs
 root_node => decision_node w.p. 1.0
 decision_node => decision_node w.p. 0.9
```

## Creating a policy graph

Once you have constructed an instance of [`SDDP.Graph`], you can create a
policy graph by passing the graph as the first argument.

```jldoctest
julia> graph = SDDP.Graph(
           :root_node,
           [:decision_node],
           [
               (:root_node => :decision_node, 1.0),
               (:decision_node => :decision_node, 0.9)
           ]);

julia> model = SDDP.PolicyGraph(
               graph,
               lower_bound = 0,
               optimizer = GLPK.Optimizer) do subproblem, node
           println("Called from node: ", node)
       end;
Called from node: decision_node
```

### Special cases

There are two special cases which cover the majority of models in the
literature.

- [`SDDP.LinearPolicyGraph`](@ref) is a special case where a
  [`SDDP.LinearGraph`](@ref) is passed as the first argument.

- [`SDDP.MarkovianPolicyGraph`](@ref) is a special case where a
  [`SDDP.MarkovianGraph`](@ref) is passed as the first argument.

Note that the type of the names of all nodes (including the root node) must be
the same. In this case, they are `Symbol`s.

## Simulating non-standard policy graphs

If you simulate a policy graph with a node that has outgoing arcs that sum to less than one,
you will end up with simulations of different lengths. (The most common case is an infinite
horizon stochastic program, aka a linear policy graph with a single cycle.)

To simulate a fixed number of stages, use:
```julia
simulations = SDDP.simulate(
    model,
    1,
    sampling_scheme = SDDP.InSampleMonteCarlo(
        max_depth = 10,
        terminate_on_dummy_leaf = false
    )
)
```
Here, `max_depth` controls the number of stages, and `terminate_on_dummy_leaf = false` stops
us from terminating early.

See also [Simulate using a different sampling scheme](@ref).

## Creating a Markovian graph automatically

SDDP.jl can create a Markovian graph by automatically discretizing a one-dimensional
stochastic process and fitting a Markov chain.

To access this functionality, pass a function that takes no arguments and returns a
`Vector{Float64}` to [`SDDP.MarkovianGraph`](@ref). To keyword arguments also need to be
provided: `budget` is the total number of nodes in the Markovian graph, and `scenarios` is
the number of realizations of the simulator function used to approximate the graph.

In some cases, `scenarios` may be too small to provide a reasonable fit of the stochastic
process. If so, SDDP.jl will automatically try to re-fit the Markov chain using more
scenarios.

```julia
function simulator()
    scenario = zeros(5)
    for i = 2:5
        scenario[i] = scenario[i - 1] + rand() - 0.5
    end
    return scenario
end

model = SDDP.PolicyGraph(
    SDDP.MarkovianGraph(simulator; budget = 10, scenarios = 100),
    sense = :Max,
    upper_bound = 1e3
) do subproblem, node
    (stage, price) = node
    @variable(subproblem, x >= 0, SDDP.State, initial_value = 1)
    @constraint(subproblem, x.out <= x.in)
    @stageobjective(subproblem, price * x.out)
end
```
