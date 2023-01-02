```@meta
DocTestSetup = quote
    using SDDP
end
```

# Create a belief state

`SDDP.jl` includes an implementation of the algorithm described in Dowson, O.,
Morton, D.P., & Pagnoncelli, B.K. (2020). Partially observable multistage
stochastic optimization. _Operations Research Letters_, 48(4), 505--512.

Given a [`SDDP.Graph`](@ref) object (see [Create a general policy graph](@ref)
for details), we can define the ambiguity partition using
[`SDDP.add_ambiguity_set`](@ref).

For example, first we create a Markovian graph:

```@repl
using SDDP
G = SDDP.MarkovianGraph([[0.5 0.5], [0.2 0.8; 0.8 0.2]])
```

Then we add an ambiguity set over the nodes in the each stage:

```@repl
for t in 1:2
    SDDP.add_ambiguity_set(G, [(t, 1), (t, 2)])
end
```

This results in the graph:

```@repl
G
```
