```@meta
CurrentModule = SDDP
```

# Intermediate V: belief states

`SDDP.jl` includes an implementation of the algorithm described in Dowson, O.,
Morton, D.P., & Pagnoncelli, B. (2019). Partially observable multistage
stochastic programming. [[link]](http://www.optimization-online.org/DB_HTML/2019/03/7141.html)

Proper documentation will be forthcoming.

In the mean-time, here is most of what you need to know.

## Defining the ambiguity partition

Given a [`SDDP.Graph`](@ref) object (see [Intermediate III: policy graphs](@ref)
for details), we can define the ambiguity partition using [`SDDP.add_ambiguity_set`](@ref).

For example:
```jldoctest; setup=:(using SDDP)
julia> G = SDDP.MarkovianGraph([[0.5 0.5], [0.2 0.8; 0.8 0.2]])
Root
 (0, 1)
Nodes
 (1, 1)
 (1, 2)
 (2, 1)
 (2, 2)
Arcs
 (0, 1) => (1, 1) w.p. 0.5
 (0, 1) => (1, 2) w.p. 0.5
 (1, 1) => (2, 1) w.p. 0.2
 (1, 1) => (2, 2) w.p. 0.8
 (1, 2) => (2, 1) w.p. 0.8
 (1, 2) => (2, 2) w.p. 0.2

julia> SDDP.add_ambiguity_set(G, [(1, 1), (1, 2)])

julia> SDDP.add_ambiguity_set(G, [(2, 1), (2, 2)])

julia> G
Root
 (0, 1)
Nodes
 (1, 1)
 (1, 2)
 (2, 1)
 (2, 2)
Arcs
 (0, 1) => (1, 1) w.p. 0.5
 (0, 1) => (1, 2) w.p. 0.5
 (1, 1) => (2, 1) w.p. 0.2
 (1, 1) => (2, 2) w.p. 0.8
 (1, 2) => (2, 1) w.p. 0.8
 (1, 2) => (2, 2) w.p. 0.2
Partition
 {
    (1, 1)
    (1, 2)
 } {
    (2, 1)
    (2, 2)
 }
```

In the next tutorial, [Intermediate VI: performance](@ref), we discuss how to
improve the computational performance of `SDDP.jl` models.
