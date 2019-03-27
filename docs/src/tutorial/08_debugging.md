# Basic VIII: debugging

Building multistage stochastic programming models is hard. There are a lot of
different pieces that need to be put together, and we typically have no idea of
the optimal policy, so it can be hard (impossible?) to validate the solution.

That said, here are a few tips to verify and validate models built using
`SDDP.jl`.

## Writing subproblems to file


[`SDDP.parameterize`](@ref)

[`SDDP.write_subproblem_to_file`](@ref).

```jldoctest tutorial_eight
using SDDP, GLPK

model = SDDP.LinearPolicyGraph(
            stages = 2,
            lower_bound = 0.0,
            optimizer = with_optimizer(GLPK.Optimizer)
        ) do subproblem, t
    @variable(subproblem, x, SDDP.State, initial_value = 1)
    @variable(subproblem, y)
    @constraint(subproblem, x.in == x.out + y)
    SDDP.parameterize(subproblem, [1.1, 2.2]) do ω
        @stageobjective(subproblem, ω * x.out)
        JuMP.fix(y, ω)
    end
end

# output

A policy graph with 2 nodes.
 Node indices: 1, 2
```

```jldoctest tutorial_eight; filter=[r"┌.+", r"└.+", r"\[...\]"]
julia> SDDP.parameterize(model[1], 1.1)

julia> SDDP.write_subproblem_to_file(model[1], "subproblem_1", format=:lp)
[...]

julia> read("subproblem_1.lp") |> String |> print
minimize
obj: 1.1 x2 + 1 x4
subject to
c1: -1 x3 - 1 x2 + 1 x1 == 0
Bounds
x4 >= 0
x3 == 1.1
```

It is also possible to parameterize the subproblems using values for `ω` that
are not in the original problem formulation.

```jldoctest tutorial_eight; filter=[r"\┌.+", r"\└.+", r"\[...\]"]
julia> SDDP.parameterize(model[1], 3.3)

julia> SDDP.write_subproblem_to_file(model[1], "subproblem_1", format=:lp)
[...]

julia> read("subproblem_1.lp") |> String |> print
minimize
obj: 3.3 x2 + 1 x4
subject to
c1: -1 x3 - 1 x2 + 1 x1 == 0
Bounds
x4 >= 0
x3 == 3.3

julia> rm("subproblem_1.lp")  # Clean up.
```

This concludes or series of basic introductory tutorials for `SDDP.jl`. When
you're ready, continue to our intermediate series of tutorials, beginning with
[Intermediate I: risk](@ref).
