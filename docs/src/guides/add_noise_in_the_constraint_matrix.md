# Add noise in the constraint matrix

```@meta
DocTestSetup = quote
    using SDDP, GLPK
end
```

`SDDP.jl` supports coefficients in the constraint matrix through the
[`JuMP.set_normalized_coefficient`](http://www.juliaopt.org/JuMP.jl/v0.20/constraints/#JuMP.set_normalized_coefficient)
function.

```jldoctest; filter=r" \: .+?1.0"
julia> model = SDDP.LinearPolicyGraph(
               stages=3, lower_bound = 0, optimizer = GLPK.Optimizer
               ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 0.0)
           @constraint(subproblem, emissions, 1x.out <= 1)
           SDDP.parameterize(subproblem, [0.2, 0.5, 1.0]) do ω
               JuMP.set_normalized_coefficient(emissions, x.out, ω)
               println(emissions)
           end
           @stageobjective(subproblem, -x.out)
       end
A policy graph with 3 nodes.
 Node indices: 1, 2, 3

julia> SDDP.simulate(model, 1);
emissions : x_out <= 1.0
emissions : 0.2 x_out <= 1.0
emissions : 0.5 x_out <= 1.0
```

!!! note
    JuMP will normalize constraints by moving all variables to the left-hand
    side. Thus, `@constraint(model, 0 <= 1 - x.out)` becomes `x.out <= 1`.
    `JuMP.set_normalized_coefficient` sets the coefficient on the _normalized_
    constraint.
