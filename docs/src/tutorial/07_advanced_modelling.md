# Basic VII: modelling tips

This tutorial discusses some different modelling tips.

```@meta
DocTestSetup = quote
    using Kokako, GLPK
end
```

## Multi-dimensional state variables

Just like normal JuMP variables, it is possible to create containers of state
variables.

```jldoctest; filter=r"A policy graph.+"s
julia> model = Kokako.LinearPolicyGraph(
               stages=1, lower_bound = 0, optimizer = with_optimizer(GLPK.Optimizer)
               ) do subproblem, t
           # A scalar state variable.
           @variable(subproblem, x >= 0, Kokako.State, initial_value = 0)
           println("Lower bound of outgoing x is: ", JuMP.lower_bound(x.out))
           # A vector of state variables.
           @variable(subproblem, y[i = 1:2] >= i, Kokako.State, initial_value = i)
           println("Lower bound of outgoing y[1] is: ", JuMP.lower_bound(y[1].out))
           # A JuMP.Containers.DenseAxisArray of state variables.
           @variable(subproblem,
               z[i = 3:4, j = [:A, :B]] >= i, Kokako.State, initial_value = i)
           println("Lower bound of outgoing z[3, :B] is: ", JuMP.lower_bound(z[3, :B].out))
       end;
Lower bound of outgoing x is: 0.0
Lower bound of outgoing y[1] is: 1.0
Lower bound of outgoing z[3, :B] is: 3.0
```

## Multi-dimensional noise terms

Multi-dimensional stagewise-independent random variables can be created by
forming the Cartesian product of the random variables.

```jldoctest; filter=[r"\(value = \d, coefficient = \d\)", r"1\-element.+"s]
julia> model = Kokako.LinearPolicyGraph(
               stages=3, lower_bound = 0, optimizer = with_optimizer(GLPK.Optimizer)
               ) do subproblem, t
           @variable(subproblem, x, Kokako.State, initial_value = 0.0)
           support = [(value = v, coefficient = c) for v in [1, 2] for c in [3, 4, 5]]
           probability = [pv * pc for pv in [0.5, 0.5] for pc in [0.3, 0.5, 0.2]]
           Kokako.parameterize(subproblem, support, probability) do ω
               JuMP.fix(x.out, ω.value)
               @stageobjective(subproblem, ω.coefficient * x.out)
               println("ω is: ", ω)
           end
       end;

julia> Kokako.simulate(model, 1);
ω is: (value = 1, coefficient = 4)
ω is: (value = 1, coefficient = 3)
ω is: (value = 2, coefficient = 4)
```

## Noise in the constraint matrix

`SDDP.jl` supports coefficients in the constraint matrix through the
[`JuMP.set_coefficient`](http://www.juliaopt.org/JuMP.jl/v0.19/constraints/#JuMP.set_coefficient)
function.

```jldoctest; filter=r" \: .+?1.0"
julia> model = Kokako.LinearPolicyGraph(
               stages=3, lower_bound = 0, optimizer = with_optimizer(GLPK.Optimizer)
               ) do subproblem, t
           @variable(subproblem, x, Kokako.State, initial_value = 0.0)
           @constraint(subproblem, emissions, 1x.out <= 1)
           Kokako.parameterize(subproblem, [0.2, 0.5, 1.0]) do ω
               JuMP.set_coefficient(emissions, x.out, ω)
               println(emissions)
           end
           @stageobjective(subproblem, -x.out)
       end
A policy graph with 3 nodes.
 Node indices: 1, 2, 3

julia> Kokako.simulate(model, 1);
emissions : x_out <= 1.0
emissions : 0.2 x_out <= 1.0
emissions : 0.5 x_out <= 1.0
```

!!! note
    JuMP will canonicalize constraints by moving all variables to the left-hand
    side. Thus, `@constraint(model, 0 <= 1 - x.out)` becomes `x.out <= 1`.
    `JuMP.set_coefficient` sets the coefficient on the _canonicalized_
    constraint.

This concludes or series of basic introductory tutorials for `SDDP.jl`. When
you're ready, continue to our intermediate series of tutorials, beginning with
[Intermediate I: risk](@ref).
