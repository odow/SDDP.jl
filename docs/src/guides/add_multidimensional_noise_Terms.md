# Add multi-dimensional noise terms

```@meta
DocTestSetup = quote
    using SDDP, GLPK
end
```

Multi-dimensional stagewise-independent random variables can be created by
forming the Cartesian product of the random variables.

```jldoctest; filter=[r"\(value = \d, coefficient = \d\)", r"1\-element.+"s]
julia> model = SDDP.LinearPolicyGraph(
               stages=3, lower_bound = 0, optimizer = GLPK.Optimizer
               ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 0.0)
           support = [(value = v, coefficient = c) for v in [1, 2] for c in [3, 4, 5]]
           probability = [pv * pc for pv in [0.5, 0.5] for pc in [0.3, 0.5, 0.2]]
           SDDP.parameterize(subproblem, support, probability) do ω
               JuMP.fix(x.out, ω.value)
               @stageobjective(subproblem, ω.coefficient * x.out)
               println("ω is: ", ω)
           end
       end;

julia> SDDP.simulate(model, 1);
ω is: (value = 1, coefficient = 4)
ω is: (value = 1, coefficient = 3)
ω is: (value = 2, coefficient = 4)
```
