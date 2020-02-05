# Add a multi-dimensional state variable

```@meta
DocTestSetup = quote
    using SDDP, GLPK
end
```

Just like normal JuMP variables, it is possible to create containers of state
variables.

```jldoctest; filter=r"A policy graph.+"s
julia> model = SDDP.LinearPolicyGraph(
           stages=1, lower_bound = 0, optimizer = GLPK.Optimizer
       ) do subproblem, t
           # A scalar state variable.
           @variable(subproblem, x >= 0, SDDP.State, initial_value = 0)
           println("Lower bound of outgoing x is: ", JuMP.lower_bound(x.out))
           # A vector of state variables.
           @variable(subproblem, y[i = 1:2] >= i, SDDP.State, initial_value = i)
           println("Lower bound of outgoing y[1] is: ", JuMP.lower_bound(y[1].out))
           # A JuMP.Containers.DenseAxisArray of state variables.
           @variable(subproblem,
               z[i = 3:4, j = [:A, :B]] >= i, SDDP.State, initial_value = i)
           println("Lower bound of outgoing z[3, :B] is: ", JuMP.lower_bound(z[3, :B].out))
       end;
Lower bound of outgoing x is: 0.0
Lower bound of outgoing y[1] is: 1.0
Lower bound of outgoing z[3, :B] is: 3.0
```
