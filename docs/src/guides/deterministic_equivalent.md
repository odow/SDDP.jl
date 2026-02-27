# How to solve the deterministic equivalent

Sometimes, it can be helpful to solve the deterministic equivalent of a
problem in order to obtain an exact solution to the problem. To obtain a JuMP
model that represents the deterministic equivalent, use [`SDDP.deterministic_equivalent`](@ref).
The returned model is just a normal JuMP model. Use JuMP to optimize it and
query the solution.

```jldoctest; filter=r"5.4725[0]+[0-9]"
julia> using SDDP, HiGHS

julia> model = SDDP.LinearPolicyGraph(
                   stages = 2,
                   lower_bound = 0.0,
                   optimizer = HiGHS.Optimizer,
               ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 1)
           @variable(subproblem, y)
           @constraint(subproblem, balance, x.in == x.out + y)
           SDDP.parameterize(subproblem, [1.1, 2.2]) do ω
               @stageobjective(subproblem, ω * x.out)
               fix(y, ω)
           end
       end
A policy graph with 2 nodes.
 Node indices: 1, 2

julia> det_equiv = SDDP.deterministic_equivalent(model, HiGHS.Optimizer)
A JuMP Model
├ solver: HiGHS
├ objective_sense: MIN_SENSE
│ └ objective_function_type: AffExpr
├ num_variables: 24
├ num_constraints: 28
│ ├ AffExpr in MOI.EqualTo{Float64}: 10
│ ├ VariableRef in MOI.EqualTo{Float64}: 8
│ ├ VariableRef in MOI.GreaterThan{Float64}: 6
│ └ VariableRef in MOI.LessThan{Float64}: 4
└ Names registered in the model: none

julia> set_silent(det_equiv)

julia> optimize!(det_equiv)

julia> objective_value(det_equiv)
-5.472500000000001
```

!!! warning
    The deterministic equivalent scales poorly with problem size. Only use this
    on small problems!
