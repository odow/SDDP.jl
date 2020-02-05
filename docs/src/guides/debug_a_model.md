# Debug a model

Building multistage stochastic programming models is hard. There are a lot of
different pieces that need to be put together, and we typically have no idea of
the optimal policy, so it can be hard (impossible?) to validate the solution.

That said, here are a few tips to verify and validate models built using
`SDDP.jl`.

## Writing subproblems to file

The first step to debug a model is to write out the subproblems to a file in
order to check that you are actually building what you think you are building.

This can be achieved with the help of two functions: [`SDDP.parameterize`](@ref)
and [`SDDP.write_subproblem_to_file`](@ref). The first lets you parameterize a
node given a noise, and the second writes out the subproblem to a file.

Here is an example model:

```jldoctest tutorial_eight
using SDDP, GLPK

model = SDDP.LinearPolicyGraph(
            stages = 2,
            lower_bound = 0.0,
            optimizer = GLPK.Optimizer,
            direct_mode = false
        ) do subproblem, t
    @variable(subproblem, x, SDDP.State, initial_value = 1)
    @variable(subproblem, y)
    @constraint(subproblem, balance, x.in == x.out + y)
    SDDP.parameterize(subproblem, [1.1, 2.2]) do ω
        @stageobjective(subproblem, ω * x.out)
        JuMP.fix(y, ω)
    end
end

# output

A policy graph with 2 nodes.
 Node indices: 1, 2
```

Initially, `model` hasn't been parameterized with a concrete realizations of
`ω`. Let's do so now by parameterizing the first subproblem with `ω=1.1`.
```jldoctest tutorial_eight
julia> SDDP.parameterize(model[1], 1.1)
```
Easy! To parameterize the second stage problem, we would have used `model[2]`.

Now to write out the problem to a file. We'll get a few warnings because some
variables and constraints don't have names. They don't matter, so ignore them.

```jldoctest tutorial_eight; filter=r"MathOptFormat\ .+?MathOptFormat\.jl"
julia> SDDP.write_subproblem_to_file(model[1], "subproblem.lp")

julia> read("subproblem.lp") |> String |> print
minimize
obj: 1.1 x_out + 1 x2
subject to
balance: 1 x_in - 1 x_out - 1 y = 0
Bounds
x2 >= 0
y = 1.1
x_in free
x_out free
End
```

It is easy to see that `ω` has been set in the objective, and as the fixed value
for `y`.

It is also possible to parameterize the subproblems using values for `ω` that
are not in the original problem formulation.

```jldoctest tutorial_eight; filter=r"MathOptFormat\ .+?MathOptFormat\.jl"
julia> SDDP.parameterize(model[1], 3.3)

julia> SDDP.write_subproblem_to_file(model[1], "subproblem.lp")

julia> read("subproblem.lp") |> String |> print
minimize
obj: 3.3 x_out + 1 x2
subject to
balance: 1 x_in - 1 x_out - 1 y = 0
Bounds
x2 >= 0
y = 3.3
x_in free
x_out free
End

julia> rm("subproblem.lp")  # Clean up.
```

## Solve the determinstic equivalent

Sometimes, it can be helpful to solve the deterministic equivalent of a
problem in order to obtain an exact solution to the problem. To obtain a JuMP
model that represents the deterministic equivalent, use [`SDDP.deterministic_equivalent`](@ref).
The returned model is just a normal JuMP model. Use JuMP to optimize it and
query the solution.

```jldoctest tutorial_eight; filter=r"5.4725[0]+[0-9]"
julia> det_equiv = SDDP.deterministic_equivalent(model, GLPK.Optimizer)
A JuMP Model
Feasibility problem with:
Variables: 24
`GenericAffExpr{Float64,VariableRef}`-in-`MathOptInterface.EqualTo{Float64}`: 10 constraints
`VariableRef`-in-`MathOptInterface.EqualTo{Float64}`: 8 constraints
`VariableRef`-in-`MathOptInterface.GreaterThan{Float64}`: 6 constraints
`VariableRef`-in-`MathOptInterface.LessThan{Float64}`: 4 constraints
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: GLPK

julia> optimize!(det_equiv)

julia> objective_value(det_equiv)
-5.472500000000001
```

!!! warning
    The determinstic equivalent scales poorly with problem size. Only use this
    on small problems!
