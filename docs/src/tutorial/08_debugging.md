# Basic VIII: debugging

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
            optimizer = with_optimizer(GLPK.Optimizer),
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
julia> SDDP.write_subproblem_to_file(model[1], "subproblem", format=:lp)
┌ Warning: Blank name detected for variable MathOptInterface.VariableIndex(4). Renamed to x4.
└ @ MathOptFormat ~/.julia/packages/MathOptFormat/iRtuE/src/MathOptFormat.jl:95
┌ Warning: Blank name detected for constraint MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.EqualTo{Float64}}(1). Renamed to c1.
└ @ MathOptFormat ~/.julia/packages/MathOptFormat/iRtuE/src/MathOptFormat.jl:54
┌ Warning: Blank name detected for constraint MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}}(2). Renamed to c2.
└ @ MathOptFormat ~/.julia/packages/MathOptFormat/iRtuE/src/MathOptFormat.jl:54
```

We can check the file by reading it back in again.

```jldoctest tutorial_eight
julia> read("subproblem.lp") |> String |> print
minimize
obj: 1.1 x_out + 1 x4
subject to
balance: 1 x_in - 1 x_out - 1 y == 0
Bounds
x4 >= 0
y == 1.1
```

It is easy to see that `ω` has been set in the objective, and as the fixed value
for `y`.

It is also possible to parameterize the subproblems using values for `ω` that
are not in the original problem formulation.

```jldoctest tutorial_eight; filter=r"MathOptFormat\ .+?MathOptFormat\.jl"
julia> SDDP.parameterize(model[1], 3.3)

julia> SDDP.write_subproblem_to_file(model[1], "subproblem", format=:lp)
┌ Warning: Blank name detected for variable MathOptInterface.VariableIndex(4). Renamed to x4.
└ @ MathOptFormat ~/.julia/packages/MathOptFormat/iRtuE/src/MathOptFormat.jl:95
┌ Warning: Blank name detected for constraint MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.EqualTo{Float64}}(1). Renamed to c1.
└ @ MathOptFormat ~/.julia/packages/MathOptFormat/iRtuE/src/MathOptFormat.jl:54
┌ Warning: Blank name detected for constraint MathOptInterface.ConstraintIndex{MathOptInterface.SingleVariable,MathOptInterface.GreaterThan{Float64}}(2). Renamed to c2.
└ @ MathOptFormat ~/.julia/packages/MathOptFormat/iRtuE/src/MathOptFormat.jl:54

julia> read("subproblem.lp") |> String |> print
minimize
obj: 3.3 x_out + 1 x4
subject to
balance: 1 x_in - 1 x_out - 1 y == 0
Bounds
x4 >= 0
y == 3.3

julia> rm("subproblem.lp")  # Clean up.
```

This concludes or series of basic introductory tutorials for `SDDP.jl`. When
you're ready, continue to our intermediate series of tutorials, beginning with
[Intermediate I: risk](@ref).
