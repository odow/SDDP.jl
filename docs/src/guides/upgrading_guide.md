# Upgrade from the old SDDP.jl

`SDDP.jl` under went a major re-write to be compatible with JuMP v0.19 and Julia
v1.0.

Some of the highlights of the new release include

 - Support for "multi-cut"s
 - Support for stagewise-independent noise in the constraint matrix
 - A way to simulate out-of-sample realizations of the stagewise-independent
   noise terms
 - Extensible ways to get information such as dual variables out of a simulation
 - Support for infinite horizon multistage stochastic programs
 - Extensible stopping rules
 - Extensible sampling schemes
 - Improved cut selection routines
 - Better checks for numerical issues
 - A much tidier (and simpler) implementation, with ample commenting throughout
   the code base

## Syntax changes

The complete re-write has resulted in a painful upgrading process as users must
simultaneously upgrade from Julia 0.6 to Julia 1.0, from JuMP 0.18 to JuMP
0.19, and from the old syntax of `SDDP.jl` to the new. In this section, we
outline some of the larger changes. For more information, we recommend reading
the updated tutorials, beginning with [Basic I: first steps](@ref)

#### `SDDPModel`

`SDDPModel` has been replaced in favor of a more general approach to formulating
multistage stochastic optimization problems.

For basic use-cases, `SDDPModel` has been replaced by `LinearPolicyGraph`.

```julia
model = SDDPModel(stages = 3) do subproblem, t
    # subproblem definition.
end

# becomes

model = SDDP.LinearPolicyGraph(stages = 3) do subproblem, t
    # subproblem definition.
end
```

If you used the `markov_transition` feature, `SDDPModel` has been replaced by
`MarkovianPolicyGraph`.

```julia
model = SDDPModel(
      markov_transition = Array{Float64, 2}[
          [ 1.0 ]',
          [ 0.75 0.25 ],
          [ 0.75 0.25 ; 0.25 0.75 ]
      ]) do subproblem, t, i
  # subproblem definition.
end

# becomes

model = SDDP.MarkovianPolicyGraph(
        transition_matrices = Array{Float64, 2}[
            [ 1.0 ]',
            [ 0.75 0.25 ],
            [ 0.75 0.25 ; 0.25 0.75 ]
        ]) do subproblem, node
    t, i = node
    # subproblem definition.
end
```

#### `solver = `

JuMP 0.19 changed the way that solvers are passed to JuMP models.

```julia
SDDPModel(solver = GurobiSolver(OutputFlag = 0))

# becomes

LinearPolicyGraph(optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag = 0))
```

#### `@state`

We changed how to specify state variables.

```julia
@state(subproblem, 0 <= x <= 1, x0==2)

# becomes

@variable(subproblem, 0 <= x <= 1, SDDP.State, initial_value = 2)
```

In addition, instead of having to create an incoming state `x0` and an outgoing
state `x`, we now refer to `x.in` and `x.out`. Here is another example:
```julia
@state(subproblem, 0 <= x[i=1:2] <= 1, x0==i)

# becomes

@variable(subproblem, 0 <= x[i=1:2] <= 1, SDDP.State, initial_value = i)
```

```julia
x0[1]

# becomes

x[1].in
```

```julia
x[2]

# becomes

x[2].out
```

#### `@rhsnoise`

We removed `@rhsnoise`. This results in more lines of code, but more
flexibility.

```julia
@rhsnoise(subproblem, ω = [1, 2, 3], 2x <= ω)
setnoiseprobability!(subproblem, [0.5, 0.2, 0.3])

# becomes

@variable(subproblem, ω)
@constraint(subproblem, 2x <= ω)
SDDP.parameterize(subproblem, [1, 2, 3], [0.5, 0.2, 0.3]) do ϕ
    JuMP.fix(ω, ϕ)
end
```

If you had multiple `@rhsnoise` constraints, use:
```julia
@rhsnoise(subproblem, ω = [1, 2, 3], 2x <= ω)
@rhsnoise(subproblem, ω = [4, 5, 6], y >= 3ω)
setnoiseprobability!(subproblem, [0.5, 0.2, 0.3])

# becomes

@variable(subproblem, ω[1:2])
@constraint(subproblem, 2x <= ω[1])
@constraint(subproblem, y >= 3 * ω[2])
SDDP.parameterize(subproblem, [(1, 4), (2, 5), (3, 6)], [0.5, 0.2, 0.3]) do ϕ
    JuMP.fix(ω[1], ϕ[1])
    JuMP.fix(ω[2], ϕ[2])
end
```

#### `@stageobjective`

`@stageobjective` no longer accepts a random list of parameters. Use
[`SDDP.parameterize`](@ref) instead.
```julia
@stageobjective(subproblem, ω = [1, 2, 3], ω * x)
setnoiseprobability!(subproblem, [0.5, 0.2, 0.3])

# becomes

SDDP.parameterize(subproblem, [1, 2, 3], [0.5, 0.2, 0.3]) do ω
    @stageobjective(subproblem, ω * x)
end
```

#### `SDDP.solve`

`SDDP.solve` has been replaced by [`SDDP.train`](@ref). See the docs for a
complete list of the new options as most things have changed.

!!! note
    Parallel training has _not_ (yet) been implemented.

#### Plotting

Much of the syntax for plotting has changed. See [Basic V: plotting](@ref) for
the new syntax.

#### Price interpolation

The syntax for models with stagewise-dependent objective processes has
completely changed. See [Advanced I: objective states](@ref) for details.
