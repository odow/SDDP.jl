# SDDP
## Installation
This package is unregistered so you will need to `Pkg.clone` it as follow:
```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```

## Usage
See the `/examples` folder for example usage but briefly

### Initialising the model object
The first step is to initialise the SDDP model object. We do this using the following syntax:

If we have more than one markov state:
```julia
m = SDDPModel([;kwargs...]) do sp, stage, markov_state
  # Stage problem definition where `sp` is a `JuMP.Model`object,
end
```


Otherwise if we have a single markov state
```julia
m = SDDPModel([;kwargs...]) do sp, stage
  # Stage problem definition
end
```

`SDDPModel` takes the following keyword arguments.
+ `sense::Symbol` = `:Min` or `:Max`
+ `stages::Int`
+ `objective_bound`
+ `markov_transition`
+ `scenario_probability`
+ `risk_measure`
+ `cut_oracle`
+ `solver::MathProgBase.AbstractMathProgSolver`

### Describing the Stage Problems
We now need to define the stage problems.

#### State Variables
We can define a new state variable in the stage problem `sp` using the `@state` macro:

It consists of three arguments

1. the stage problem model object.
2. the value of the state variable at the END of the stage. It can consist of any valid JuMP `@variable` syntax.
3. the value of the state variable at the BEGINNING of the stage. It must consist of a keyword argument.

```julia
@state(sp, x >= 0.5, x0=1)
```
Both `x` and `x0` are JuMP variables.

Alternatively, we can use indexing just as we would in a JuMP `@variable` macro:
```julia
X0 = [3., 2.]
@state(sp, x[i=1:2], x0=X0[i])
```
In this case, both `x` and `x0` are JuMP dicts that can be indexed with the keys `1` and `2`.
All the indices must be specified in the second argument, but they can be referred to in the third argument. The indexing of `x0` will be identical to that of `x.`

#### Stage Profit
Define the stage profit for the stage problem `sp` using the `@stageprofit` macro. If our stage objective is `min cᵀx+θ` then:
```julia
@stageprofit(sp, dot(c, x))
```
This can be any valid JuMP syntax. The value/cost to go is handled automatically by the solver.

#### Scenario constraints
You can add scenarios to the stage problem with the `@scenarioconstraint` macro. It has three arguments
1. stage problem model object
2. keyword definition of scenario set. Length of the set must be identical to the value of the `scenarios` keyword in `SDDPModel()`
3. JuMP constraint. Can be any valid JuMP constraint syntax. It must include the keyword from the second argument

```julia
@scenario(sp, RHS=rand(3), a*x <= b + RHS)

RHS = rand(3)
@scenario(sp, i=1:3, a*x <= b + c*RHS[i])
```

If you add multiple scenario constraints the length of each scenario set must be identical. Scenario One corresponds to `RHS=Ω[1,i]` for all `i`.

```julia
Ω = rand(3,4)
for i=1:4
  @scenario(sp, RHS=Ω[:,i], a*x[i] <= b + RHS)
end
```

There is also a plural version of the macro similar to JuMP's `@constraints`:

```julia
Ω = rand(3,2)
@scenarios(sp, i=1:3, begin
  x[i] <= Ω[i, 1]

  y[i] <= Ω[i, 2]
end)
