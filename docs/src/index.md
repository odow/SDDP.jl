```@meta
CurrentModule = SDDP
```

## Overview
*SDDP.jl - Stochastic Dual Dynamic Programming in Julia.*

Solving a mathematical optimization problem requires four steps:

1. Formulating the problem;
2. Communicating the problem to the solver;
3. Solving the problem efficiently;
4. Understanding the solution.

For deterministic problems, [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl)
provides a first-in-class solution to each of these steps.

SDDP.jl builds upon JuMP to provide a high-level interface to the current
state-of-the-art in order to make solving large multi-stage convex stochastic
optimization problems easy.

In this manual, we detail the many features of SDDP.jl through the classic
example of balancing a portfolio of stocks and bonds over time.

## Getting started

This package is unregistered so you will need to `Pkg.clone` it as follows:

```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```

If you want to use the parallel features of SDDP.jl, you should start Julia with
some worker processes (`julia -p N`), or add by running `julia> addprocs(N)` in
a running Julia session.

## Formulating the problem

### The Asset Management Problem

The goal of the asset management problem is to choose an investment portfolio
that is composed of stocks and bonds in order to meet a target wealth goal at
the end of the time horizon. After five, and ten years, the agent observes the
portfolio and is able to re-balance their wealth between the two asset classes.
As an extension to the original problem, we introduce two new random variables.
The first that represents a source of additional wealth in years 5 and 10. The
second is an immediate reward that the agent incurs for holding stocks at the
end of years 5 and 10. This can be though of as a dividend that cannot be
reinvested.

## Communicating the problem to the solver

The second step in the optimization process is communicating the problem to the
solver. To do this, we are going to build each subproblem as a JuMP model, and
provide some metadata that describes how the JuMP subproblems inter-relate.

### The Model Constructor

The core type of SDDP.jl is the `SDDPModel` object. It can be constructed with
```julia
m = SDDPModel( ... metadata as keyword arguments ... ) do sp, t, i
    ... JuMP subproblem definition ...
end
```
We draw the readers attention to three sections in the SDDPModel constructor.

### State Variables

We can define a new state variable in the stage problem `sp` using the `@state`
macro:

```julia
@state(sp, x >= 0.5, x0==1)
```

We can also use indexing just as we would in a JuMP `@variable` macro:
```julia
X0 = [3.0, 2.0]
@state(sp, x[i=1:2], x0==X0[i])
```
In this case, both `x` and `x0` are JuMP dicts that can be indexed with the keys
`1` and `2`. All the indices must be specified in the second argument, but they
can be referred to in the third argument. The indexing of `x0` will be identical
to that of `x.`

There is also a plural version of the `@state` macro
```julia
@states(sp, begin
    x >= 0.0, x0==1
    y >= 0.0, y0==1
end)
```

### Standard JuMP machinery

Remember that `sp` is just a normal JuMP model, and so (almost) anything you can
do in JuMP, you can do in SDDP.jl. The one exception is the objective, which we
detail in the next section.

However,, control variables are just normal JuMP variables and can be created
using `@variable` or `@variables`. Dynamical constraints, and feasiblity sets
can be specified using `@constraint` or `@constraints`.

### The stage objective

```julia
@stageobjective(sp, obj)
```

```julia
@stageobjective(sp, kw=realizations, obj)
setnoiseprobability!(sp, [0.2, 0.3, 0.5])
```

### Dynamics with Linear Noise

```julia
@rhsnoise(sp, w=[1,2,3], x <= w)
setnoiseprobability!(sp, [0.2, 0.3, 0.5])
```

### Asset Management Example

We now have all the information necessary to define the Asset Management example
in SDDP.jl:

```julia
using SDDP, JuMP, Clp

m = SDDPModel(
               # we are minimizing
                sense = :Min,
               # there are 4 stages
               stages = 4,
               # a large negative value
      objective_bound = -1000.0,
               # a MathOptBase compatible solver
               solver = ClpSolver(),
               # transition probabilities of the lattice
    markov_transition = Array{Float64, 2}[
                        [1.0]',
                        [0.5 0.5],
                        [0.5 0.5; 0.5 0.5],
                        [0.5 0.5; 0.5 0.5]
                    ],
               # risk measures for each stage
         risk_measure = [
                        Expectation(),
                        Expectation(),
                        NestedAVaR(lambda = 0.5, beta=0.5),
                        Expectation()
                    ]
                            ) do sp, t, i
    # Data used in the problem
    ωs = [1.25, 1.06]
    ωb = [1.14, 1.12]
    Φ  = [-1, 5]
    Ψ  = [0.02, 0.0]

    # state variables
    @states(sp, begin
        xs >= 0, xsbar==0
        xb >= 0, xbbar==0
    end)

    if t == 1 # we can switch based on the stage
        # a constraint without noise
        @constraint(sp, xs + xb == 55 + xsbar + xbbar)
        # an objective without noise
        @stageobjective(sp, 0)
    elseif t == 2 || t == 3
        # a constraint with noisein the RHS
        @rhsnoise(sp, φ=Φ, ωs[i] * xsbar + ωb[i] * xbbar + φ == xs + xb)
        # an objective with noise
        @stageobjective(sp, ψ = Ψ, -ψ * xs)
        # noise will be sampled as (Φ[1], Ψ[1]) w.p. 0.6, (Φ[2], Ψ[2]) w.p. 0.4
        setnoiseprobability!(sp, [0.6, 0.4])
    else # when t == 4
        # some control variables
        @variable(sp, u >= 0)
        @variable(sp, v >= 0)
        # dynamics constraint
        @constraint(sp, ωs[i] * xsbar + ωb[i] * xbbar + u - v == 80)
        # an objective without noise
        @stageobjective(sp, 4u - v)
    end
end
```

## Solving the problem efficiently

## Understanding the solution

## Extras for experts

### New risk measures

SDDP.jl makes it easy to create new risk measures. First, create a new subtype
of the abstract type `SDDP.AbstractRiskMeasure`:

```julia
immutable MyNewRiskMeasure <: SDDP.AbstractRiskMeasure
end
```

Then, overload the method `SDDP.modifyprobability!` for your new type.
`SDDP.modifyprobability!` has the following signature:

```julia
SDDP.modifyprobability!(
        measure::AbstractRiskMeasure,
        riskadjusted_distribution,
        original_distribution::Vector{Float64},
        observations::Vector{Float64},
        m::SDDPModel,
        sp::JuMP.Model
)
```
where `original_distribution` contains the risk netural probability of each
outcome in `observations` occurring (so that the probability of `observations[i]`
is `original_distribution[i]`). The method should modify (in-place) the elements
of `riskadjusted_distribution` to represent the risk-adjusted probabilities of
the distribution.

To illustrate this, we shall define the worst-case riskmeasure (which places all
the probability on the worst outcome).

```julia
immutable WorstCase <: SDDP.AbstractRiskMeasure end
function SDDP.modifyprobability!(measure::WorstCase,
        riskadjusted_distribution,
        original_distribution::Vector{Float64},
        observations::Vector{Float64},
        m::SDDPModel,
        sp::JuMP.Model
    )
    if JuMP.getobjectivesense(sp) == :Min
        # if minimizing, worst is the maximum outcome
        idx = indmax(observations)
    else
        # if maximizing, worst is the minimum outcome
        idx = indmin(observations)
    end
    # zero all probabilities
    riskadjusted_distribution .= 0.0
    # set worst to 1.0
    riskadjusted_distribution[idx] = 1.0
    # return
    return nothing
end
```

The risk measure `WorstCase()` can now be used in any SDDP model.

### New cut oracles
