# Tutorial Eight: odds and ends

In our previous tutorials, we have discussed the formulation, solution, and
visualization of a multistage stochastic optimization problem. By now you should
have all the tools necessary to solve your own problems using SDDP.jl. However,
there are a few odds and ends that we need to address.

## The `@state` macro

In [Tutorial One: first steps](@ref), we introduced a single state variable.
However, much more complex syntax is supported. This syntax hijacks JuMP's
`@variable` macro syntax, so it will be familiar to users of JuMP, but may look
unusual at first glance.
```julia
RESERVOIRS = [ :upper, :middle, :lower ]
V0 = Dict(
    :upper  = 10.0,
    :middle =  5.0,
    :lower  =  2.5
)
@state(m, volume[reservoir=RESERVOIRS], initial_volume == V0[reservoir])
```
This will create JuMP variables that can be accessed as `volume[:upper]` and
`initial_volume[:lower]`. There is also [`@states`](@ref), an analogue to JuMP's
`@variables` macro. It can be used as follows:
```julia
@states(sp, begin
    x >= 0,   x0==1
    y[i=1:2], y0==i
end)
```

## Multiple RHS noise constraints

In [Tutorial Two: RHS noise](@ref), we added a single constraint with noise in
the RHS term; however, you probably want to add many of these constraints. There
are two ways to do this. First, you can just add multiple calls like:
```julia
@rhsnoise(sp, w=[1,2,3], x <= w)
@rhsnoise(sp, w=[4,5,6], y >= w)
```
If you have multiple calls to `@rhsnoise`, they must have the same number of
elements in their sample space. For example, the following will not work:
```julia
@rhsnoise(sp, w=[1,2,3], x <= w)    # 3 elements
@rhsnoise(sp, w=[4,5,6,7], y >= w)  # 4 elements
```
Another option is to use the [`@rhsnoises`](@ref) macro. It is very similar to
[`@states`](@ref) and JuMP's `@consraints` macro:
```julia
@rhsnoises(sp, w=[1,2,3], begin
    x <= w
    y >= w + 3
end)
```

## `if` statements in more detail

In [Tutorial Four: Markovian policy graphs](@ref), we used an `if` statement to
control the call to `setnoiseprobability!` depending on the Markov state.
This can be used much more generally. For example:
```julia
m = SDDPModel(
    #...arguments omitted ...
        ) do sp, t, i
    @state(m, x>=0, x0==0)
    if t == 1
        @stageobjective(m, x)
    else
        if i == 1
            @variable(m, u >= 1)
            @constraint(m, x0 + u == x)
        else
            @rhsnoise(m, w=[1,2], x0 + w == x)
        end
        @stageobjective(m, x0 + x)
    end
end
```

You could, of course, do the following as well:
```julia
function first_stage(sp::JuMP.Model, x0, x)
    @stageobjective(m, x)
end
function second_stage(sp::JuMP.Model, x0, x)
    if i == 1
        @variable(m, u >= 1)
        @constraint(m, x0 + u == x)
    else
        @rhsnoise(m, w=[1,2], x0 + w == x)
    end
    @stageobjective(m, x0 + x)
end
m = SDDPModel(
    #...arguments omitted ...
        ) do sp, t, i
    @state(m, x>=0, x0==0)
    if t == 1
        first_stage(m, x0, x)
    else
        second_stage(m, x0, x)
    end
end
```

## Stage-dependent risk measures

In [Tutorial Five: risk](@ref), we discussed adding risk measures to our model.
We used the same risk measure for every stage. However, often you may want to
have time-varying risk measures. This can be accomplished very easily: instead
of passing a risk measure to the `risk_measure` keyword argument, we pass a
function that takes two inputs: `t::Int`, the index of the stage; and `i::Int`,
the index of the Markov state. For example:
```julia
function build_my_risk_measure(t::Int, i::Int)
    if t == 1
        return Expectation()
    elseif i == 1
        return AVaR(0.5)
    else
        return 0.5 * Expectation() + 0.5 * WorstCase()
    end
end

m = SDDPModel(
    # ... arguments omitted ...
    risk_measure = build_my_risk_measure
                                ) do sp, t, i
    # ... formulation omitted ...
end
```

## Reading and writing cuts to file

It's possible to save the cuts that are discovered during a solve to a file so
that they can be read back in or used for other analysis. This can be done using
the `cut_output_file` to [`solve`](@ref):
```julia
m = build_model()
SDDP.solve(m,
    iteration_limit = 10,
    cut_output_file = "cuts.csv"
)
```
`cuts.csv` is a csv file with `N+3` columns, where `N` is the number of state
variables in the model. The columns are: (1) the index of the stage; (2) the
index of the Markov state; (3) the intercept; and (4+), the coefficients of the
cuts for each of the state variables.

This cut file can be read back into a model using [`loadcuts!`](@ref):
```julia
m2 = build_model()
loadcuts!(m2, "cuts.csv")
```

Another option is to use the [`SDDP.writecuts!`](@ref) method after the model
has been solved:
```julia
m = build_model()
SDDP.solve(m, iteration_limit = 10)
SDDP.writecuts!("cuts.csv", m)
```

## Timer outputs

SDDP.jl embeds the great [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl)
package to help profile where time is spent during the solution process. You can
print the timing statistics by setting `print_level=2` in a [`solve`](@ref)
call. This produces a log like:

```
-------------------------------------------------------------------------------
                          SDDP.jl © Oscar Dowson, 2017-2018
-------------------------------------------------------------------------------
    Solver:
        Serial solver
    Model:
        Stages:         3
        States:         1
        Subproblems:    5
        Value Function: Default
-------------------------------------------------------------------------------
              Objective              |  Cut  Passes    Simulations   Total
     Simulation       Bound   % Gap  |   #     Time     #    Time    Time
-------------------------------------------------------------------------------
       27.000K         5.818K        |     1    0.0      0    0.0    0.0
        2.000K         6.952K        |     2    0.0      0    0.0    0.0
        2.000K         6.952K        |     3    0.0      0    0.0    0.0
       11.000K         7.135K        |     4    0.0      0    0.0    0.0
        2.000K         7.135K        |     5    0.0      0    0.0    0.0
        2.000K         7.135K        |     6    0.0      0    0.0    0.0
        5.000K         7.135K        |     7    0.0      0    0.0    0.0
        2.000K         7.135K        |     8    0.0      0    0.0    0.0
        5.000K         7.135K        |     9    0.0      0    0.0    0.0
       12.500K         7.135K        |    10    0.0      0    0.0    0.0
-------------------------------------------------------------------------------
 ─────────────────────────────────────────────────────────────────────────────────
         Timing statistics                Time                   Allocations
                                  ──────────────────────   ───────────────────────
         Tot / % measured:            40.8ms / 98.0%           0.97MiB / 100%

 Section                  ncalls     time   %tot     avg     alloc   %tot      avg
 ─────────────────────────────────────────────────────────────────────────────────
 Solve                         1   40.0ms   100%  40.0ms   0.96MiB  100%   0.96MiB
   Iteration Phase            10   38.3ms  95.6%  3.83ms    950KiB  96.3%  95.0KiB
     Backward Pass            10   30.5ms  76.3%  3.05ms    784KiB  79.4%  78.4KiB
       JuMP.solve            150   23.8ms  59.4%   159μs    423KiB  42.9%  2.82KiB
         optimize!           150   19.7ms  49.2%   131μs         -  0.00%        -
         prep JuMP model     150   2.69ms  6.72%  17.9μs    193KiB  19.6%  1.29KiB
         getsolution         150   1.09ms  2.71%  7.23μs    209KiB  21.1%  1.39KiB
       Cut addition           30   1.12ms  2.81%  37.4μs   69.1KiB  7.01%  2.30KiB
         risk measure         30   27.5μs  0.07%   917ns   8.91KiB  0.90%        -
     Forward Pass             10   7.68ms  19.2%   768μs    163KiB  16.5%  16.3KiB
       JuMP.solve             30   5.89ms  14.7%   196μs   93.4KiB  9.47%  3.11KiB
         optimize!            30   4.59ms  11.5%   153μs         -  0.00%        -
         prep JuMP model      30    909μs  2.27%  30.3μs   45.6KiB  4.62%  1.52KiB
         getsolution          30    303μs  0.76%  10.1μs   41.6KiB  4.21%  1.39KiB
 ─────────────────────────────────────────────────────────────────────────────────
    Other Statistics:
        Iterations:         10
        Termination Status: iteration_limit
===============================================================================
```

## Discount factors

If a model has a large number of stages, a common modelling trick is to add a
discount factor to the cost-to-go function. The stage-objective can then be
written as ``c^\\top + \\rho \\theta``, where `\\theta` is the cost-to-go
variable. However, since SDDP.jl does not expose the cost-to-go variable to the
user, this modelling trick must be accomplished as follows:
```julia
SDDPModel() do sp, t
    ρ = 0.95
    # ... lines omitted ...
    @stageobjective(sp, ρ^(t - 1) * dot(c, x))
end
```

That concludes our eighth tutorial for SDDP.jl. In our next tutorial,
[Tutorial Nine: nonlinear models](@ref), we discuss how SDDP.jl can be used to
solve problems that have nonlinear transition functions.
