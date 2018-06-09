# Tutorial Five: risk

Over the previous four tutorials, we formulated a simple hydrothermal scheduling
problem. Now, in this tutorial, we introduce some *risk* into the model using
nested risk measures.

Recall that our model for the hydrothermal scheduling problem from
[Tutorial Four: Markovian policy graphs](@ref) is:

```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0,
      markov_transition = Array{Float64, 2}[
          [ 1.0 ]',
          [ 0.75 0.25 ],
          [ 0.75 0.25 ; 0.25 0.75 ]
      ]
                                        ) do sp, t, i
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
        mupliplier * fuel_cost[t] * thermal_generation
    )
    if i == 1  # wet climate state
        setnoiseprobability!(sp, [1/6, 1/3, 0.5])
    else       # dry climate state
        setnoiseprobability!(sp, [0.5, 1/3, 1/6])
    end
end
```

## Formulating the problem

For this problem, we are going to use a convex combination of the expectation
(\$\\mathbb{E}\$) and average value-at-risk measures (AV@R\${}_{1-\\beta}\$).
In particular, we use AV@R at the β=0.1 quantile (i.e. the worst 10% of
outcomes). This can be constructed as:
```julia
risk_measure = 0.5 * Expectation() + 0.5 * AVaR(0.1)
```

Since this is a commonly used risk measure, a slightly more computationally
efficient form is [`EAVaR`](@ref):
```julia
risk_measure = EAVaR(lambda=0.5, beta=0.1)
```
This is short-hand for `lambda * Expectation() + (1-lambda) * AVaR(beta)`. As
`lambda` and `beta` tend toward `1.0`, the measure becomes more risk-neutral
(i.e. less risk averse).

Risk measures are set in the model using the `risk_measure` keyword in the
[`SDDPModel`](@ref) constructor. For example, our model is now:

```julia
m = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0.0,
      markov_transition = Array{Float64, 2}[
          [ 1.0 ]',
          [ 0.75 0.25 ],
          [ 0.75 0.25 ; 0.25 0.75 ]
      ],
           risk_measure = EAVaR(lambda=0.5, beta=0.1)
                                        ) do sp, t, i
    @state(sp, 0 <= outgoing_volume <= 200, incoming_volume == 200)
    @variables(sp, begin
        thermal_generation >= 0
        hydro_generation   >= 0
        hydro_spill        >= 0
    end)
    @rhsnoise(sp, inflow = [0.0, 50.0, 100.0],
        outgoing_volume - (incoming_volume - hydro_generation - hydro_spill) == inflow
    )
    @constraints(sp, begin
        thermal_generation + hydro_generation == 150
    end)
    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, mupliplier = [1.2, 1.0, 0.8],
        mupliplier * fuel_cost[t] * thermal_generation
    )
    if i == 1  # wet climate state
        setnoiseprobability!(sp, [1/6, 1/3, 0.5])
    else       # dry climate state
        setnoiseprobability!(sp, [0.5, 1/3, 1/6])
    end
end
```

## Solving the problem

Deciding when to terminate a risk-averse SDDP model is an unresolved problem in
the literature. In addition to the termination methods discussed in previous
tutorials, the user can also terminate the solve using the keys `[CRTL]+[C]`.

In addition, to demonstrate that we cannot use Monte Carlo simulation to
estimate the upper bound of a risk-averse model, we perform a Monte Carlo
simulation of the policy every two iterations. If left unchecked, this solve
will not terminate as we have set `terminate=false`.

```julia
status = solve(m,
    simulation = MonteCarloSimulation(
        frequency  = 2,
        confidence = 0.95,
        terminate  = false,
        min        = 50,
        step       = 50,
        max        = 100,
    )
)
```
After 7 iterations, we interrupt the solve using the `[CRTL]+[C]` keys. (Since
this is a trivial model to solve, we had to be very quick to terminate!) The
return value `status` is `:interrupted`, and the log is:
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
       21.000K         8.520K        |     1    0.0      0    0.0    0.0
   6.857K    9.363K   12.465K -45.0  |     2    0.0     50    0.0    0.0
        7.000K        12.465K        |     3    0.0     50    0.0    0.1
   7.911K   11.189K   12.465K -36.5  |     4    0.0    100    0.1    0.1
        8.000K        12.477K        |     5    0.0    100    0.1    0.1
   7.490K   10.230K   12.477K -40.0  |     6    0.0    150    0.1    0.1
        2.000K        12.477K        |     7    0.0    150    0.1    0.1
WARNING: Terminating solve due to user interaction
-------------------------------------------------------------------------------
    Other Statistics:
        Iterations:         7
        Termination Status: interrupted
===============================================================================
```

## Extra for experts: new risk measures

One of the cool features of SDDP.jl is how easy it is to create new risk
measures. To illustrate this, we consider implementing the worst-case risk
measure. First, we need to create a new concrete subtype of the abstract type
`AbstractRiskMeasure` defined by SDDP.jl:

```julia
"""
    TheWorstCase()

Create an instance of the worst-case risk measures. This places all of the
weight on the maximum outcome if minimizing, and the minimum outcome if
maximizing.
"""
struct TheWorstCase <: SDDP.AbstractRiskMeasure end
```

Then, we need to overload the [`SDDP.modifyprobability!`](@ref) function
provided by SDDP.jl. This function takes six arguments:
 1. an instance of the risk measure (e.g. `TheWorstCase()`);
 2. a vector of the risk-adjusted probability distribution that the function
    modifies in-place;
 3. a vector of the original probability distribution;
 4. a vector of the observations;
 5. the SDDPModel `m`; and
 6. the JuMP subproblem `sp`.

For example, the worst-case risk measure places all of the probability on the
worst outcome:
```julia
function SDDP.modifyprobability!(::TheWorstCase,
    risk_adjusted_distribution,
    original_distribution::Vector{Float64},
    observations::Vector{Float64},
    m::SDDPModel,
    sp::JuMP.Model
    )
    if getsense(sp) == :Min
        worst_index = indmax(observations)
    else
        worst_index = indmin(observations)
    end
    risk_adjusted_distribution .= 0.0
    risk_adjusted_distribution[worst_index] = 1.0
    return nothing
end
```
!!! note
    This implementation isn't a proper implementation as it assumes that the
    worst-case outcome has a positive probability of occurring. Accounting for
    this edge-case efficiently makes the implementation to verbose for this
    simple example.

Now `TheWorstCase()` can be used like a risk measure defined by SDDP.jl. It is
even possible to compose it with other risk measures, for example:
```julia
risk_measure = 0.5 * Expectation() + 0.5 * TheWorstCase()
```

This concludes our fifth tutorial for SDDP.jl. In the next tutorial,
[Tutorial Six: cut selection](@ref), we introduce cut selection.
