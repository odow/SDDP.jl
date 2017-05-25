# SDDP

[![Build Status](https://travis-ci.com/odow/SDDP.jl.svg?token=BjRx6YCjMdN19LP812Rj&branch=master)](https://travis-ci.com/odow/SDDP.jl)

## Installation
This package is unregistered so you will need to `Pkg.clone` it as follow:
```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```

## Quick Start Guide
See the `/examples` folder for example usage

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
This function constructs an SDDPModel. `SDDPModel` takes the following keyword
arguments. Some are required, and some are optional.

#### Required Keyword arguments

 * `stages::Int`
 The number of stages in the problem. A stage is defined as each step in time at
 which a decion can be made. Defaults to `1`.

 * `objective_bound::Float64`
 A valid bound on the initial value/cost to go. i.e. for maximisation this may be some large positive number, for minimisation this may be some large negative number.

 * `solver::MathProgBase.AbstractMathProgSolver`
 MathProgBase compliant solver that returns duals from a linear program. If this isn't specified then you must use `JuMP.setsolver(m, solver)` in the stage definition.

#### Optional Keyword arguments
 * `sense`
 Must be either `:Max` or `:Min`. Defaults to `:Min`.

 * `cut_oracle::SDDP.AbstractCutOracle`
 The cut oracle is responsible for collecting and storing the cuts that define
 a value function. The cut oracle may decide that only a subset of the total
 discovered cuts are relevant, which improves solution speed by reducing the size
 of the subproblems that need solving. Currently must be one of
    * `DefaultCutOracle()`
    The default. All cuts are added to the subproblem.

    * `DematosCutOracle()`
    Implements De Matos Cut Selection (also called Level-One Cut Selection) as
    proposed in de Matos, V., Philpott, A., and Finardi, E.. “Improving the
    Performance of Stochastic Dual Dynamic Programming.” Journal of
    Computational and Applied Mathematics 290 (December 2015): 196–208.

 * `risk_measure::SDDP.AbstractRiskMeasure`
 A `SDDP.AbstractRiskMeasure` to use. Currently only `Expectation()` and
 `NestedAVaR` are valid.
    * `NestedAVaR(;lambda=0.0, beta=0.0)`
    A risk measure that is a convex combination of Expectation and Average Value @
Risk (also called Conditional Value @ Risk): λ * E[x] + (1 - λ) * AV@R(1-β)[x]

    It has the keyword arguments:

         * `lambda`
        Convex weight on the expectation (`(1-lambda)` weight is put on the AV@R component.
        Inreasing values of `lambda` are less risk averse (more weight on expecattion)

         * `beta`
         The quantile at which to calculate the Average Value @ Risk. Increasing values
         of `beta` are less risk averse. If `beta=0`, then the AV@R component is the
         worst case risk measure.


 * `scenario_probability`
 Support vector for the scenarios. Defaults to uniform distribution. This can
 either be

    1. a vector that has the same length as the number of scenarios (identical
    in each stage). In which case `scenario_probability[s]` is the probability
    of scenario `s` occuring in each stage;

    2. a vector containing a probability vector for each stage. In which case
    `scenario_probability[t][s]` is the probability of scenario `s` occuring in
    stage `t`;

    3. a vector (with the same length as the number of stages) of vectors (where
    `scenario_probability[t]` has the same length as the number of markov states
    in stage `t`) of probability vectors. In which case
    `scenario_probability[t][i][s]` is the probability of scenario `s` occurring
    in markov state `i` in stage `t`.

 * `markov_transition`

# Returns
    * `m`: the `SDDPModel`

### Solve

To solve the SDDP model `m` we use `status = solve(m::SDDPModel; kwargs...)`.
This accepts a number of keyword arguments to control the solution process.

#### Positional arguments
 * `m`: the SDDPModel to solve

#### Keyword arguments
 * `max_iterations::Int`:
    The maximum number of cuts to add to a single stage problem before terminating.
    Defaults to `10`.
 * `time_limit::Real`:
    The maximum number of seconds (in real time) to compute for before termination.
    Defaults to `Inf`.
 * `simulation::MonteCarloSimulation`:
    We control the behaviour of the policy simulation phase of the algorithm using
    the `MonteCarloSimulation(;kwargs...)` constructor. This just groups a
    series of related keyword arguments. The keywords are
    * `frequency::Int`
    The frequency (by iteration) with which to run the policy simulation phase of
    the algorithm in order to construct a statistical bound for the policy. Defaults
    to `0` (never run).
    * `min::Float64`
    Minimum number of simulations to conduct before constructing a confidence interval
    for the bound. Defaults to `20`.
    * `step::Float64`
    Number of additional simulations to conduct before constructing a new confidence
    interval for the bound. Defaults to `1`.
    * `max::Float64`
    Maximum number of simulations to conduct in the policy simulation phase. Defaults
    to `min`.
    * `confidence::Float64`
    Confidence level of the confidence interval. Defaults to `0.95` (95% CI).
    * `termination::Bool`
    Whether to terminate the solution algorithm with the status `:converged` if the
    deterministic bound is with in the statistical bound after `max` simulations.
    Defaults to `false`.
 * `bound_convergence`:
    We may also wish to terminate the algorithm if the deterministic bound stalls
    for a specified number of iterations (regardless of whether the policy has
    converged). This can be controlled by the `BoundConvergence(;kwargs...)`
    constructor. It has the following keywords:
    * `iterations::Int`
    Terminate if the maximum deviation in the deterministic bound from the mean
    over the last `iterations` number of iterations is less than `rtol` (in
    relative terms) or `atol` (in absolute terms).
    * `rtol::Float64`
    Maximum allowed relative deviation from the mean.
    Defaults to `0.0`
    * `atol::Float64`
    Maximum allowed absolute deviation from the mean.
    Defaults to `0.0`
 * `cut_selection_frequency::Int`:
    Frequency (by iteration) with which to rebuild subproblems using a subset of
    cuts. Frequent cut selection (i.e. `cut_selection_frequency` is small) reduces
    the size of the subproblems that are solved, but incurrs the overhead of rebuilding
    the subproblems. However, infrequent cut selection (i.e.
    `cut_selection_frequency` is large) allows the subproblems to grow large (many
    constraints) leading to an increase in the solve time of individual subproblems.
    Defaults to `0` (never run).
 * `print_level::Int`:
     0 - off: nothing logged to screen (still written to log file if specified).
     1 - on: solve iterations written to screen.
     Defaults to `1`
 * `log_file::String`:
    Relative filename to write the log to disk. Defaults to `""` (no log written)
 * `solve_type`:
    One of
    * `Asyncronous()` - solve using a parallelised algorithm
    * `Serial()` - solve using a serial algorithm
    Default chooses automatically based on the number of available processors.
 * `reduce_memory_footprint::Bool`:
    Implements the idea proposed in https://github.com/JuliaOpt/JuMP.jl/issues/969#issuecomment-282191105
    to reduce the memory consumption when running SDDP. This is an issue if you
    wish to save the model `m` to disk since it discards important information.
    Defaults to `false`.
 * `cut_output_file::String`:
    Relative filename to write discovered cuts to disk. Defaults to `""` (no cuts written)

#### Returns
 * `status::Symbol`:
    Reason for termination. One of
    * `:solving`
    * `:interrupted`
    * `:converged`
    * `:max_iterations`
    * `:bound_convergence`
    * `:time_limit`

### A Note on Value Functions

You may notice we parameterise the SDDPModel by the DefaultValueFunction. Although
this is the only value function provided in this package, it enables extensibility
for some of our research codes that are not yet at the point for public release.
