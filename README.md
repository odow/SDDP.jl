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

`SDDPModel` takes the following keyword arguments. The first two are essential:
+ `stages::Int`the number of stages in the model
+ `objective_bound::Float64` a valid bound on the initial value/cost to go. i.e. for maximisation this may be some large positive number, for minimisation this may be some large negative number.

The following arguments are optional:
- `sense`: must be either `:Max` or `:Min`. Defaults to `:Min`.
- `scenario_probability`: support vector for the scenarios. Defaults to uniform distribution. This can either be
    1. a vector that has the same length as the number of scenarios (identical in each stage). In which case `scenario_probability[s]` is the probability of scenario `s` occuring in each stage;
    2. a vector containing a probability vector for each stage. In which case `scenario_probability[t][s]` is the probability of scenario `s` occuring in stage `t`;
    3. a vector (with the same length as the number of stages) of vectors (where `scenario_probability[t]` has the same length as the number of markov states in stage `t`) of probability vectors. In which case `scenario_probability[t][i][s]` is the probability of scenario `s` occurring in markov state `i` in stage `t`.
- `markov_transition`: Transition probabilities for markov chain. Either a square matrix of size `markov_states` or a vector of such matrices with length `stages`. Defaults to uniform transition probability.
- `solver`: MathProgBase compliant solver that returns duals from a linear program. If this isn't specified then you must using `JuMP.setsolver(m, solver)` in the stage definition.
- `risk_measure`
- `cut_oracle`
- `value_function` Ignore this for now. Future proof's the extensibility of the library.

### A Note on Value Functions

You may notice we parameterise the SDDPModel by the DefaultValueFunction. Although
this is the only value function provided in this package, it enables extensibility
for some of our research codes that are not yet at the point for public release.
