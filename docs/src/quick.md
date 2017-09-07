# Quick Start

```@meta
CurrentModule = SDDP
```

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

```@docs
SDDPModel
```

### Solve

To solve the SDDP model `m` we use `status = solve(m::SDDPModel; kwargs...)`.
This accepts a number of keyword arguments to control the solution process.

```@docs
solve
```

### Simulate

### Visualise

```@docs
@visualise
```
