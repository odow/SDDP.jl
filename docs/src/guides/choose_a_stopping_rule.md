# Choose a stopping rule

The theory of SDDP tells us that the algorithm converges to an optimal policy
almost surely in a finite number of iterations. In practice, this number is very
large. Therefore, we need some way of pre-emptively terminating SDDP when the
solution is “good enough.” We call heuristics for pre-emptively terminating SDDP
_stopping rules_.

## Basic limits

The training of an SDDP policy can be terminated after a fixed number of
iterations using the `iteration_limit` keyword.

```julia
SDDP.train(model; iteration_limit = 10)
```

The training of an SDDP policy can be terminated after a fixed number of
seconds using the `time_limit` keyword.

```julia
SDDP.train(model; time_limit = 2.0)
```

## Stopping rules

In addition to the limits provided as keyword arguments, a variety of other
stopping rules are available. These can be passed to [`SDDP.train`](@ref)
as a vector to the `stopping_rules` keyword.  Training stops if any of the rules
becomes active. To stop when all of the rules become active, use
[`SDDP.StoppingChain`](@ref). For example:

```julia
# Terminate if BoundStalling becomes true
SDDP.train(
    model;
    stopping_rules = [SDDP.BoundStalling(10, 1e-4)],
)

# Terminate if BoundStalling OR TimeLimit becomes true
SDDP.train(
    model;
    stopping_rules = [SDDP.BoundStalling(10, 1e-4), SDDP.TimeLimit(100.0)],
)

# Terminate if BoundStalling AND TimeLimit becomes true
SDDP.train(
    model;
    stopping_rules = [
        SDDP.StoppingChain(SDDP.BoundStalling(10, 1e-4), SDDP.TimeLimit(100.0)),
    ],
)
```

## Supported rules

The stopping rules implemented in SDDP.jl are:

 - [`SDDP.IterationLimit`](@ref)
 - [`SDDP.TimeLimit`](@ref)
 - [`SDDP.Statistical`](@ref)
 - [`SDDP.BoundStalling`](@ref)
 - [`SDDP.StoppingChain`](@ref)
 - [`SDDP.SimulationStoppingRule`](@ref)
 - [`SDDP.FirstStageStoppingRule`](@ref)
