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
SDDP.train(model, iteration_limit = 10)
```

The training of an SDDP policy can be terminated after a fixed number of
seconds using the `time_limit` keyword.

```julia
SDDP.train(model, time_limit = 2.0)
```

## Stopping rules

In addition to the limits provided as keyword arguments, a variety of other
stopping rules are available. These can be passed to [`SDDP.train`](@ref)
as a vector to the `stopping_rules` keyword. For example:

```julia
SDDP.train(model, stopping_rules = [SDDP.BoundStalling(10, 1e-4)])
```

Here are the stopping rules implemented in `SDDP.jl`:

```@docs
SDDP.IterationLimit
SDDP.TimeLimit
SDDP.Statistical
SDDP.BoundStalling
```
