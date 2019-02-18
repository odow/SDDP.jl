# Intermediate II: stopping rules

The theory of SDDP tells us that the algorithm converges to an optimal policy
almost surely in a finite number of iterations. In practice, this number is very
large. Therefore, we need some way of pre-emptively terminating SDDP when the
solution is “good enough.” We call heuristics for pre-emptively terminating SDDP
_stopping rules_.

## Basic limits

The training of an SDDP policy can be terminated after a fixed number of
iterations using the `iteration_limit` keyword.

```julia
Kokako.train(model, iteration_limit = 10)
```

The training of an SDDP policy can be terminated after a fixed number of
seconds using the `time_limit` keyword.

```julia
Kokako.train(model, time_limit = 2.0)
```

## Stopping rules

In addition to the limits provided as keyword arguments, a variety of other
stopping rules are available. These can be passed to [`Kokako.train`](@ref)
as a vector to the `stopping_rules` keyword. For example:

```julia
Kokako.train(model, stopping_rules = [Kokako.BoundStalling(10, 1e-4)])
```

Here are the stopping rules implemented in `Kokako.jl`:

```@docs
Kokako.IterationLimit
Kokako.TimeLimit
Kokako.Statistical
Kokako.BoundStalling
```

In the next tutorial, [Intermediate III: policy graphs](@ref), we discuss
generic extensions to [`Kokako.LinearPolicyGraph`](@ref) and
[`Kokako.MarkovianPolicyGraph`](@ref).
