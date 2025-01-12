# Add a risk measure

```@meta
DocTestSetup = quote
    using SDDP, HiGHS
end
```
## Training a risk-averse model

`SDDP.jl` supports a variety of risk measures. Two common ones are
[`SDDP.Expectation`](@ref) and [`SDDP.WorstCase`](@ref). Let's see how to
train a policy using them. There are three possible ways.

If the same risk measure is used at every node in the policy graph, we can just
pass an instance of one of the risk measures to the `risk_measure` keyword
argument of the [`SDDP.train`](@ref) function.

```julia
SDDP.train(
    model,
    risk_measure = SDDP.WorstCase(),
    iteration_limit = 10
)
```

However, if you want different risk measures at different nodes, there are two
options. First, you can pass `risk_measure` a dictionary of risk measures,
with one entry for each node. The keys of the dictionary are the indices of the
nodes.

```julia
SDDP.train(
    model,
    risk_measure = Dict(
        1 => SDDP.Expectation(),
        2 => SDDP.WorstCase()
    ),
    iteration_limit = 10
)
```

An alternative method is to pass `risk_measure` a function that takes one
argument, the index of a node, and returns an instance of a risk measure:
```julia
SDDP.train(
    model,
    risk_measure = (node_index) -> begin
        if node_index == 1
            return SDDP.Expectation()
        else
            return SDDP.WorstCase()
        end
    end,
    iteration_limit = 10
)
```

!!! note
    If you simulate the policy, the simulated value is the risk-neutral value of
    the policy.

## Supported risk measures

The following risk measures are implemented in SDDP.jl:

 - [`Expectation`](@ref)  [default]
 - [`WorstCase`](@ref)
 - [`AVaR`](@ref) or [`CVaR`](@ref)
 - [`ConvexCombination`](@ref)
 - [`EAVaR`](@ref)
 - [`ModifiedChiSquared`](@ref)
 - [`Entropic`](@ref)
 - [`Wasserstein`](@ref)
