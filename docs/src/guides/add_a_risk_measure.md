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

## Risk measures

To illustrate the risk-measures included in `SDDP.jl`, we consider a discrete
random variable with four outcomes.

The random variable is supported on the values 1, 2, 3, and 4:

```@repl intermediate_risk
noise_supports = [1, 2, 3, 4]
```

The associated probability of each outcome is as follows:

```@repl intermediate_risk
nominal_probability = [0.1, 0.2, 0.3, 0.4]
```

With each outcome ω, the agent observes a cost `Z(ω)`:
```@repl intermediate_risk
cost_realizations = [5.0, 4.0, 6.0, 2.0]
```

We assume that we are minimizing:
```@repl intermediate_risk
is_minimization = true
```

Finally, we create a vector that will be used to store the risk-adjusted
probabilities:

```@repl intermediate_risk
risk_adjusted_probability = zeros(4)
```

### Expectation

```@docs
SDDP.Expectation
```

```@repl intermediate_risk
using SDDP
SDDP.adjust_probability(
    SDDP.Expectation(),
    risk_adjusted_probability,
    nominal_probability,
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```

[`SDDP.Expectation`](@ref) is the default risk measure in `SDDP.jl`.

### Worst-case

```@docs
SDDP.WorstCase
```

```@repl intermediate_risk
SDDP.adjust_probability(
    SDDP.WorstCase(),
    risk_adjusted_probability,
    nominal_probability,
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```

### Average value at risk (AV@R)

```@docs
SDDP.AVaR
```

```@repl intermediate_risk
SDDP.adjust_probability(
    SDDP.AVaR(0.5),
    risk_adjusted_probability,
    nominal_probability,
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```

### Convex combination of risk measures

Using the axioms of coherent risk measures, it is easy to show that any convex
combination of coherent risk measures is also a coherent risk measure. Convex
combinations of risk measures can be created directly:

```@repl intermediate_risk
cvx_comb_measure = 0.5 * SDDP.Expectation() + 0.5 * SDDP.WorstCase()
SDDP.adjust_probability(
    cvx_comb_measure,
    risk_adjusted_probability,
    nominal_probability,
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```

As a special case, the [`SDDP.EAVaR`](@ref) risk-measure is a convex
combination of [`SDDP.Expectation`](@ref) and [`SDDP.AVaR`](@ref):
```@repl intermediate_risk
SDDP.EAVaR(beta=0.25, lambda=0.4)
```

```@docs
SDDP.EAVaR
```

### Distributionally robust

`SDDP.jl` supports two types of distributionally robust risk measures: the
modified Χ² method of Philpott et al. (2018), and a method based on the
Wasserstein distance metric.

#### Modified Chi-squard

```@docs
SDDP.ModifiedChiSquared
```

```@repl intermediate_risk
SDDP.adjust_probability(
    SDDP.ModifiedChiSquared(0.5),
    risk_adjusted_probability,
    [0.25, 0.25, 0.25, 0.25],
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```

#### Wasserstein

```@docs
SDDP.Wasserstein
```

```@repl intermediate_risk
import HiGHS
SDDP.adjust_probability(
    SDDP.Wasserstein(HiGHS.Optimizer; alpha=0.5) do x, y
        return abs(x - y)
    end,
    risk_adjusted_probability,
    nominal_probability,
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```

### Entropic

```@docs
SDDP.Entropic
```

```@repl intermediate_risk
SDDP.adjust_probability(
    SDDP.Entropic(0.1),
    risk_adjusted_probability,
    nominal_probability,
    noise_supports,
    cost_realizations,
    is_minimization
)
risk_adjusted_probability
```
