# Add multi-dimensional noise terms

```@meta
DocTestSetup = quote
    using SDDP, HiGHS
end
```

Multi-dimensional stagewise-independent random variables can be created by
forming the Cartesian product of the random variables.

## A simple example

If the sample space and probabilities are given as vectors for each marginal
distribution, do:

```jldoctest; filter=[r"\(value = \d, coefficient = \d\)", r"1\-element.+"s]
julia> model = SDDP.LinearPolicyGraph(
           stages = 3,
           lower_bound = 0,
           optimizer = HiGHS.Optimizer,
       ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 0.0)
           Ω = [(value = v, coefficient = c) for v in [1, 2] for c in [3, 4, 5]]
           P = [v * c for v in [0.5, 0.5] for c in [0.3, 0.5, 0.2]]
           SDDP.parameterize(subproblem, Ω, P) do ω
               JuMP.fix(x.out, ω.value)
               @stageobjective(subproblem, ω.coefficient * x.out)
           end
       end;

julia> for s in SDDP.simulate(model, 1)[1]
           println("ω is: ", s[:noise_term_term])
       end
ω is: (value = 1, coefficient = 4)
ω is: (value = 1, coefficient = 3)
ω is: (value = 2, coefficient = 4)
```

## Using Distributions.jl

For sampling multidimensional random variates, it can be useful to use the
`Product` type from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

### Finite discrete distributions

There are several ways to go about this. If the sample space is finite, and
small enough that it makes sense to enumerate each element, we can use
`Base.product` and `Distributions.support` to generate the entire sample space
`Ω` from each of the marginal distributions.

We can then evaluate the density function of the product distribution on each
element of this space to get the vector of corresponding probabilities, `P`.

```jldoctest; filter=[r"\[\d+, \d+, \d+\]", r"1\-element.+"s]
julia> import Distributions

julia> distributions = [
           Distributions.Binomial(10, 0.5),
           Distributions.Bernoulli(0.5),
           Distributions.truncated(Distributions.Poisson(5), 2, 8)
       ];

julia> supports = Distributions.support.(distributions);

julia> Ω = vec([collect(ω) for ω in Base.product(supports...)]);

julia> P = [Distributions.pdf(Distributions.Product(distributions), ω) for ω in Ω];

julia> model = SDDP.LinearPolicyGraph(
           stages = 3,
           lower_bound = 0,
           optimizer = HiGHS.Optimizer,
       ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 0.0)
           SDDP.parameterize(subproblem, Ω, P) do ω
               JuMP.fix(x.out, ω[1])
               @stageobjective(subproblem, ω[2] * x.out + ω[3])
           end
       end;

julia> for s in SDDP.simulate(model, 1)[1]
           println("ω is: ", s[:noise_term])
       end
ω is: [10, 0, 3]
ω is: [0, 1, 6]
ω is: [6, 0, 5]
```

### Sampling

For sample spaces that are too large to explicitly represent, we can instead
approximate the distribution by a sample of `N` points. Now `Ω` is a sample from
the full sample space, and `P` is the uniform distribution over those points.
Points with higher density in the full sample space will appear more frequently
in `Ω`.

```jldoctest; filter=[r"\[\d+, \d+, \d+\]", r"1\-element.+"s]
julia> import Distributions

julia> distributions = Distributions.Product([
           Distributions.Binomial(100, 0.5),
           Distributions.Geometric(1 / 20),
           Distributions.Poisson(20),
       ]);

julia> N = 100;

julia> Ω = [rand(distributions) for _ in 1:N];

julia> P = fill(1 / N, N);

julia> model = SDDP.LinearPolicyGraph(
           stages = 3,
           lower_bound = 0,
           optimizer = HiGHS.Optimizer,
       ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 0.0)
           SDDP.parameterize(subproblem, Ω, P) do ω
               JuMP.fix(x.out, ω[1])
               @stageobjective(subproblem, ω[2] * x.out + ω[3])
           end
       end;

julia> for s in SDDP.simulate(model, 1)[1]
           println("ω is: ", s[:noise_term])
       end
ω is: [54, 38, 19]
ω is: [43, 3, 13]
ω is: [43, 4, 17]
```
