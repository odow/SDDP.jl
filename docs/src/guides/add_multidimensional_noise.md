# Add multi-dimensional noise terms

```@meta
DocTestSetup = quote
    using SDDP, HiGHS, Distributions, IterTools
end
```

Multi-dimensional stagewise-independent random variables can be created by
forming the Cartesian product of the random variables.

```jldoctest; filter=[r"\(value = \d, coefficient = \d\)", r"1\-element.+"s]
julia> model = SDDP.LinearPolicyGraph(
               stages=3, lower_bound = 0, optimizer = HiGHS.Optimizer
               ) do subproblem, t
           @variable(subproblem, x, SDDP.State, initial_value = 0.0)
           support = [(value = v, coefficient = c) for v in [1, 2] for c in [3, 4, 5]]
           probability = [v * c for v in [0.5, 0.5] for c in [0.3, 0.5, 0.2]]
           SDDP.parameterize(subproblem, support, probability) do ω
               JuMP.fix(x.out, ω.value)
               @stageobjective(subproblem, ω.coefficient * x.out)
               println("ω is: ", ω)
           end
       end;

julia> SDDP.simulate(model, 1);
ω is: (value = 1, coefficient = 4)
ω is: (value = 1, coefficient = 3)
ω is: (value = 2, coefficient = 4)
```

For sampling multidimensional random variates, it can be useful to use the `Product` type from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package.

There are several ways to go about this. If the sample space is small enough that it makes sense to enumerate each element, we can use the `IterTools.product` method to generate the entire sample space `Ω` from each of the marginal distributions. We can then evaluate the density function of the product distribution on each element of this space to get the vector of corresponding probabilities, `P`. These two vectors may then be passed to `SDDP.parameterize`.

Because `Distributions.support` which is used to get the support of each marginal distribution only works with distributions types that have finite support, the code will error with distributions that are incompatible with SDDP.jl.

```jldoctest; filter=[r"\(\d, \d, \d\)", r"1\-element.+"s]
julia> distributions = Product([Binomial(10,0.5), Bernoulli(0.5), truncated(Poisson(5),2,8)]);

julia> Ω = IterTools.product(collect.(support.(distributions.v))...) |> collect |> vec;
julia> P = [pdf(distributions, collect(ω)) for ω in Ω];

julia> model = SDDP.LinearPolicyGraph(
                      stages=3, lower_bound = 0, optimizer = HiGHS.Optimizer
                      ) do subproblem, t
                  @variable(subproblem, x, SDDP.State, initial_value = 0.0)
                  SDDP.parameterize(subproblem, Ω, P) do ω
                      JuMP.fix(x.out, ω[1])
                      @stageobjective(subproblem, ω[2] * x.out + ω[3])
                      println("ω is: ", ω)
                  end
              end;

julia> SDDP.simulate(model, 1);
ω is: (5, 0, 6)
ω is: (5, 0, 7)
ω is: (4, 0, 2)
```

For sample spaces that are too large to explicitly represent, we can instead approximate the distribution by a sample of `N` points. Now `Ω` is a sample from the full sample space, and `P` is the uniform distribution over those points. Points with higher density in the full sample space will appear more frequently in `Ω`.

```jldoctest; filter=[r"\[\d, \d, \d\]", r"1\-element.+"s]
julia> distributions = Product([Binomial(100,0.5), Geometric(1/20), Poisson(20)]);
julia> N = 100;

julia> Ω = [rand(distributions) for _ in 1:N];
julia> P = fill(1 / N, N);

julia> model = SDDP.LinearPolicyGraph(
                      stages=3, lower_bound = 0, optimizer = HiGHS.Optimizer
                      ) do subproblem, t
                  @variable(subproblem, x, SDDP.State, initial_value = 0.0)
                  SDDP.parameterize(subproblem, Ω, P) do ω
                      JuMP.fix(x.out, ω[1])
                      @stageobjective(subproblem, ω[2] * x.out + ω[3])
                      println("ω is: ", ω)
                  end
              end;

julia> SDDP.simulate(model, 1);
ω is: [5, 0, 6]
ω is: [5, 0, 7]
ω is: [4, 0, 2]
```