# API Reference

## Policy graphs

```@docs
Kokako.LinearPolicyGraph
Kokako.MarkovianPolicyGraph
Kokako.PolicyGraph
```

## Subproblem definition

```@docs
@stageobjective
Kokako.parameterize
```

## Training the policy

```@docs
Kokako.train
Kokako.termination_status
```

### Risk measures

```@docs
Kokako.Expectation
Kokako.WorstCase
Kokako.AVaR
Kokako.EAVaR
Kokako.ModifiedChiSquared
Kokako.Wasserstein
```

## Simulating the policy

```@docs
Kokako.simulate
Kokako.calculate_bound
Kokako.Historical
```

## Visualizing the policy

```@docs
Kokako.SpaghettiPlot
Kokako.add_spaghetti
```
