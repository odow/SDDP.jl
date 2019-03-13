# API Reference

## Policy graphs

```@docs
SDDP.Graph
SDDP.add_node
SDDP.add_edge
SDDP.LinearGraph
SDDP.MarkovianGraph
SDDP.LinearPolicyGraph
SDDP.MarkovianPolicyGraph
SDDP.PolicyGraph
```

## Subproblem definition

```@docs
@stageobjective
SDDP.parameterize
SDDP.add_objective_state
SDDP.objective_state
```

## Training the policy

```@docs
SDDP.numerical_stability_report
SDDP.train
SDDP.termination_status
```

## Simulating the policy

```@docs
SDDP.simulate
SDDP.calculate_bound
SDDP.Historical
```

## Visualizing the policy

```@docs
SDDP.SpaghettiPlot
SDDP.add_spaghetti
SDDP.publication_plot
```
