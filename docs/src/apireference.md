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
SDDP.write_cuts_to_file
SDDP.read_cuts_from_file
```

### Stopping rules

```@docs
SDDP.AbstractStoppingRule
SDDP.stopping_rule_status
SDDP.convergence_test
```

### Sampling schemes

```@docs
SDDP.AbstractSamplingScheme
SDDP.sample_scenario
```

### Sampling schemes

```@docs
SDDP.AbstractRiskMeasure
SDDP.adjust_probability
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
## Debugging the model

```@docs
SDDP.write_subproblem_to_file
```
