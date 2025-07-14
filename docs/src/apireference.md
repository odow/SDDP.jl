# [API Reference](@id api_reference_list)

This page lists the public API of SDDP.jl. Any functions in SDDP that are not
listed here are considered part of the private API and may change in any future
release.

!!! info
    This page is a semi-structured list of the SDDP.jl API. For a more
    structured overview, read the How-to guides or Tutorial parts of this
    documentation.

Load SDDP using:
```julia
using SDDP
```

SDDP exports only [`@stageobjective`](@ref). Therefore, all other calls must be
prefixed with `SDDP.`.

## `Graph`
```@docs
SDDP.Graph
```

## `add_node`
```@docs
SDDP.add_node
```

## `add_edge`
```@docs
SDDP.add_edge
```

## `add_ambiguity_set`
```@docs
SDDP.add_ambiguity_set
```

## `LinearGraph`
```@docs
SDDP.LinearGraph
```

## `MarkovianGraph`
```@docs
SDDP.MarkovianGraph
```

## `UnicyclicGraph`
```@docs
SDDP.UnicyclicGraph
```

## `LinearPolicyGraph`
```@docs
SDDP.LinearPolicyGraph
```

## `MarkovianPolicyGraph`
```@docs
SDDP.MarkovianPolicyGraph
```

## `PolicyGraph`
```@docs
SDDP.PolicyGraph
```

## `@stageobjective`
```@docs
@stageobjective
```

## `parameterize`
```@docs
SDDP.parameterize
```

## `add_objective_state`
```@docs
SDDP.add_objective_state
```

## `objective_state`
```@docs
SDDP.objective_state
```

## `Noise`
```@docs
SDDP.Noise
```

## `numerical_stability_report`
```@docs
SDDP.numerical_stability_report
```

## `train`
```@docs
SDDP.train
```

## `termination_status`
```@docs
SDDP.termination_status
```

## `write_cuts_to_file`
```@docs
SDDP.write_cuts_to_file
```

## `read_cuts_from_file`
```@docs
SDDP.read_cuts_from_file
```

## `write_log_to_csv`
```@docs
SDDP.write_log_to_csv
```

## `set_numerical_difficulty_callback`
```@docs
SDDP.set_numerical_difficulty_callback
```

## `AbstractStoppingRule`
```@docs
SDDP.AbstractStoppingRule
```

## `stopping_rule_status`
```@docs
SDDP.stopping_rule_status
```

## `convergence_test`
```@docs
SDDP.convergence_test
```

## `IterationLimit`
```@docs
SDDP.IterationLimit
```

## `TimeLimit`
```@docs
SDDP.TimeLimit
```

## `Statistical`
```@docs
SDDP.Statistical
```

## `BoundStalling`
```@docs
SDDP.BoundStalling
```

## `StoppingChain`
```@docs
SDDP.StoppingChain
```

## `SimulationStoppingRule`
```@docs
SDDP.SimulationStoppingRule
```

## `FirstStageStoppingRule`
```@docs
SDDP.FirstStageStoppingRule
```

## `AbstractSamplingScheme`
```@docs
SDDP.AbstractSamplingScheme
```

## `sample_scenario`
```@docs
SDDP.sample_scenario
```

## `InSampleMonteCarlo`
```@docs
SDDP.InSampleMonteCarlo
```

## `OutOfSampleMonteCarlo`
```@docs
SDDP.OutOfSampleMonteCarlo
```

## `Historical`
```@docs
SDDP.Historical
```

## `PSRSamplingScheme`
```@docs
SDDP.PSRSamplingScheme
```

## `SimulatorSamplingScheme`
```@docs
SDDP.SimulatorSamplingScheme
```

## `AbstractParallelScheme`
```@docs
SDDP.AbstractParallelScheme
```

## `Serial`
```@docs
SDDP.Serial
```

## `Threaded`
```@docs
SDDP.Threaded
```

## `Asynchronous`
```@docs
SDDP.Asynchronous
```

## `AbstractForwardPass`
```@docs
SDDP.AbstractForwardPass
```

## `DefaultForwardPass`
```@docs
SDDP.DefaultForwardPass
```

## `RevisitingForwardPass`
```@docs
SDDP.RevisitingForwardPass
```

## `RiskAdjustedForwardPass`
```@docs
SDDP.RiskAdjustedForwardPass
```

## `AlternativeForwardPass`
```@docs
SDDP.AlternativeForwardPass
```

## `AlternativePostIterationCallback`
```@docs
SDDP.AlternativePostIterationCallback
```

## `RegularizedForwardPass`
```@docs
SDDP.RegularizedForwardPass
```

## `ImportanceSamplingForwardPass`
```@docs
SDDP.ImportanceSamplingForwardPass
```

## `AbstractRiskMeasure`
```@docs
SDDP.AbstractRiskMeasure
```

## `adjust_probability`
```@docs
SDDP.adjust_probability
```

## `Expectation`
```@docs
SDDP.Expectation
```

## `WorstCase`
```@docs
SDDP.WorstCase
```

## `AVaR`
```@docs
SDDP.AVaR
```

## `CVaR`
```@docs
SDDP.CVaR
```

## `ConvexCombination`
```@docs
SDDP.ConvexCombination
```

## `EAVaR`
```@docs
SDDP.EAVaR
```

## `ModifiedChiSquared`
```@docs
SDDP.ModifiedChiSquared
```

## `Entropic`
```@docs
SDDP.Entropic
```

## `Wasserstein`
```@docs
SDDP.Wasserstein
```

## `AbstractDualityHandler`
```@docs
SDDP.AbstractDualityHandler
```

## `ContinuousConicDuality`
```@docs
SDDP.ContinuousConicDuality
```

## `LagrangianDuality`
```@docs
SDDP.LagrangianDuality
```

## `StrengthenedConicDuality`
```@docs
SDDP.StrengthenedConicDuality
```

## `BanditDuality`
```@docs
SDDP.BanditDuality
```

## `FixedDiscreteDuality`
```@docs
SDDP.FixedDiscreteDuality
```

## `simulate`
```@docs
SDDP.simulate
```

## `calculate_bound`
```@docs
SDDP.calculate_bound
```

## `add_all_cuts`
```@docs
SDDP.add_all_cuts
```

## Decision rules

## `DecisionRule`
```@docs
SDDP.DecisionRule
```

## `evaluate`
```@docs
SDDP.evaluate
```

## `SpaghettiPlot`
```@docs
SDDP.SpaghettiPlot
```

## `add_spaghetti`
```@docs
SDDP.add_spaghetti
```

## `publication_plot`
```@docs
SDDP.publication_plot
```

## `ValueFunction`
```@docs
SDDP.ValueFunction
```

## `evaluate`
```@docs
SDDP.evaluate(::SDDP.ValueFunction, ::Dict{Symbol,Float64})
```

## `plot`
```@docs
SDDP.plot
```

## `write_subproblem_to_file`
```@docs
SDDP.write_subproblem_to_file
```

## `deterministic_equivalent`
```@docs
SDDP.deterministic_equivalent
```

## `write_to_file`
```@docs
SDDP.write_to_file
```

## `read_from_file`
```@docs
SDDP.read_from_file
```

## `write`
```@docs
Base.write(::IO, ::SDDP.PolicyGraph)
```

## `read`
```@docs
Base.read(::IO, ::Type{SDDP.PolicyGraph})
```

## `evaluate`
```@docs
SDDP.evaluate(::SDDP.PolicyGraph{T}, ::SDDP.ValidationScenarios{T}) where {T}
```

## `ValidationScenarios`
```@docs
SDDP.ValidationScenarios
```

## `ValidationScenario`
```@docs
SDDP.ValidationScenario
```
