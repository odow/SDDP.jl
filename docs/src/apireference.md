```@meta
CurrentModule = SDDP
```

# API Reference

## Defining the model

<!-- SDDPModel -->
```@docs
```

### States
```@docs
@state
@states
```
### Noises
```@docs
@noise
@noises
setnoiseprobability!
```

### Objective
<!-- stageobjective! -->
```@docs
```

### Risk Measures
```@docs
AbstractRiskMeasure
Expectation
NestedAVaR
modifyprobability!
```

### Cut Oracles
```@docs
AbstractCutOracle
DefaultCutOracle
LevelOneCutOracle
storecut!
validcuts
```

## Solve
<!-- SDDP.solve -->
```@docs
MonteCarloSimulation
BoundConvergence
```
## Results
<!-- simulate -->
<!-- @visualise -->
```@docs
getbound
```

## Save/Load
<!-- loadcuts! -->
```@docs
```
