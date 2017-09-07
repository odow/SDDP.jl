```@meta
CurrentModule = SDDP
```

# API Reference

## Defining the model

```@docs
SDDPModel
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
```@docs
@stageobjective
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
```@docs
simulate
@visualise
getbound
```

## Save/Load
<!-- loadcuts! -->
```@docs
```
