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
stageobjective!
```

### Risk Measures
```@docs
AbstractRiskMeasure
Expectation
NestedAVaR
SDDP.modifyprobability!
```

### Cut Oracles
```@docs
AbstractCutOracle
DefaultCutOracle
LevelOneCutOracle
SDDP.storecut!!
SDDP.validcuts
```

## Solve
```@docs
solve
MonteCarloSimulation
BoundConvergence
Serial
Asyncronous
```
## Results
```@docs
getbound
simulate
@visualise
```

## Save/Load
```@docs
loadcuts!
```
