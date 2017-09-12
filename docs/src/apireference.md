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
```@docs
solve
MonteCarloSimulation
BoundConvergence
```
## Results
```@docs
simulate
getbound
newplot
addplot!
show
plotvaluefunction
```

## Save/Load
<!-- loadcuts! -->
```@docs
```
