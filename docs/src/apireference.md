```@meta
CurrentModule = SDDP
```

# API Reference

## Communicating the problem to the solver

```@docs
SDDPModel
@state
@states
@rhsnoise
@rhsnoises
setnoiseprobability!
@stageobjective
```

### Risk Measures
```@docs
AbstractRiskMeasure
modifyprobability!
AVaR
ConvexCombination
Expectation
DRO
NestedAVaR
WorstCase
```

### Cut Oracles
```@docs
AbstractCutOracle
storecut!
validcuts
DefaultCutOracle
LevelOneCutOracle
```

## Solving the problem efficiently
```@docs
solve
MonteCarloSimulation
BoundConvergence
```
## Understanding the solution
```@docs
simulate
getbound
newplot
addplot!
show(::SDDP.SimulationPlot)
plotvaluefunction
```

## Read and write the model to disk

```@docs
loadcuts!
savemodel!
loadmodel
```
