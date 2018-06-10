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
EAVaR
Expectation
DRO
WorstCase
```

### Cut Oracles
```@docs
AbstractCutOracle
storecut!
validcuts
allcuts
DefaultCutOracle
LevelOneCutOracle
```

### Price Interpolation
```@docs
StaticPriceInterpolation
DynamicPriceInterpolation
DiscreteDistribution
```


## Solving the problem efficiently
```@docs
solve
MonteCarloSimulation
BoundStalling
Asynchronous
Serial
```
## Understanding the solution
```@docs
simulate
getbound
newplot
addplot!
show(::SDDP.SimulationPlot)
plotvaluefunction
getsubproblem
```

## Read and write the model to disk

```@docs
writecuts!
loadcuts!
savemodel!
loadmodel
```
