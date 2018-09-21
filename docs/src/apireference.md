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
addconstraintnoise!
```

### Risk Measures
```@docs
AbstractRiskMeasure
modify_probability
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
store_cut
valid_cuts
all_cuts
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
