```@meta
CurrentModule = SDDP
```

# API Reference

## Communicating the problem to the solver

```@docs
SDDPModel
Expectation
NestedAVaR
DefaultCutOracle
LevelOneCutOracle
@state
@states
@rhsnoise
@rhsnoises
setnoiseprobability!
@stageobjective
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
show
plotvaluefunction
```

## Read and write the model to disk

```@docs
    loadcuts!
    savemodel!
    loadmodel
```

## Extras for Experts
```@docs
AbstractCutOracle
storecut!
validcuts
```

```@docs
AbstractRiskMeasure
modifyprobability!
```
