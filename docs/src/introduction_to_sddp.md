```@meta
CurrentModule = SDDP
```

## Introduction to SDDP

To start, lets discuss the types of problems SDDP.jl can solve, since it has a
few features that are non-standard, and it is missing some features that are
standard.

SDDP.jl can solve multi-stage convex stochastic optimizations problems with
 - a finite discrete number of states;
 - continuous state and control variables;
 - Hazard-Decision (Wait-and-See) uncertainty realization;
 - stagewise independent uncertainty in the RHS of the constraints that is
    drawn from a finite discrete distribution;
 - stagewise independent uncertainty in the objective function that is
    drawn from a finite discrete distribution;
 - a markov chain for temporal dependence. The markov chain forms a directed,
    acyclic, feed-forward graph with a finite (and at least one) number of
    markov states in each stage.

!!! note
    Stagewise independent uncertainty in the constraint coefficients is **not**
    supported. You should reformulate the problem, or model the uncertainty as
    a markov chain.
