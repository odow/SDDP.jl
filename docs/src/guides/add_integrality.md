# Integrality

There's nothing special about binary and integer variables in SDDP.jl. Use them
at will!

If you want finer control, you can pass an [`SDDP.AbstractDualityHandler`](@ref)
to the `duality_handler` argument of [`SDDP.train`](@ref).

See the [Duality handlers](@ref) section for the list of options you can pass.

!!! info
    Wondering where "SDDiP" is? SDDiP is vanilla SDDP, except that we use
    [`SDDP.LagrangianDuality`](@ref) to compute the dual variables.

## Breaking changes in SDDP.jl v0.4

SDDP.jl v0.4.0 introduced a number of breaking changes in how we deal with
binary and integer variables.

### Breaking changes

 * We have renamed `SDDiP` to [`SDDP.LagrangianDuality`](@ref).
 * We have renamed `ContinuousRelaxation` to [`SDDP.ContinuousConicDuality`](@ref).
 * Instead of passing the argument to [`SDDP.PolicyGraph`](@ref), you now pass
   it to [`SDDP.train`](@ref), e.g.,
   `SDDP.train(model; duality_handler = SDDP.LagrangianDuality())`
 * We no longer turn continuous and integer states into a binary expansion. If
   you want to binarize your states, do it manually.

### Why did we do this?

SDDiP (the algorithm presented in the paper) is really two parts:

 1. If you have an integer program, you can compute the dual of the fishing
    constraint using Lagrangian duality; and
 2. If you have pure binary state variables, then cuts constructed from the
    Lagrangian duals result in an optimal policy.

However, these two points are quite independent. If you have integer or
continuous state variables, you can still use Lagrangian duality!

The new system is more flexible because the duality handler is a property of the
solution process, not the model. This allows us to use Lagrangian duality to
solve any dual problem, and it leaves the decision of binarizing the state
variables up to the user. (Hint: we don't think you should do it!)

### Other additions

We also added support for strengthened Benders cuts, which we call
[`SDDP.StrengthenedConicDuality`](@ref).

### Future plans

We have a number of future plans in the works, including better Lagrangian
solution methods and better ways of integrating the different types of duality
handlers (e.g., start with [`SDDP.ContinuousConicDuality`](@ref), then shift to
[`SDDP.StrengthenedConicDuality`](@ref), then [`SDDP.LagrangianDuality`](@ref)).

If these sorts of things interest you, the code is now much more hackable, so
please reach out or read [Issue #246](https://github.com/odow/SDDP.jl/issues/246).

Alternatively, if you have interesting examples using SDDiP that you find are
too slow, please send me the examples so we can use them as benchmarks in future
improvements.
