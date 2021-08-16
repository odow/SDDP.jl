# Integrality

There's nothing special about binary and integer variables in SDDP.jl. Use them
at will!

If you want finer control, you can pass an [`SDDP.AbstractDualityHandler`](@ref)
to the `duality_handler` argument of [`SDDP.train`](@ref).

See the [Duality handlers](@ref) section for the list of options you can pass.

!!! info
    Wondering where "SDDiP" is? SDDiP is vanilla SDDP, except that we use
    [`SDDP.LagrangianDuality`](@ref) to compute the dual variables.
