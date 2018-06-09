```@meta
CurrentModule = SDDP
```
# SDDP.jl Documentation

SDDP.jl is a package for solving large multistage convex stochastic optimization
problems using _stochastic dual dynamic programming_. In this manual, we're
going to assume a reasonable amount of background knowledge about stochastic
optimization, the SDDP algorithm, Julia, and JuMP.

!!! note
    If you don't have that background, you may want to brush up on some
    [Readings](@ref).

## Getting started

This package is unregistered so you will need to `Pkg.clone` it as follows:

```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```

If you want to use the parallel features of SDDP.jl, you should start Julia with
some worker processes (`julia -p N`), or add by running `julia> addprocs(N)` in
a running Julia session.

Once you've got SDDP.jl installed, you should read some tutorials, beginning with
[Tutorial One: first steps](@ref).

## Citing SDDP.jl

If you use SDDP.jl, we ask that you please cite the following [paper](http://www.optimization-online.org/DB_FILE/2017/12/6388.pdf):
```
@article{dowson_sddp.jl,
	title = {{SDDP}.jl: a {Julia} package for stochastic dual dynamic programming},
	url = {http://www.optimization-online.org/DB_HTML/2017/12/6388.html},
	journal = {Optimization Online},
	author = {Dowson, Oscar and Kapelevich, Lea},
	year = {2017}
}
```

## FAQ

**Q.** How do I make the constraint coefficients random?

**A.** Due to the design of JuMP, it's difficult to efficiently modify constraint
coefficients. Therefore, you can only vary the right hand-side of a constraint
using the `@rhsnoise` macro.

As a work around, we suggest you either reformulate the model so the uncertainty
appears in the RHS, or model the uncertainty as a Markov process.
[Tutorial Four: Markovian policy graphs](@ref) explains how to implement this.
You might also want to take a look at the [asset management example](https://github.com/odow/SDDP.jl/blob/master/examples/AssetManagement/asset_management.jl)
to see an example of this. Make sure you keep in mind that a new value function
is built at each Markov state which increases the computation time and memory
requirements.
