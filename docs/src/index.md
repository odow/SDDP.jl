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

!!! note
    You can find the old, terribly incomplete documentation at [Old Manual](@ref).

## Getting started

This package is unregistered so you will need to `Pkg.clone` it as follows:

```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```

If you want to use the parallel features of SDDP.jl, you should start Julia with
some worker processes (`julia -p N`), or add by running `julia> addprocs(N)` in
a running Julia session.

Once you've got SDDP.jl installed, you should read some tutorials, beginnng with
[First steps](@ref).

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
