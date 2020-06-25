```@meta
CurrentModule = SDDP
```

# SDDP.jl

!!! warning
    `SDDP.jl` under went a major re-write to be compatible with JuMP v0.19 and
    Julia v1.0. The [Upgrade from the old SDDP.jl](@ref) guide has advice on
    how to upgrade your existing `SDDP.jl` models.

`SDDP.jl` is a package for solving large multistage convex stochastic
programming problems using *stochastic dual dynamic programming*. In this
manual, we're going to assume a reasonable amount of background knowledge about
stochastic optimization, the SDDP algorithm, Julia, and JuMP.

!!! tip
    If you haven't used JuMP before, we recommend that you read the
    [JuMP documentation](http://www.juliaopt.org/JuMP.jl/latest/) and try
    building and solving JuMP models _before_ trying `SDDP.jl`.

## Installation

You can install `SDDP.jl` as follows:

```julia
julia> ] add https://github.com/odow/SDDP.jl.git
```

### Want the old version?

Still using Julia 0.6 and things broke when you went `Pkg.update()`? Run
```julia
julia> Pkg.checkout("SDDP", "release-v0")
```

## Tutorials

Once you've got `SDDP.jl` installed, you should read some tutorials, beginning
with [Basic I: first steps](@ref).

## Citing `SDDP.jl`

If you use `SDDP.jl`, we ask that you please cite the following:
```
@article{dowson_sddp.jl,
	title = {{SDDP}.jl: a {Julia} package for stochastic dual dynamic programming},
	journal = {INFORMS Journal on Computing},
	author = {Dowson, O. and Kapelevich, L.},
	note = {in press},
	year = {2020}
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_FILE/2017/12/6388.pdf).

If you use the infinite horizon functionality, we ask that you please cite the
following:
```
@article{dowson_policy_graph,
	title = {The policy graph decomposition of multistage stochastic
      optimization problems},
	doi = {https://doi.org/10.1002/net.21932},
	journal = {Networks},
	author = {Dowson, O.},
	volume = {76},
	issue = {1},
	pages = {3-23},
	year = {2020}
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_HTML/2018/11/6914.html).

If you use the partially observable functionality, we ask that you please cite the
following:
```
@article{dowson_pomsp,
	title = {Partially observable multistage stochastic programming},
	doi = {https://doi.org/10.1016/j.orl.2020.06.005},
	journal = {Operations Research Letters},
	author = {Dowson, O., Morton, D.P., and Pagnoncelli, B.K.},
	volume = {48},
	issue = {4},
	pages = {505-512},
	year = {2020}
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_HTML/2019/03/7141.html).
