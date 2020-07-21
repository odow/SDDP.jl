```@meta
CurrentModule = SDDP
```

# SDDP.jl

[`SDDP.jl`](https://github.com/odow/SDDP.jl) is a package for solving large
multistage convex stochastic programming problems using *stochastic dual dynamic
programming*. In this manual, we're going to assume a reasonable amount of
background knowledge about stochastic optimization, the SDDP algorithm, Julia,
and JuMP.

!!! tip
    If you haven't used JuMP before, we recommend that you read the
    [JuMP documentation](http://www.juliaopt.org/JuMP.jl/latest/) and try
    building and solving JuMP models _before_ trying `SDDP.jl`.

## Installation

You can install `SDDP.jl` as follows:

```julia
julia> ] add https://github.com/odow/SDDP.jl.git
```

## Tutorials

Once you've got `SDDP.jl` installed, you should read some tutorials, beginning
with [Basic I: first steps](@ref).

## How-to guides

If you just want help on a specific topic, check out one of the how-to guides. A
good one to get started on is [Debug a model](@ref).

## Examples

`SDDP.jl` also contains a number of examples. A good one to get started on is
the [Hydro-thermal scheduling](@ref) problem. In particular, it shows how to
solve an infinite horizon problem.

There is also a whole folder of coded examples in the [examples directory](https://github.com/odow/SDDP.jl/tree/master/examples)
of the [Github page](https://github.com/odow/SDDP.jl).

## API Reference

If you just want help on a specific function, see the [API Reference](@ref api_reference_list)
page.

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
