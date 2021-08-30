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
julia> ] add SDDP
```

## Tutorials

Once you've got `SDDP.jl` installed, you should read some tutorials, beginning
with [An introduction to SDDP.jl](@ref).

If you want a basic understanding of the algorithm behind SDDP.jl, start with
[Introductory theory](@ref).
## How-to guides

If you just want help on a specific topic, check out one of the how-to guides. A
good one to get started on is [Debug a model](@ref).

## Examples

`SDDP.jl` also contains a number of examples. A good one to get started on is
the [Hydro-thermal scheduling](@ref) problem. In particular, it shows how to
solve an infinite horizon problem.

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
	doi = {https://doi.org/10.1287/ijoc.2020.0987},
	year = {2021},
	volume = {33},
	issue = {1},
	pages = {27-33},
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_FILE/2017/12/6388.pdf).

If you use the infinite horizon functionality, we ask that you please cite the
following:
```
@article{dowson_policy_graph,
	title = {The policy graph decomposition of multistage stochastic optimization problems},
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

If you use the partially observable functionality, we ask that you please cite
the following:
```
@article{dowson_pomsp,
	title = {Partially observable multistage stochastic programming},
	doi = {https://doi.org/10.1016/j.orl.2020.06.005},
	journal = {Operations Research Letters},
	author = {Dowson, O. and Morton, D.P. and Pagnoncelli, B.K.},
	volume = {48},
	issue = {4},
	pages = {505-512},
	year = {2020}
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_HTML/2019/03/7141.html).

If you use the objective state functionality, we ask that you please cite the
following:
```
@article{downward_objective,
	title = {Stochastic dual dynamic programming with stagewise-dependent objective uncertainty},
	doi = {https://doi.org/10.1016/j.orl.2019.11.002},
	journal = {Operations Research Letters},
	author = {Downward, A. and Dowson, O. and Baucke, R.},
	volume = {48},
	issue = {1},
	pages = {33-39},
	year = {2020}
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_FILE/2018/02/6454.pdf).

If you use the entropic risk measure, we ask that you please cite the following:
```
@article{dowson_entropic,
	title = {Multistage stochastic programs with the entropic risk measure},
	journal = {Optimization Online},
	author = {Dowson, O. and Morton, D.P. and Pagnoncelli, B.K.},
	url = {http://www.optimization-online.org/DB_HTML/2020/08/7984.html},
}
```
Here is a [preprint](http://www.optimization-online.org/DB_HTML/2020/08/7984.html).
