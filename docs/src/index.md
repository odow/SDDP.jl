```@meta
CurrentModule = SDDP
```

```@raw html
<img src="assets/logo_without_text.svg" alt="logo" width="150px"/>
```

# Introduction

[![Build Status](https://github.com/odow/SDDP.jl/workflows/CI/badge.svg?branch=master)](https://github.com/odow/SDDP.jl/actions?query=workflow%3ACI)
[![code coverage](https://codecov.io/gh/odow/SDDP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/odow/SDDP.jl)

Welcome to [SDDP.jl](https://github.com/odow/SDDP.jl), a package for solving large
convex multistage stochastic programming problems using stochastic dual dynamic
programming.

SDDP.jl is built on [JuMP](https://jump.dev), so it supports a number of
open-source and commercial solvers, making it a powerful and flexible tool for
stochastic optimization.

The implementation of the stochastic dual dynamic programming algorithm in
SDDP.jl is state of the art, and it includes support for a number of advanced
features not commonly found in other implementations. This includes support for:

 * infinite horizon problems
 * convex risk measures
 * mixed-integer state and control variables
 * partially observable stochastic processes.

## Installation

Install `SDDP.jl` as follows:

```julia
julia> import Pkg

julia> Pkg.add("SDDP")
```

## License

`SDDP.jl` is licensed under the [MPL 2.0 license](https://github.com/odow/SDDP.jl/blob/master/LICENSE.md).

## Resources for getting started

There are a few ways to get started with SDDP.jl:

 * Become familiar with JuMP by reading the [JuMP documentation](http://jump.dev/JuMP.jl/stable/)
 * Read the introductory tutorial [An introduction to SDDP.jl](@ref)
 * Browse some of the examples, such as [Example: deterministic to stochastic](@ref)

## Getting help

If you need help, please [open a GitHub issue](https://github.com/odow/SDDP.jl/issues/new).

## How the documentation is structured

Having a high-level overview of how this documentation is structured will help
you know where to look for certain things.

* **Tutorials** contains step-by-step explanations of how to use SDDP.jl. Once
  you've got `SDDP.jl` installed, start by reading [An introduction to SDDP.jl](@ref).

* **Guides** contains "how-to" snippets that demonstrate specific topics within
  SDDP.jl. A good one to get started on is [Debug a model](@ref).

* **Explanation** contains step-by-step explanations of the theory and
  algorithms that underpin SDDP.jl. If you want a basic understanding of the
  algorithm behind SDDP.jl, start with [Introductory theory](@ref).

* **Examples** contain worked examples of various problems solved using SDDP.jl.
  A good one to get started on is the [Hydro-thermal scheduling](@ref) problem.
  In particular, it shows how to solve an infinite horizon problem.

* The **API Reference** contains a complete list of the functions you can use in
  SDDP.jl. Look here if you want to know how to use a particular function.

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
	title = {Incorporating convex risk measures into multistage stochastic programming algorithms},
	doi = {https://doi.org/10.1007/s10479-022-04977-w},
	journal = {Annals of Operations Research},
	author = {Dowson, O. and Morton, D.P. and Pagnoncelli, B.K.},
	year = {2022},
}
```
Here is an earlier [preprint](http://www.optimization-online.org/DB_HTML/2020/08/7984.html).
