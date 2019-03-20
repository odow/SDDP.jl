```@meta
CurrentModule = SDDP
```

# SDDP.jl

!!! warn
    `SDDP.jl` under went a major re-write to be compatible with JuMP v0.19 and
    Julia v1.0. The [Upgrading guide](@ref) has advice on how to upgrade your
    existing `SDDP.jl` models.

`SDDP.jl` is a package for solving large multistage convex stochastic
programming problems using *stochastic dual dynamic programming*. In this
manual, we're going to assume a reasonable amount of background knowledge about
stochastic optimization, the SDDP algorithm, Julia, and JuMP.

!!! info
    If you haven't used JuMP before, we recommend that you read the
    [JuMP documentation](http://www.juliaopt.org/JuMP.jl/latest/) and try
    building and solving JuMP models _before_ trying `SDDP.jl`.

## Installation

You can install `SDDP.jl` as follows:

```julia
import Pkg
Pkg.add("https://github.com/odow/SDDP.jl.git")
```

If you get an error like:
`ERROR: Unsatisfiable requirements detected for package MathOptFormat`,
run

```julia
import Pkg
Pkg.add("https://github.com/odow/MathOptFormat.jl.git")
Pkg.add("https://github.com/odow/SDDP.jl.git")
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

If you use `SDDP.jl`, we ask that you please cite the following
[paper](http://www.optimization-online.org/DB_FILE/2017/12/6388.pdf):
```
@article{dowson_sddp.jl,
	title = {{SDDP}.jl: a {Julia} package for stochastic dual dynamic programming},
	url = {http://www.optimization-online.org/DB_HTML/2017/12/6388.html},
	journal = {Optimization Online},
	author = {Dowson, Oscar and Kapelevich, Lea},
	year = {2017}
}
```

If you use the infinite horizon functionality, we ask that you please cite the
following [paper](http://www.optimization-online.org/DB_HTML/2018/11/6914.html):
```
@article{dowson_policy_graph,
	title = {The policy graph decomposition of multistage stochastic
      optimization problems},
	url = {http://www.optimization-online.org/DB_HTML/2018/11/6914.html},
	journal = {Optimization Online},
	author = {Dowson, Oscar},
	year = {2018}
}
```
