```@meta
CurrentModule = Kokako
```

# Kokako.jl

!!! note
    SDDP.jl is currently undergoing a re-write in this repository under the name
    Kokako.jl. Once completed, this package will be renamed back to SDDP.jl.


Kokako.jl is a package for solving large multistage convex stochastic
programming problems using *stochastic dual dynamic programming*. In this
manual, we're going to assume a reasonable amount of background knowledge about
stochastic optimization, the SDDP algorithm, Julia, and JuMP.

## Installation

You can install `Kokako.jl` as follows:

```julia
import Pkg
Pkg.add("https://github.com/odow/Kokako.jl.git")
```

Once you've got Kokako installed, you should read some tutorials, beginning
with [Tutorial One: first steps](@ref).

## Citing SDDP.jl

If you use SDDP.jl, we ask that you please cite the following
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
