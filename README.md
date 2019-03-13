# SDDP

| **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:--------------------:|:----------------:|
| [![][docs-latest-img]][docs-latest-url] | [![Build Status][build-img]][build-url] | [![Codecov branch][codecov-img]][codecov-url]

## Documentation

You can find the documentation at https://odow.github.io/SDDP.jl/latest/.

**If you are struggling to figure out how to use something, raise a Github issue!**

## Other Packages

`SDDP.jl` isn't the only Julia package for solving multi-stage stochastic programs.
You may want to checkout [StructDualDynProg.jl](https://github.com/blegat/StructDualDynProg.jl)
or [StochDynamicProgramming.jl](https://github.com/JuliaOpt/StochDynamicProgramming.jl)
to see if they better suit your needs.

## Citing SDDP.jl

If you use SDDP.jl, we ask that you please cite the following [paper](http://www.optimization-online.org/DB_FILE/2017/12/6388.pdf):
```
@article{dowson_sddp.jl,
	title = {{SDDP}.jl: a {Julia} package for {Stochastic} {Dual} {Dynamic} {Programming}},
	url = {http://www.optimization-online.org/DB_HTML/2017/12/6388.html},
	journal = {Optimization Online},
	author = {Dowson, Oscar and Kapelevich, Lea},
	year = {2017}
}
```

[build-img]: https://travis-ci.org/odow/SDDP.jl.svg?branch=master
[build-url]: https://travis-ci.org/odow/SDDP.jl

[codecov-img]: https://codecov.io/github/odow/SDDP.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/odow/SDDP.jl?branch=master

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://odow.github.io/SDDP.jl/latest/
