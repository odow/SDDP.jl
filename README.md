# SDDP

| **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:--------------------:|:----------------:|
| [![][docs-latest-img]][docs-latest-url] | [![Build Status][build-img]][build-url] | [![Codecov branch][codecov-img]][codecov-url]

## Documentation

You can find the documentation at https://odow.github.io/SDDP.jl/latest/.

**If you are struggling to figure out how to use something, raise a Github issue!**

## Examples

We need your examples! We're trying to collate a large array of examples to test the
correctness (and later, performance) of the package. Either make a PR or go to the
examples folder and click [`Upload Files`](https://github.com/odow/SDDP.jl/upload/master/examples) and Github will walk you through the process.
Bonus points for models where you know the optimal first stage objective value.

## Bugs

We need your bug reports! We've only stressed a few code paths on real-world models.
If you run into any problems, [file an issue here](https://github.com/odow/SDDP.jl/issues/new).

## Other Packages

`SDDP.jl` isn't the only Julia package for solving multi-stage stochastic programs.
You may want to checkout [StructDualDynProg.jl](https://github.com/blegat/StructDualDynProg.jl)
or [StochDynamicProgramming.jl](https://github.com/JuliaOpt/StochDynamicProgramming.jl)
to see if they better suit your needs.

## SDDiP

[@lkapelevich](https://github.com/lkapelevich) wrote an extension for SDDP.jl to
solve multi-stage stochastic programs with binary state variables. Check it out
at https://github.com/lkapelevich/SDDiP.jl!

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
