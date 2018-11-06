# Kōkako.jl
*\[kor-kah-co\]*

| **Build Status** | **Coverage** |
|:--------------------:|:----------------:|
| [![Build Status][build-img]][build-url] | [![Codecov branch][codecov-img]][codecov-url]


[build-img]: https://travis-ci.com/odow/Kokako.jl.svg?branch=master
[build-url]: https://travis-ci.com/odow/Kokako.jl

[codecov-img]: https://codecov.io/github/odow/Kokako.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/odow/Kokako.jl?branch=master

__Note: this package is under active development. Things may change. If you currently use 
[SDDP.jl](https://github.com/odow/SDDP.jl), we recommend you keep using it.__

Kōkako.jl is a Julia package for solving multistage stochastic optimization
problems. It is very similar to the package [SDDP.jl](https://github.com/odow/SDDP.jl),
but has been written for Julia 1.0 and JuMP 0.19.

<p align="center">
  <img width="20%" src="https://upload.wikimedia.org/wikipedia/commons/a/a2/K%C5%8Dkako.jpg">
</p>

### Citing Kokako.jl

If you use Kokako.jl, we ask that you please cite the following paper:
```
@article{dowson_policy_graph,
	title = {The policy graph decomposition of multistage stochastic optimization problems},
	url = {},
	journal = {Optimization Online},
	author = {Dowson, Oscar},
	year = {2018}
}
```
Photo credit from [Wikipedia](https://en.wikipedia.org/wiki/K%c5%8dkako#/media/File:K%C5%8Dkako.jpg).
