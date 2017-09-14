# SDDP

| **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:--------------------:|:----------------:|
| [![][docs-latest-img]][docs-latest-url] | [![Build Status][build-img]][build-url] | [![Codecov branch][codecov-img]][codecov-url]

## Installation
This package is unregistered so you will need to `Pkg.clone` it as follows:
```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```
## Development Notes

This package is under development and some features may be changed or added.
Most of the time, this will just be a renaming of minor parts of the code.

For a time line of features, see [NEWS.md](https://github.com/odow/SDDP.jl/NEWS.md).
Items that require modification of existing codes are prefixed with `!!`.

## Documentation

The documentation is still very incomplete, and the internals of the library
need a tidy and a refactor, however the user-facing API from the examples should
be stable enough to use.

However, you can find some documentation at https://odow.github.io/SDDP.jl/build/index.html

Or have a read of the [draft tutorial/paper](https://github.com/odow/SDDP.jl/raw/master/draft_paper.pdf)
which goes into a bit more depth.

## SDDiP

[@lkapelevich](https://github.com/lkapelevich) wrote an extension for SDDP.jl to
solve multi-stage stochastic programs with binary state variables. Check it out
at https://github.com/lkapelevich/SDDiP.jl!

[build-img]: https://travis-ci.org/odow/SDDP.jl.svg?branch=master
[build-url]: https://travis-ci.org/odow/SDDP.jl

[codecov-img]: https://codecov.io/github/odow/SDDP.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/odow/SDDP.jl?branch=master

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://odow.github.io/SDDP.jl/build/index.html
