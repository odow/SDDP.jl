```@meta
CurrentModule = SDDP
```

# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.3 (April 2, 2023)

### Other

 - Fixed typo in [Example: deterministic to stochastic](@ref) tutorial (#578)
 - Fixed typo in documentation of [`SDDP.simulate`](@ref) (#577)

## v1.1.2 (March 18, 2023)

### Other

 - Added [Example: deterministic to stochastic](@ref) tutorial (#572)

## v1.1.1 (March 16, 2023)

### Other

 - Fixed email in `Project.toml`
 - Added notebook to documentation tutorials (#571)

## v1.1.0 (January 12, 2023)

### Added

 * Added the `node_name_parser` argument to [`SDDP.write_cuts_to_file`](@ref)
   and added the option to skip nodes in [`SDDP.read_cuts_from_file`](@ref)
   (#565)

## v1.0.0 (January 3, 2023)

Although we're bumping MAJOR version, this is a non-breaking release. Going
forward:

 * New features will bump the MINOR version
 * Bug fixes, maintenance, and documentation updates will bump the PATCH
   version
 * We will support only the Long Term Support (currently v1.6.7) and the latest
   patch (currently v1.8.4) releases of Julia. Updates to the LTS version will
   bump the MINOR version
 * Updates to the compat bounds of package dependencies will bump the PATCH
   version.

We do not intend any breaking changes to the public API, which would require a
new MAJOR release. The public API is everything defined in the documentation.
Anything not in the documentation is considered private and may change in any
PATCH release.

### Added

 * Added `num_nodes` argument to [`SDDP.UnicyclicGraph`](@ref) (#562)
 * Added support for passing an optimizer to [`SDDP.Asynchronous`](@ref) (#545)

### Other

 * Updated [Plotting tools](@ref) to use live plots (#563)
 * Added [vale](https://vale.sh) as a linter (#565)
 * Improved documentation for initializing a parallel scheme (#566)

## v0.4.9 (January 3, 2023)

### Added

 * Added [`SDDP.UnicyclicGraph`](@ref) (#556)

### Other

 * Added tutorial on Markov Decision Processes (#556)
 * Added two-stage newsvendor tutorial (#557)
 * Refactored the layout of the documentation (#554) (#555)
 * Updated copyright to 2023 (#558)
 * Fixed errors in the documentation (#561)

## v0.4.8 (December 19, 2022)

### Added

 * Added `terminate_on_cycle` option to [`SDDP.Historical`](@ref) (#549)
 * Added `include_last_node` option to [`SDDP.DefaultForwardPass`](@ref) (#547)

### Fixed

 * Reverted then fixed (#531) because it failed to account for problems with
   integer variables (#546) (#551)

## v0.4.7 (December 17, 2022)

### Added

 * Added `initial_node` support to `InSampleMonteCarlo` and
   `OutOfSampleMonteCarlo` (#535)

### Fixed

 * Rethrow `InterruptException` when solver is interrupted (#534)
 * Fixed numerical recovery when we need dual solutions (#531) (Thanks @bfpc)
 * Fixed re-using the `dashboard = true` option between solves (#538)
 * Fixed bug when no `@stageobjective` is set (now defaults to `0.0`) (#539)
 * Fixed errors thrown when invalid inputs are provided to `add_objective_state`
   (#540)

### Other

 * Drop support for Julia versions prior to 1.6 (#533)
 * Updated versions of dependencies (#522) (#533)
 * Switched to HiGHS in the documentation and tests (#533)
 * Added license headers (#519)
 * Fixed link in air conditioning example (#521) (Thanks @conema)
 * Clarified variable naming in deterministic equivalent (#525) (Thanks @lucasprocessi)
 * Added this change log (#536)
 * Cuts are now written to `model.cuts.json` when numerical instability is
   discovered. This can aid debugging because it allows to you reload the cuts
   as of the iteration that caused the numerical issue (#537)

## v0.4.6 (March 25, 2022)

### Other

 * Updated to JuMP v1.0 (#517)

## v0.4.5 (March 9, 2022)

### Fixed

 * Fixed issue with `set_silent` in a subproblem (#510)

### Other

 * Fixed many typos (#500) (#501) (#506) (#511) (Thanks @bfpc)
 * Update to JuMP v0.23 (#514)
 * Added auto-regressive tutorial (#507)

## v0.4.4 (December 11, 2021)

### Added

 * Added `BanditDuality` (#471)
 * Added benchmark scripts (#475) (#476) (#490)
 * `write_cuts_to_file` now saves visited states (#468)

### Fixed

 * Fixed `BoundStalling` in a deterministic policy (#470) (#474)
 * Fixed magnitude warning with `zero` coefficients (#483)

### Other

 * Improvements to LagrangianDuality (#481) (#482) (#487)
 * Improvements to `StrengthenedConicDuality` (#486)
 * Switch to functional form for the tests (#478)
 * Fixed typos (#472) (Thanks @vfdev-5)
 * Update to JuMP v0.22 (#498)

## v0.4.3 (August 31, 2021)

### Added

 * Added biobjective solver (#462)
 * Added `forward_pass_callback` (#466)

### Other

 * Update tutorials and documentation (#459) (#465)
 * Organize how paper materials are stored (#464)

## v0.4.2 (August 24, 2021)

### Fixed

 * Fixed a bug in Lagrangian duality (#457)

## v0.4.1 (August 23, 2021)

### Other

  * Minor changes to our implementation of `LagrangianDuality` (#454) (#455)

## v0.4.0 (August 17, 2021)

A large refactoring for how we handle stochastic integer programs. This added
support for things like `ContinuousConicDuality` and `LagrangianDuality`. It was
breaking because we removed the `integrality_handler` argument to `PolicyGraph`.

## v0.3.0 (October 17, 2019)

Breaking changes to update to JuMP v0.21 (#300).

## v0.2.0 (December 16, 2019)

This version added the asynchronous parallel implementation with a few minor
breaking changes in how we iterated internally. It didn't break basic
user-facing models, only implementations that implemented some of the extension
features.

## v0.1.0 (October 17, 2019)

A complete rewrite of SDDP.jl based on the policy graph framework. This was
essentially a new package. It has minimal code in common with the previous
implementation.

## v0.0.1 (April 18, 2018)

Initial release. Development had been underway since January 22, 2016 in the
[StochDualDynamicProgram.jl](https://github.com/odow/StochDualDynamicProgram.jl)
repository. The last development commit there was April 5, 2017. Work then
continued in this repository for a year before the first tagged release.
