```@meta
CurrentModule = SDDP
```

# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.10.0 (November 19, 2024)

### Added

 - Added `root_node_risk_measure` kwarg to [`train`](@ref) (#804)

### Fixed

 - Fixed a bug with cut sharing in a graph with zero-probability arcs (#797)

### Other

 - Added a new tutorial [Example: inventory management](@ref) (#795)
 - Added a stochastic lead time example to docs (#800)

## v1.9.0 (October 17, 2024)

### Added

 - Added `write_only_selected_cuts` and `cut_selection` keyword arguments to
   [`write_cuts_to_file`](@ref) and [`read_cuts_from_file`](@ref) to skip
   potentially expensive operations (#781) (#784)
 - Added [`set_numerical_difficulty_callback`](@ref) to modify the subproblem on
   numerical difficulty (#790)

### Fixed

 - Fixed the tests to skip threading tests if running in serial (#770)
 - Fixed [`BanditDuality`](@ref) to handle the case where the standard deviation
   is `NaN` (#779)
 - Fixed an error when lagged state variables are encountered in `MSPFormat`
   (#786)
 - Fixed `publication_plot` with replications of different lengths (#788)
 - Fixed CTRL+C interrupting the code at unsafe points (#789)

### Other

 - Documentation improvements (#771) (#772)
 - Updated printing because of changes in JuMP (#773)

## v1.8.1 (August 5, 2024)

### Fixed

 - Fixed various issues with `SDDP.Threaded()` (#761)
 - Fixed a deprecation warning for sorting a dictionary (#763)

### Other

 - Updated copyright notices (#762)
 - Updated `.JuliaFormatter.toml` (#764)

## v1.8.0 (July 24, 2024)

### Added

 - Added `SDDP.Threaded()`, which is an experimental parallel scheme that
   supports solving problems using multiple threads. Some parts of SDDP.jl may
   not be thread-safe, and this can cause incorrect results, segfaults, or other
   errors. Please use with care and report any issues by opening a GitHub issue.
   (#758)

### Other

 - Documentation improvements and fixes (#747) (#759)

## v1.7.0 (June 4, 2024)

### Added

 - Added `sample_backward_noise_terms_with_state` for creating backward pass
   sampling schemes that depend on the current primal state. (#742) (Thanks
   @arthur-brigatto)

### Fixed

 - Fixed error message when `publication_plot` has non-finite data (#738)

### Other

 - Updated the logo constructor (#730)

## v1.6.7 (February 1, 2024)

### Fixed

 - Fixed non-constant state dimension in the `MSPFormat` reader (#695)
 - Fixed SimulatorSamplingScheme for deterministic nodes (#710)
 - Fixed line search in BFGS (#711)
 - Fixed handling of `NEARLY_FEASIBLE_POINT` status (#726)

### Other

 - Documentation improvements (#692) (#694) (#706) (#716) (#727)
 - Updated to StochOptFormat v1.0 (#705)
 - Added an experimental `OuterApproximation` algorithm (#709)
 - Updated `.gitignore` (#717)
 - Added code for MDP paper (#720) (#721)
 - Added Google analytics (#723)

## v1.6.6 (September 29, 2023)

### Other

 - Updated [Example: two-stage newsvendor](@ref) tutorial (#689)
 - Added a warning for people using [`SDDP.Statistical`](@ref) (#687)

## v1.6.5 (September 25, 2023)

### Fixed

 - Fixed duplicate nodes in [`MarkovianGraph`](@ref) (#681)

### Other

 - Updated tutorials (#677) (#678) (#682) (#683)
 - Fixed documentation preview (#679)

## v1.6.4 (September 23, 2023)

### Fixed

 - Fixed error for invalid `log_frequency` values (#665)
 - Fixed objective sense in [`deterministic_equivalent`](@ref) (#673)

### Other

 - Documentation updates (#658) (#666) (#671)
 - Switch to GitHub action for deploying docs (#668) (#670)
 - Update to Documenter@1 (#669)

## v1.6.3 (September 8, 2023)

### Fixed

 - Fixed default stopping rule with `iteration_limit` or `time_limit` set (#662)

### Other

 - Various documentation improvements (#651) (#657) (#659) (#660)

## v1.6.2 (August 24, 2023)

### Fixed

 - `MSPFormat` now detect and exploit stagewise independent lattices (#653)
 - Fixed `set_optimizer` for models read from file (#654)

### Other

 - Fixed typo in `pglib_opf.jl` (#647)
 - Fixed documentation build and added color (#652)

## v1.6.1 (July 20, 2023)

### Fixed

 - Fixed bugs in `MSPFormat` reader (#638) (#639)

### Other

 - Clarified `OutOfSampleMonteCarlo` docstring (#643)

## v1.6.0 (July 3, 2023)

### Added

 - Added [`RegularizedForwardPass`](@ref) (#624)
 - Added [`FirstStageStoppingRule`](@ref) (#634)

### Other

 - Removed an unbound type parameter (#632)
 - Fixed typo in docstring (#633)
 - Added [Here-and-now and hazard-decision](@ref) tutorial (#635)

## v1.5.1 (June 30, 2023)

This release contains a number of minor code changes, but it has a large impact
on the content that is printed to screen. In particular, we now log
periodically, instead of each iteration, and a "good" stopping rule is used as
the default if none are specified. Try using `SDDP.train(model)` to see the
difference.

### Other

 - Fixed various typos in the documentation (#617)
 - Fixed printing test after changes in JuMP (#618)
 - Set [`SimulationStoppingRule`](@ref) as the default stopping rule (#619)
 - Changed the default logging frequency. Pass `log_every_seconds = 0.0` to
   [`train`](@ref) to revert to the old behavior. (#620)
 - Added example usage with Distributions.jl (@slwu89) (#622)
 - Removed the numerical issue `@warn` (#627)
 - Improved the quality of docstrings (#630)

## v1.5.0 (May 14, 2023)

### Added

 - Added the ability to use a different model for the forward pass. This is a
   novel feature that lets you train better policies when the model is
   non-convex or does not have a well-defined dual. See the [Alternative forward models](@ref)
   tutorial in which we train convex and non-convex formulations of the optimal
   power flow problem. (#611)

### Other

 - Updated missing `changelog` entries (#608)
 - Removed global variables (#610)
 - Converted the `Options` struct to keyword arguments. This struct was a
   private implementation detail, but the change is breaking if you developed an
   extension to SDDP that touched these internals. (#612)
 - Fixed some typos (#613)

## v1.4.0 (May 8, 2023)

### Added

 - Added [`SDDP.SimulationStoppingRule`](@ref) (#598)
 - Added `sampling_scheme` argument to [`SDDP.write_to_file`](@ref) (#607)

### Fixed

 - Fixed parsing of some `MSPFormat` files (#602) (#604)
 - Fixed printing in header (#605)

## v1.3.0 (May 3, 2023)

### Added

 - Added experimental support for `SDDP.MSPFormat.read_from_file` (#593)

### Other

 - Updated to StochOptFormat v0.3 (#600)

## v1.2.1 (May 1, 2023)

### Fixed

 - Fixed `log_every_seconds` (#597)

## v1.2.0 (May 1, 2023)

### Added

 - Added [`SDDP.SimulatorSamplingScheme`](@ref) (#594)
 - Added `log_every_seconds` argument to [`SDDP.train`](@ref) (#595)

### Other

 - Tweaked how the log is printed (#588)
 - Updated to StochOptFormat v0.2 (#592)

## v1.1.4 (April 10, 2023)

### Fixed

 - Logs are now flushed every iteration (#584)

### Other

 - Added docstrings to various functions (#581)
 - Minor documentation updates (#580)
 - Clarified integrality documentation (#582)
 - Updated the README (#585)
 - Number of numerical issues is now printed to the log (#586)

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

 - Added the `node_name_parser` argument to [`SDDP.write_cuts_to_file`](@ref)
   and added the option to skip nodes in [`SDDP.read_cuts_from_file`](@ref)
   (#565)

## v1.0.0 (January 3, 2023)

Although we're bumping MAJOR version, this is a non-breaking release. Going
forward:

 - New features will bump the MINOR version
 - Bug fixes, maintenance, and documentation updates will bump the PATCH
   version
 - We will support only the Long Term Support (currently v1.6.7) and the latest
   patch (currently v1.8.4) releases of Julia. Updates to the LTS version will
   bump the MINOR version
 - Updates to the compat bounds of package dependencies will bump the PATCH
   version.

We do not intend any breaking changes to the public API, which would require a
new MAJOR release. The public API is everything defined in the documentation.
Anything not in the documentation is considered private and may change in any
PATCH release.

### Added

 - Added `num_nodes` argument to [`SDDP.UnicyclicGraph`](@ref) (#562)
 - Added support for passing an optimizer to [`SDDP.Asynchronous`](@ref) (#545)

### Other

 - Updated [Plotting tools](@ref) to use live plots (#563)
 - Added [vale](https://vale.sh) as a linter (#565)
 - Improved documentation for initializing a parallel scheme (#566)

## v0.4.9 (January 3, 2023)

### Added

 - Added [`SDDP.UnicyclicGraph`](@ref) (#556)

### Other

 - Added tutorial on Markov Decision Processes (#556)
 - Added two-stage newsvendor tutorial (#557)
 - Refactored the layout of the documentation (#554) (#555)
 - Updated copyright to 2023 (#558)
 - Fixed errors in the documentation (#561)

## v0.4.8 (December 19, 2022)

### Added

 - Added `terminate_on_cycle` option to [`SDDP.Historical`](@ref) (#549)
 - Added `include_last_node` option to [`SDDP.DefaultForwardPass`](@ref) (#547)

### Fixed

 - Reverted then fixed (#531) because it failed to account for problems with
   integer variables (#546) (#551)

## v0.4.7 (December 17, 2022)

### Added

 - Added `initial_node` support to `InSampleMonteCarlo` and
   `OutOfSampleMonteCarlo` (#535)

### Fixed

 - Rethrow `InterruptException` when solver is interrupted (#534)
 - Fixed numerical recovery when we need dual solutions (#531) (Thanks @bfpc)
 - Fixed re-using the `dashboard = true` option between solves (#538)
 - Fixed bug when no `@stageobjective` is set (now defaults to `0.0`) (#539)
 - Fixed errors thrown when invalid inputs are provided to `add_objective_state`
   (#540)

### Other

 - Drop support for Julia versions prior to 1.6 (#533)
 - Updated versions of dependencies (#522) (#533)
 - Switched to HiGHS in the documentation and tests (#533)
 - Added license headers (#519)
 - Fixed link in air conditioning example (#521) (Thanks @conema)
 - Clarified variable naming in deterministic equivalent (#525) (Thanks @lucasprocessi)
 - Added this change log (#536)
 - Cuts are now written to `model.cuts.json` when numerical instability is
   discovered. This can aid debugging because it allows to you reload the cuts
   as of the iteration that caused the numerical issue (#537)

## v0.4.6 (March 25, 2022)

### Other

 - Updated to JuMP v1.0 (#517)

## v0.4.5 (March 9, 2022)

### Fixed

 - Fixed issue with `set_silent` in a subproblem (#510)

### Other

 - Fixed many typos (#500) (#501) (#506) (#511) (Thanks @bfpc)
 - Update to JuMP v0.23 (#514)
 - Added auto-regressive tutorial (#507)

## v0.4.4 (December 11, 2021)

### Added

 - Added `BanditDuality` (#471)
 - Added benchmark scripts (#475) (#476) (#490)
 - `write_cuts_to_file` now saves visited states (#468)

### Fixed

 - Fixed `BoundStalling` in a deterministic policy (#470) (#474)
 - Fixed magnitude warning with `zero` coefficients (#483)

### Other

 - Improvements to LagrangianDuality (#481) (#482) (#487)
 - Improvements to `StrengthenedConicDuality` (#486)
 - Switch to functional form for the tests (#478)
 - Fixed typos (#472) (Thanks @vfdev-5)
 - Update to JuMP v0.22 (#498)

## v0.4.3 (August 31, 2021)

### Added

 - Added biobjective solver (#462)
 - Added `forward_pass_callback` (#466)

### Other

 - Update tutorials and documentation (#459) (#465)
 - Organize how paper materials are stored (#464)

## v0.4.2 (August 24, 2021)

### Fixed

 - Fixed a bug in Lagrangian duality (#457)

## v0.4.1 (August 23, 2021)

### Other

  - Minor changes to our implementation of `LagrangianDuality` (#454) (#455)

## v0.4.0 (August 17, 2021)

### Breaking

 - A large refactoring for how we handle stochastic integer programs. This added
   support for things like [`SDDP.ContinuousConicDuality`](@ref) and
   [`SDDP.LagrangianDuality`](@ref). It was breaking because we removed the
   `integrality_handler` argument to `PolicyGraph`. (#449) (#453)

### Other

 - Documentation improvements (#447) (#448) (#450)

## v0.3.17 (July 6, 2021)

### Added

 - Added [`SDDP.PSRSamplingScheme`](@ref) (#426)

### Other

 - Display more model attributes (#438)
 - Documentation improvements (#433) (#437) (#439)

## v0.3.16 (June 17, 2021)

### Added

 - Added [`SDDP.RiskAdjustedForwardPass`](@ref) (#413)
 - Allow [`SDDP.Historical`](@ref) to sample sequentially (#420)

### Other

 - Update risk measure docstrings (#418)

## v0.3.15 (June 1, 2021)

### Added

 - Added [`SDDP.StoppingChain`](@ref)

### Fixed

 - Fixed scoping bug in [`SDDP.@stageobjective`](@ref) (#407)
 - Fixed a bug when the initial point is infeasible (#411)
 - Set subproblems to silent by default (#409)

### Other

 - Add JuliaFormatter (#412)
 - Documentation improvements (#406) (#408)

## v0.3.14 (March 30, 2021)

### Fixed

 - Fixed `O(N^2)` behavior in `get_same_children` (#393)

## v0.3.13 (March 27, 2021)

### Fixed

 - Fixed bug in `print.jl`
 - Fixed compat of `Reexport` (#388)

## v0.3.12 (March 22, 2021)

### Added

 - Added problem statistics to header (#385) (#386)

### Fixed

 - Fixed subtypes in `visualization` (#384)

## v0.3.11 (March 22, 2021)

### Fixed

 - Fixed constructor in direct mode (#383)

### Other

 - Fix documentation (#379)

## v0.3.10 (February 23, 2021)

### Fixed

 - Fixed `seriescolor` in publication plot (#376)

## v0.3.9 (February 20, 2021)

### Added

 - Add option to simulate with different incoming state (#372)
 - Added warning for cuts with high dynamic range (#373)

### Fixed

 - Fixed `seriesalpha` in publication plot (#375)

## v0.3.8 (January 19, 2021)

### Other

 - Documentation improvements (#367) (#369) (#370)

## v0.3.7 (January 8, 2021)

### Other

 - Documentation improvements (#362) (#363) (#365) (#366)
 - Bump copyright (#364)

## v0.3.6 (December 17, 2020)

### Other

 - Fix typos (#358)
 - Collapse navigation bar in docs (#359)
 - Update `TagBot.yml` (#361)

## v0.3.5 (November 18, 2020)

### Other

 - Update citations (#348)
 - Switch to GitHub actions (#355)

## v0.3.4 (August 25, 2020)

### Added

 - Added non-uniform distributionally robust risk measure (#328)
 - Added numerical recovery functions (#330)
 - Added experimental StochOptFormat (#332) (#336) (#337) (#341) (#343) (#344)
 - Added entropic risk measure (#347)

### Other

 - Documentation improvements (#327) (#333) (#339) (#340)

## v0.3.3 (June 19, 2020)

### Added

 - Added asynchronous support for price and belief states (#325)
 - Added `ForwardPass` plug-in system (#320)

### Fixed

 - Fix check for probabilities in Markovian graph (#322)

## v0.3.2 (April 6, 2020)

### Added

 - Added `log_frequency` argument to [`SDDP.train`](@ref) (#307)

### Other

 - Improve error message in deterministic equivalent (#312)
 - Update to `RecipesBase` 1.0 (#313)

## v0.3.1 (February 26, 2020)

### Fixed

 - Fixed filename in `integrality_handlers.jl` (#304)

## v0.3.0 (February 20, 2020)

### Breaking

 - Breaking changes to update to JuMP v0.21 (#300).

## v0.2.4 (February 7, 2020)

### Added

 - Added a counter for the number of total subproblem solves (#301)

### Other

 - Update formatter (#298)
 - Added tests (#299)

## v0.2.3 (January 24, 2020)

### Added

 - Added support for convex risk measures (#294)

### Fixed

 - Fixed bug when subproblem is infeasible (#296)
 - Fixed bug in deterministic equivalent (#297)

### Other

 - Added example from IJOC paper (#293)

## v0.2.2 (January 10, 2020)

### Fixed

 - Fixed flakey time limit in tests (#291)

### Other

 - Removed MathOptFormat.jl (#289)
 - Update copyright (#290)

## v0.2.1 (December 19, 2019)

### Added

 - Added support for approximating a Markov lattice (#282) (#285)
 - Add tools for visualizing the value function (#272) (#286)
 - Write `.mof.json` files on error (#284)

### Other

 - Improve documentation (#281) (#283)
 - Update tests for Julia 1.3 (#287)

## v0.2.0 (December 16, 2019)

This version added the asynchronous parallel implementation with a few minor
breaking changes in how we iterated internally. It didn't break basic
user-facing models, only implementations that implemented some of the extension
features. It probably could have been a v1.1 release.

### Added

 - Added asynchronous parallel implementation (#277)
 - Added roll-out algorithm for cyclic graphs (#279)

### Other

 - Improved error messages in `PolicyGraph` (#271)
 - Added JuliaFormatter (#273) (#276)
 - Fixed compat bounds (#274) (#278)
 - Added documentation for simulating non-standard graphs (#280)

## v0.1.0 (October 17, 2019)

A complete rewrite of SDDP.jl based on the policy graph framework. This was
essentially a new package. It has minimal code in common with the previous
implementation.

Development started on September 28, 2018 in [Kokako.jl](https://github.com/odow/Kokako.jl),
and the code was merged into `SDDP.jl` on March 14, 2019.

The pull request [SDDP.jl#180](https://github.com/odow/SDDP.jl/pull/180) lists
the 29 issues that the rewrite closed.

## v0.0.1 (April 18, 2018)

Initial release. Development had been underway since January 22, 2016 in the
[StochDualDynamicProgram.jl](https://github.com/odow/StochDualDynamicProgram.jl)
repository. The last development commit there was April 5, 2017. Work then
continued in this repository for a year before the first tagged release.
