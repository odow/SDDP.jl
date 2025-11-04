# # Advanced III: integrality

# The fundamental reason why SDDP works is convexity. In the classical
# algorithm, this ruled out the use of integer variables. However, various
# extensions to the algorithm have been proposed, and these can be run by
# `SDDP.jl` using the `integrality_handler` keyword argument to
# [`SDDP.PolicyGraph`](@ref).

# ## Continuous relaxation

# If your model includes binary or integer variables (e.g.,
# [`air_conditioning.jl`](https://github.com/odow/SDDP.jl/blob/master/examples/air_conditioning.jl)),
# passing `integrality_handler = SDDP.ContinuousRelaxation()` will make
# `SDDP.jl` contruct a (sub-optimal) policy using the continuous relaxation of
# the problem. But, when you simulate this policy, `SDDP.jl` will solve the
# original mixed-integer problem.

# ## SDDiP

# If you pass `integrality_handler = SDDP.SDDiP()`, SDDP.jl will use the
# stochastic dual dynamic integer programming method of Zou, Ahmed, and Sun. See
# [`SDDP.SDDiP`](@ref) for more.
