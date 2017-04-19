#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################
using JuMP, Clp

@testset "Two stage" begin
    m = SDDPModel(
            sense  = :Min
            # ... or ...
            semse  = :Max

            stages = 2,

            markov_states = 2,
            # ... or ...
            markov_states = [2, 3],

            transition = Uniform
            # ... or ...
            transition = [0.5 0.5; 0.5 0.5] # time invariant
            # ... or ...
            transition = [ [0.4 0.4 0.2; 0.5 0.4 0.1] ] # time varying (stages-1) transitions
            # ... or ...
            transition = [[0.5, 0.5], [0.4 0.4 0.2; 0.5 0.4 0.1]] # stages transitions, first element is initial probability

            risk_measure = Expectation()
            # ... or ...
            risk_measure =

                        ) do subproblem, stage

        @state(subproblem, x >= 0, x0 == 5)

        @variables(subproblem, begin
            0 <= u <= 2
                 v
         end)

        @constraint(subproblem, x == x0 - u + v)

        @scenario(subproblem, w=[0.25, 0.5], v <= w)
        setscenarioprobability!(subproblem, [0.2, 0.8])

        stageobjective!(subproblem, stage * u)

    end

    solve(m,
        # solver
        solver::MathProgBase.AbstratMathProgSolver, # LP solver
        # logging
        log_level::Int,     # 0 = quiet, 1 = basic, 2 = iteration info
        log_output_file::String, # write log to file
        cut_output_file::String, # write cuts to file
        # convergence
        timeout::Float64,   # Time limit for algorithm (seconds)
        max_cuts::Int,      # Maximum cuts for algorithm
        rel_gap = BoundConvergence(
                iterations::Int,
                gap::Float64
            ),  # Relative bound decreases by less than rel_gap[2] for rel_gap[1] iterations
        policy_simulator = MonteCarloSimulator(
                frequency::Int,
                samples::Union{Int, AbstractVector{Int}}, # number of forward samples
                confidence_level::Float64,                # 0 - 1
                terminate::Bool     # terminate if statistical bound has converged
            ),
        # potential improvements
        cut_oracle::AbstractCutOracle = DefaultCutOracle(), # DeMatos()
        scenario_incrementation::AbstractVector{Int},
        init_expected_value_iterations::Int,
        sampling, # default, uniform, weighted
        parallel::Bool, # solve serial or asyncronous version
        solve
        )

end
#
# @testset "Rib Price Oracle" begin
#     @testset "Simple RHS" begin
#         m = SDDP.Subproblem()
#         ext = SDDP.ext(m)
#         @variable(m, x)
#         stageobjective!(m, x)
#         @test length(ext.valuefunctions) == 1
#     end
# end
