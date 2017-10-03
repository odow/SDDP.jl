#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
Modified version of the Asset Management problem taken from

    J. R. Birge,  F. Louveaux,  Introduction to Stochastic Programming,
    Springer Series in Operations Research and Financial Engineering,
    Springer New York, New York, NY, 2011
==#

using SDDP, JuMP, Clp, Base.Test

ws = [1.25, 1.06]
wb = [1.14, 1.12]
Phi = [-1, 5]
Psi = [0.02, 0.0]

m = SDDPModel(
    sense = :Min,
    stages = 4,
    objective_bound = -1000.0,
    solver = ClpSolver(),
    markov_transition = Array{Float64, 2}[
        [1.0]',
        [0.5 0.5],
        [0.5 0.5; 0.5 0.5],
        [0.5 0.5; 0.5 0.5]
    ],
    risk_measure = [Expectation(), Expectation(), NestedAVaR(lambda = 0.5, beta=0.5), Expectation()]
                            ) do sp, t, i
    @state(sp, xs >= 0, xsbar==0)
    @state(sp, xb >= 0, xbbar==0)
    if t == 1
        @constraint(sp, xs + xb == 55 + xsbar + xbbar)
        @stageobjective(sp, 0)
    elseif t == 2 || t == 3
        @rhsnoise(sp, phi=Phi, ws[i] * xsbar +
            wb[i] * xbbar + phi == xs + xb)
        @stageobjective(sp, psi = Psi, -psi * xs)
        setnoiseprobability!(sp, [0.6, 0.4])
    else
        @variable(sp, u  >= 0)
        @variable(sp, v >= 0)
        @constraint(sp, ws[i] * xsbar + wb[i] * xbbar +
            u - v == 80)
        @stageobjective(sp, 4u - v)
    end
end

srand(111)
@time status = solve(m,
    max_iterations = 100,
    simulation = MonteCarloSimulation(
        frequency = 5,
        min  = 100,
        step = 100,
        max  = 500,
        termination = false
    ),
    bound_convergence = BoundConvergence(
        iterations = 5,
        rtol       = 0.0,
        atol       = 0.001
    ),
    log_file="asset.log"
)
rm("asset.log")
@test status == :bound_convergence
@test isapprox(getbound(m), -1.278, atol=1e-3)

# results = simulate(m, 100, [:xs, :xb])
#
# p = SDDP.newplot()
# SDDP.addplot!(p, 1:100, 1:4, (i, t)->round(results[i][:stageobjective][t], 2),
#     title="Objective", ylabel="Profit (\$)", cumulative=true)
# SDDP.addplot!(p, 1:100, 1:4, (i, t)->round(results[i][:xs][t], 2),
#     title="Stocks")
# SDDP.addplot!(p, 1:100, 1:4, (i, t)->round(results[i][:xb][t], 2),
#     title="Bonds")
# SDDP.show(p)
