#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp
using Base.Test

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
    risk_measure = [Expectation(), Expectation(), EAVaR(lambda = 0.5, beta=0.5), Expectation()]
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
status = solve(m, iteration_limit = 30, print_level = 0)

#=
    These lines are commented out for testing

    To plot the value function against bonds when stock = 50
=#

# SDDP.plotvaluefunction(m, 3, 1, 50.0, 0.0:2.0:100, label1 = "Bonds")

#=
    To plot the value function against stocks and bonds
=#

# SDDP.plotvaluefunction(m, 3, 1, 0.0:2.0:100, 0.0:2.0:100, label1 = "Stocks", label2="Bonds")
