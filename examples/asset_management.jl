#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
The Asset Management problem taken from

    J. R. Birge,  F. Louveaux,  Introduction to Stochastic Programming,
    Springer Series in Operations Research and Financial Engineering,
    Springer New York, New York, NY, 2011
==#

using SDDP, JuMP, Clp
using Base.Test

rstock = [1.25, 1.06]
rbonds = [1.14, 1.12]

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
    ]
) do sp, t, i

    @state(sp, stock_out >= 0, stock_in==0)
    @state(sp, bonds_out >= 0, bonds_in==0)
    if t == 1
        @constraint(sp, stock_out + bonds_out == 55)
        stageobjective!(sp, 0)
    elseif t > 1 && t < 4
        @constraint(sp, rstock[i] * stock_in + rbonds[i] * bonds_in == stock_out + bonds_out)
        stageobjective!(sp, 0)
    else
        @variable(sp, over  >= 0)
        @variable(sp, short >= 0)
        @constraint(sp, rstock[i] * stock_in + rbonds[i] * bonds_in - over + short == 80)
        stageobjective!(sp, -over + 4*short)
    end
end

srand(111)
@time solve(m, max_iterations = 25,
print_level = 0)
@test isapprox(getbound(m), 1.514, atol=1e-4)
