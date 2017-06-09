#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the Production Management example from FAST
    https://github.com/leopoldcambier/FAST/blob/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/production management multiple stages/
==#

using SDDP, JuMP, Clp, Base.Test

DEMAND = [2, 10]
H = 3
N = 2
C = [0.2, 0.7]
S = 2 + [0.33, 0.54]

m = SDDPModel(
                sense  = :Min,
                stages = H,
                solver = ClpSolver(),
                objective_bound = -50
                    ) do sp, t

    @state(sp,    x[i=1:N] >= 0, x0==0)
    @variable(sp, s[i=1:N] >= 0)
    @constraint(sp, s .<= x0)

    if t == 1
        @constraint(sp, sum(s) ==  0)
    else
        @noise(sp, d = DEMAND, sum(s) <= d)
    end

    stageobjective!(sp, dot(C, x) - dot(S, s))

end

status = solve(m, max_iterations = 10)

@test isapprox(getbound(m), -23.96, atol=1e-2)
