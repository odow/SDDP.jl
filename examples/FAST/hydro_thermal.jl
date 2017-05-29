#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the Hydro-thermal example from FAST
    https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/examples/hydro%20thermal
==#

using SDDP, JuMP, Clp, Base.Test

RAINFALL = [2, 10]

m = SDDPModel(
                sense  = :Min,
                stages = 2,
                solver = ClpSolver(),
                objective_bound = 0.0
                    ) do sp, t
    @state(sp, 0 <= x <= 8, x0==0)
    @variable(sp, y >= 0)
    @variable(sp, p >= 0)
    @constraint(sp, p + y >= 6)
    if t == 1
        @constraint(sp, x <= x0 + mean(RAINFALL) - y)
    else
        @scenario(sp, r = RAINFALL, x <= x0 + r - y)
    end
    stageobjective!(sp, 5 * p)
end

status = solve(m, max_iterations = 10)

@test getbound(m) == 10
