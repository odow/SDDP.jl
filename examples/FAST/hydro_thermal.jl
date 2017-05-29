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

    C = 5
    V = 8
    d = 6

    @state(sp, 0 <= x <= V, x0==0)
    @variables(sp, begin
        y >= 0
        p >= 0
        rainfall
    end)
    @constraints(sp, begin
        p + y >= d
        x <= x0 + rainfall - y
    end)

    if t == 1
        @constraint(sp, rainfall ==  mean(RAINFALL))
    else
        @scenario(sp, r = RAINFALL, rainfall == r)
    end

    stageobjective!(sp, C * p)

end

status = solve(m, max_iterations = 10)

@test getbound(m) == 10
