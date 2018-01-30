#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#==
    An implementation of the QuickStart example from FAST
    https://github.com/leopoldcambier/FAST/tree/daea3d80a5ebb2c52f78670e34db56d53ca2e778/demo

==#

using SDDP, JuMP, Clp, Base.Test

m = SDDPModel(
                sense  = :Min,
                stages = 2,
                solver = ClpSolver(),
                objective_bound = -5
                    ) do sp, t
    @state(sp, x >= 0, x0==0)
    if t == 1
        @stageobjective(sp, x)
    else
        @variable(sp, s >= 0)
        @constraint(sp, s <= x0)
        @rhsnoise(sp, d = [2, 3], s <= d)
        @stageobjective(sp, -2s)
    end
end

status = solve(m, max_iterations = 3,
print_level = 0)

@test getbound(m) == -2
