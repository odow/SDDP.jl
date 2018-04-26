#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    Example 1.3.2 from
        Bertsekas, D. (2005). Dynamic Programming and Optimal Control:
        Volume I (3rd ed.). Bellmont, MA: Athena Scientific.
=#

using SDDP, JuMP, Ipopt, Base.Test

m = SDDPModel(
        stages          = 3,
        sense           = :Min,
        solver          = IpoptSolver(print_level=0),
        objective_bound = 0.0
            ) do sp, t
    @state(sp, 0 <= x <= 2, x0 == 0)
    @variables(sp, begin
        0 <= u <= 2
             w
    end)
    @rhsnoise(sp, ω=[0.0, 1.0, 2.0], w == ω)
    setnoiseprobability!(sp, [0.1, 0.7, 0.2])
    @constraints(sp, begin
        x0 + u     <= 2
        x0 + u - w <= x
    end)
    @stageobjective(sp, u + (x0 + u - w)^2)
end

solve(m, iteration_limit=20, print_level=0)
@test isapprox(getbound(m), 1.375)
