#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using JuMP, SDDP, Ipopt

m = SDDPModel(
    stages          = 2,
    sense           = :Max,
    objective_bound = 2.0,
    solver          = IpoptSolver(print_level=0)
        ) do sp, t

    @state(sp, x′ >= 1, x==1)
    if t == 1
        @rhsnoise(sp, w=[1,2,3], x′ == w)
        @stageobjective(sp, 0.0)
    else
        @variable(sp, y)
        @NLconstraint(sp, y <= log(x))
        @stageobjective(sp, y)
    end
end
solve(m, iteration_limit=10, print_level=0)
@test isapprox(getbound(m), log(6) / 3, atol=1e-4)
