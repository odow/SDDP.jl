#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

# For repeatability
srand(11111)

m = SDDPModel(
    sense             = :Max,
    stages            = 2,
    objective_bound   = 5,
    # markov_transition = markov_transition,
    solver            = ClpSolver(),
    value_function    = InterpolatedValueFunction(
                            # dynamics can't depend on other things
                            dynamics       = (p, w, t, i) -> p + w,
                            initial_price  = 1.50,
                            rib_locations  = [1.0, 2.0],
                            noise          = Noise([-0.25, 0.25])
                        )
                                            ) do sp, t

    # create state variables
    @states(sp, begin
        0 <= x  <= 1.5, x0 == 1
    end)

    # auxillary variables
    @variables(sp, begin
        u >= 0
    end)

    # constraints
    @constraints(sp, begin
        x  == x0 - u
    end)

    stageobjective!(sp,
        price -> price * u
    )
end

SDDP.solve(m,
    max_iterations = 5,
    simulation_frequency = 2,
    simulation_max = 20
)

@test isapprox(m.log[end].bound, 1.5)
