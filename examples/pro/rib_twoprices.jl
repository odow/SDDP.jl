#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

# For repeatability
srand(11111)

function pricedynamics(p::NTuple{2, Float64}, w::NTuple{2, Float64}, t::Int, i::Int)
    (p[1] + w[1], p[2] + w[2])
end

rib_locations = [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0), (2.0, 2.0)]

p0 = (1.50, 1.50)

noises = [(-0.25, -0.25), (-0.25, 0.25), (0.25, -0.25), (0.25, 0.25)]

m = SDDPModel(
    sense             = :Max,
    stages            = 2,
    objective_bound   = 5,
    # markov_transition = markov_transition,
    solver            = ClpSolver(),
    value_function    = InterpolatedValueFunction(
                            # dynamics can't depend on other things
                            dynamics       = pricedynamics,
                            initial_price  = p0,
                            rib_locations  = rib_locations,
                            noise          = Noise(noises)
                        )
                                            ) do sp, t

    # create state variables
    @states(sp, begin
        0 <= x  <= 1.5, x0 == 1
        0 <= y  <= 1.5, y0 == 1
    end)

    # auxillary variables
    @variables(sp, begin
        u >= 0
        v >= 0
    end)

    # constraints
    @constraints(sp, begin
        x  == x0 - u
        y  == y0 - v
    end)

    stageobjective!(sp,
        price -> price[1] * u + price[2] * v
    )
end

SDDP.solve(m,
    max_iterations = 10,
    simulation = MonteCarloSimulation(
        frequency = 2,
        max       = 20
    )
)

@test isapprox(m.log[end].bound, 3.0)
