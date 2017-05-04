#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

# For repeatability
srand(11111)

sampled_errors = [-0.1290, -0.1010, -0.0814, -0.0661, -0.0530, -0.0412, -0.0303, -0.0199, -0.00987, 0.0, 0.00987, 0.0199, 0.0303, 0.0412, 0.0530, 0.0661, 0.0814, 0.1010, 0.1290]

σ² = linspace(1, 0, 28) # Decreasing variance in changes in price over time
transaction_cost = 0.01

ribs = collect(linspace(3, 9, 5))

box(x, a, b) = min(b, max(a, x))
price_dynamics(p, w, t, i) = box(1.01 * exp(log(p) + σ²[t]*w), 3, 9)

m = SDDPModel(
    sense             = :Max,
    stages            = 28,
    objective_bound   = 20,
    solver            = ClpSolver(),
    value_function    = InterpolatedValueFunction(
                            # dynamics can't depend on other things
                            dynamics       = price_dynamics,
                            initial_price  = 4.50,
                            rib_locations  = ribs,
                            noise          = Noise(sampled_errors)
                        )
                                            ) do sp, t

    # create state variables
    @states(sp, begin
        0 <= contracts_sold  <= 1.5, contracts_sold0 == 0
        0 <= production <= 1.5, production0 == 0
    end)

    # auxillary variables
    @variables(sp, begin
        0 <= sell <= 1.2
        output >= 0
    end)

    # constraints
    @constraints(sp, begin
        contracts_sold  == contracts_sold0 + sell
        production == production0 + output
        contracts_sold0 + sell <= 1.2
    end)

    # a constraint with varying RHS (but we leverage the JuMP tooling to evaluate that)
    @scenario(sp,
        alpha = linspace(0., 0.05, 5),
        output <= alpha
    )

    if t < 28
        stageobjective!(sp,
            price -> (sell * (price - transaction_cost)) # returns AffExpr for stage objective
        )
    else
        stageobjective!(sp,
            price -> (production - contracts_sold0) * price # objective
        )
    end
end

SDDP.solve(m,
    max_iterations = 15,
    simulation = MonteCarloSimulation(
        frequency = 15,
        min       = 500,
        max       = 500
    )
)

@test getbound(m) >= 4.1
@test getbound(m) <= 4.3
