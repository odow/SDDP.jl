#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

function contracting_example(DISCRETIZATION = 1)
    srand(10)

    MIN_PRICE     = 3.0
    INITIAL_PRICE = 4.50
    MAX_PRICE     = 9.0
    NOISES        = DiscreteDistribution([-0.1290, -0.1010, -0.0814, -0.0661, -0.0530,
        -0.0412, -0.0303, -0.0199, -0.00987, 0.0, 0.00987, 0.0199, 0.0303, 0.0412,
        0.0530, 0.0661, 0.0814, 0.1010, 0.1290])

    function pricedynamics(price, noise, stage, markov)
         # Decreasing variance in changes in price over time
        σ² = linspace(1, 0, 28)
        next_price = 1.01 * exp(log(price) + σ²[stage]*noise)
        min(MAX_PRICE,max(MIN_PRICE, next_price))
    end

    value_function = if DISCRETIZATION == 1
        DynamicPriceInterpolation(
            dynamics       = pricedynamics,
            initial_price  = INITIAL_PRICE,
            min_price      = MIN_PRICE,
            max_price      = MAX_PRICE,
            noise          = NOISES,
            cut_oracle = SDDP.NanniciniOracle(typeof(INITIAL_PRICE), 20)
        )
    else
        StaticPriceInterpolation(
            dynamics       = pricedynamics,
            initial_price  = INITIAL_PRICE,
            rib_locations  =  collect(linspace(MIN_PRICE, MAX_PRICE, DISCRETIZATION)),
            noise          = NOISES
        )
    end

    m = SDDPModel(
        sense             = :Max,
        stages            = 28,
        objective_bound   = 5,
        solver            = ClpSolver(),
        value_function    = value_function
                                            ) do sp, t

        transaction_cost = 0.01

        # create state variables
        @states(sp, begin
            0 <= contracts_sold  <= 1.5, contracts_sold0 == 0
            0 <= production <= 1.5, production0 == 0
        end)

        # auxillary variables
        @variables(sp, begin
            0 <= sell <= 1.2
        end)

        # constraints
        @constraints(sp, begin
            contracts_sold  == contracts_sold0 + sell
            contracts_sold0 + sell <= 1.2
        end)

        @rhsnoise(sp, output = linspace(0., 0.05, 5), production == production0 + output)

        if t < 28
            @stageobjective(sp, price -> (sell * (price - transaction_cost)))
        else
            @stageobjective(sp, price -> (production - contracts_sold0) * price)
        end
    end
    return m
end

# dynamic interpolation
m = contracting_example()
srand(123)
SDDP.solve(m, max_iterations = 50, cut_selection_frequency=10, print_level=2)
@test SDDP.getbound(m) .<= 5.0

# 3 fixed ribs
m = contracting_example(3)
srand(123)
SDDP.solve(m, max_iterations = 10)
@test SDDP.getbound(m) .<= 4.071

# 5 fixed ribs
m = contracting_example(5)
srand(123)
SDDP.solve(m, max_iterations = 10)
@test SDDP.getbound(m) .<= 4.12
