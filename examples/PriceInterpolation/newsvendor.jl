#  Copyright 2017, Oscar Dowson

using SDDP, JuMP, Clp, Base.Test

function newsvendor_example(DISCRETIZATION = 1)

    srand(10)

    MIN_PRICE     = 0.75
    INITIAL_PRICE = 1.50
    MAX_PRICE     = 2.25
    NOISES        = DiscreteDistribution([-0.25, -0.125, 0.125, 0.25])

    function pricedynamics(price, noise, stage, markov)
        min(MAX_PRICE,max(MIN_PRICE, price + noise))
    end

    value_function = if DISCRETIZATION == 1
        DynamicPriceInterpolation(
            dynamics       = pricedynamics,
            initial_price  = INITIAL_PRICE,
            min_price      = MIN_PRICE,
            max_price      = MAX_PRICE,
            noise          = NOISES
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
        stages            = 3,
        objective_bound   = 10,
        solver            = ClpSolver(),
        value_function    = value_function
                                            ) do sp, t

        # create state variables
        @states(sp, begin
            0 <= x  <= 3, x0 == 2
        end)

        # auxillary variables
        @variables(sp, begin
            0 <= u <= 1.0
        end)

        # constraints
        @rhsnoise(sp, w = linspace(0, 0.5, 20), x  == x0 - u + w)

        # stageobjective!(sp, 1.5 * u)
        @stageobjective(sp,
            price -> price * u
        )
    end
    return m
end

m = newsvendor_example()
status = SDDP.solve(m, max_iterations = 50, time_limit     = 20.0)
@test status == :max_iterations
@test getbound(m) .<= 4.098
# SDDP.plotvaluefunction(m, 2, 1, linspace(0.0, 1, 50), linspace(1, 2, 50), label1 = "x", label2="price")
# SDDP.plotvaluefunction(m, 2, 1, 0.0:0.25:3, label1 = "x")
