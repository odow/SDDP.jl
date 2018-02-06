#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
    Example: newsvendor.

    This example is based on the classical newsvendor problem, but features
    and AR(1) spot-price.

    V(x[t-1], ω[t]) =         max p[t] × u[t]
                       subject to x[t] = x[t-1] - u[t] + ω[t]
                                  u[t] ∈ [0, 1]
                                  x[t] ≥ 0
                                  p[t] = p[t-1] + ϕ[t]

    x[0] = 2.0
    p[0] = 1.5
    ω[t] ~ {0, 0.05, 0.10, ..., 0.45, 0.5} with uniform probability.
    ϕ[t] ~ {-0.25, -0.125, 0.125, 0.25} with uniform probability.
=#

using SDDP, JuMP, Clp, Base.Test

function newsvendor_example(DISCRETIZATION = 1)

    srand(10)

    MIN_PRICE     = 0.75
    INITIAL_PRICE = 1.50
    MAX_PRICE     = 2.25
    NOISES        = DiscreteDistribution([-0.25, -0.125, 0.125, 0.25])

    function pricedynamics(price, noise, stage, markov)
        clamp(price + noise, MIN_PRICE, MAX_PRICE)
    end

    value_function = if DISCRETIZATION == 1
        DynamicPriceInterpolation(
            dynamics       = pricedynamics,
            initial_price  = INITIAL_PRICE,
            min_price      = MIN_PRICE,
            max_price      = MAX_PRICE,
            noise          = NOISES,
        lipschitz_constant = 2.0
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
        @state(sp, x >= 0, x0 == 2)
        @variable(sp, 0 <= u <= 1)
        @rhsnoise(sp, ω = 0:0.05:0.5, x  == x0 - u + ω)
        @stageobjective(sp, price -> price * u)
    end
    return m
end

for discretization in [1, 3]
    m = newsvendor_example(discretization)
    srand(123)
    status = SDDP.solve(m,
        max_iterations = 50,
        cut_output_file="price.cuts",
        print_level=0,
        # test cut selection with default oracle
        cut_selection_frequency = 25
    )
    @test status == :max_iterations
    @test getbound(m) .<= 4.098

    # test round-tripping the cuts
    m2 = newsvendor_example(discretization)
    try
        SDDP.loadcuts!(m2,"price.cuts")
    finally
        rm("price.cuts")
    end
    # we have to solve once more to initialize the second model
    SDDP.solve(m, max_iterations=1, print_level=0)
    SDDP.solve(m2, max_iterations=1, print_level=0)
    @test isapprox(getbound(m2), getbound(m), atol=1e-4)

    # test value function plotting without launching a plot
    sp = SDDP.getsubproblem(m, 2, 1)
    vf = SDDP.valueoracle(sp)
    x, y = SDDP.processvaluefunctiondata(vf, false, 0:0.1:3, 1:0.05:2)
    @test size(x) == (2,651)
    @test size(y) == (651,)
    # SDDP.plotvaluefunction(m, 2, 1, 0:0.1:3, 1:0.05:2)

    x, y = SDDP.processvaluefunctiondata(vf, false, 0:0.1:3, 1.5)
    @test size(x) == (1,31)
    @test size(y) == (31,)
    # SDDP.plotvaluefunction(m, 2, 1, 0:0.1:3, 1.5)

    x, y = SDDP.processvaluefunctiondata(vf, false, 1.5, 1:0.05:2)
    @test size(x) == (1,21)
    @test size(y) == (21,)
    # SDDP.plotvaluefunction(m, 2, 1, 1.5, 1:0.05:2)

    @test_throws Exception SDDP.processvaluefunctiondata(vf, false, 0:0.1:3)

end
