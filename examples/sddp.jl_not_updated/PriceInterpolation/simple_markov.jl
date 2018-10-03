#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

function buildvaluefunction(stage::Int, markov::Int)
    if stage == 1
        NOISES = DiscreteDistribution([-0.25, -0.125, 0.125, 0.25])
    else
        if markov == 1
            NOISES = DiscreteDistribution([-0.25, -0.125])
        else
            NOISES = DiscreteDistribution([0.125, 0.25])
        end
    end
    return DynamicPriceInterpolation(
        dynamics       = (price, noise) -> price + noise,
        initial_price  = 1.5,
        min_price      = 0.75,
        max_price      = 2.25,
        noise          = NOISES,
    lipschitz_constant = 2.0
    )
end

m = SDDPModel(
    sense             = :Max,
    stages            = 3,
    objective_bound   = 10,
    solver            = ClpSolver(),
    # The transition matrix
    #    x - x
    #   / \ /
    # x    X
    #   \ / \
    #    x - x
    markov_transition = Array{Float64, 2}[
        [1.0]',
        [0.5 0.5],
        [0.5 0.5; 0.5 0.5]
    ],
    value_function    = buildvaluefunction
                                        ) do sp, t, i
    @state(sp, x >= 0, x0 == 2)
    @variable(sp, 0 <= u <= 1)
    @rhsnoise(sp, ω = 0:0.05:0.5, x  == x0 - u + ω)
    @stageobjective(sp, price -> price * u)
end

solve(m, iteration_limit = 50, print_level = 0)
s = simulate(m, 100, [:price, :x, :u])
for si in s
    if si[:markov][2] == 1
        @test si[:price][2] < si[:price][1]
    else
        @test si[:price][2] > si[:price][1]
    end
end
