#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#==
    This example comes from
        https://github.com/blegat/StochasticDualDynamicProgramming.jl/blob/fe5ef82db6befd7c8f11c023a639098ecb85737d/test/prob5.2_2stages.jl
==#

using SDDP, JuMP, Base.Test, Clp

n = 4
m = 3
ic = [16, 5, 32, 2]
C = [25, 80, 6.5, 160]
T = [8760, 7000, 1500] / 8760
D2 = [diff([0, 3919, 7329, 10315])  diff([0, 7086, 9004, 11169])]
p2 = [0.9, 0.1]

mod = SDDPModel(
                  sense = :Min,
                 stages = 2,
                 solver = ClpSolver(),
        objective_bound = 0,
        # noise_probability = [ Float64[], p2, p2 ]
                                ) do sp, t

    @state(sp, x[i=1:n] >= 0, x0 == 0)

    @variables(sp, begin
        y[1:n, 1:m] >= 0
        v[1:n]      >= 0
        penalty     >= 0 # to relax constraint
    end)

    @constraints(sp, begin
        x .== x0 + v
        [i=1:n], sum(y[i, :]) <= x0[i]
    end)

    stageobjective!(sp, dot(ic, v) +  dot(C, y * T) + 1e6 * penalty)

    if t != 1 # no uncertainty in first stage
        for j in 1:m
            @rhsnoise(sp, s=1:size(D2, 2), sum(y[:,j]) + penalty >= D2[j,s])
        end
        setnoiseprobability!(sp, p2)
    end
    if t == 2
        @constraint(sp, sum(v) == 0) # no investment in last stage
    end
end

@time status = SDDP.solve(mod,
    max_iterations = 50,
    print_level = 0,
    simulation = MonteCarloSimulation(
        frequency = 10,
        min       = 100,
        step      = 1,
        max       = 100
    )
)

@test isapprox(getbound(mod), 340315.52, atol=0.1)
sim = simulate(mod, 1, [:x, :penalty])
@test length(sim) == 1
@test isapprox(sim[1][:x][1], [5085,1311,3919,854])
