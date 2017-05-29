#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

#==
    This example comes from
        https://github.com/blegat/StochasticDualDynamicProgramming.jl/blob/b80eb4d5dfa6f7f4570f225d54adebdd6156c52f/test/prob5.2standalone.jl

    Note: we get a different objective to what is stated, but we show that it is
    equal to solving the deterministic equivalent so we must be correct...
==#

using SDDP, JuMP, Base.Test, Clp

n = 4
m = 3
# Investment cost
ic = [16, 5, 32, 2]
# Fuel cost
C = [25, 80, 6.5, 160]
# Duration
T = [8760, 7000, 1500] / 8760
# Height
D2 = [diff([0, 3919, 7329, 10315])  diff([0, 7086, 9004, 11169])]


mod = Model(solver=ClpSolver())
@variables(mod, begin
    x[1:n] >= 0
    y2[1:n,1:m, 1:2]  >= 0
    y3[1:n,1:m, 1:2, 1:2] >= 0
    penalty>=0
end)
@constraints(mod, begin
    [i=1:n, s=1:2], sum(y2[i,:, s]) <= x[i]
    [i=1:n, s1=1:2, s2=1:2], sum(y3[i,:, s1,s2]) <= x[i]
    [j=1:m, s=1:2], sum(y2[:,j, s])  + penalty >= D2[j,s]
    [j=1:m, s1=1:2,s2=1:2], sum(y3[:,j, s1, s2])  + penalty >= D2[j,s2]
end)
obj = dot(ic, x)
for s1 in 1:2
    append!(obj, 0.5 * dot(C, y2[:,:,s1] * T))
    for s2 in 1:2
        append!(obj, 0.25 * dot(C, y3[:,:,s1,s2] * T))
    end
end
@objective(mod, Min, obj + 1e3 * penalty)
solve(mod)
det_obj = getobjectivevalue(mod)
det_x   = getvalue(x)

mod = SDDPModel(
                  sense = :Min,
                 stages = 3,
                 solver = ClpSolver(),
        objective_bound = 0
                                ) do sp, t

    @state(sp, x[i=1:n] >= 0, x0 == 0)

    @variables(sp, begin
        y[1:n, 1:m] >= 0
        v[1:n] >= 0
        penalty >= 0 # to relax constraint
    end)

    @constraints(sp, begin
        x .== x0 + v
        [i=1:n], sum(y[i, :]) <= x0[i]
    end)

    stageobjective!(sp, dot(ic, v) +  dot(C, y * T) + 1e4 * penalty)

    if t != 1 # no uncertainty in first stage
        for j in 1:m
            @scenario(sp, s=1:size(D2, 2), sum(y[:,j]) + penalty >= D2[j,s])
        end
    end
    if t == 3
        # no investment in last stage
        @constraint(sp, sum(v) == 0)
    end
end

@time status = SDDP.solve(mod,
    max_iterations = 50,
    simulation = MonteCarloSimulation(
        frequency = 10,
        min       = 100,
        step      = 1,
        max       = 100
    )
)


@test isapprox(getbound(mod), det_obj, atol=0.1)

sim = simulate(mod, 1, [:x, :penalty])
@test isapprox(sim[1][:x][1], det_x)
