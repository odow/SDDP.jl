#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using JuMP, SDDP, Ipopt, Base.Test

m = SDDPModel(
    stages          = 3,
    sense           = :Min,
    solver          = IpoptSolver(print_level=0),
    objective_bound = 0.0
                ) do sp, t
    @states(sp, begin
        200 >= volume′ >= 0,  volume==200.0
               inflow′ >= 10, inflow==50.0
    end)
    @variables(sp, begin
        hydro_generation   >= 0
        hydro_spill        >= 0
        thermal_generation >= 0
        inflow_noise_term  >= 0.8
    end)
    @constraints(sp, begin
        volume′ == volume - hydro_generation - hydro_spill + inflow′
        hydro_generation + thermal_generation >= 150.0
    end)

    @rhsnoise(sp, ω = [0.9, 1.0, 1.1], inflow_noise_term == ω)
    @NLconstraint(sp, log(inflow′) == log(inflow) + log(inflow_noise_term))

    fuel_cost = [50.0, 100.0, 150.0]
    @stageobjective(sp, fuel_cost[t] * thermal_generation )
end

srand(123)
solve(m, iteration_limit=10, print_level=0)
@test isapprox(getbound(m), 5585.8, atol=1e-1)
s = simulate(m, 1, [:inflow′])
ω = [0.9, 1.0, 1.1]
@test isapprox(
    log(s[1][:inflow′][3]),
    log(s[1][:inflow′][2]) + log(ω[s[1][:noise][3]]),
    atol=1e-4
)
