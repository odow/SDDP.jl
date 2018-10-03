#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using SDDP, JuMP, Clp, Base.Test

function solve_model(noise_probability)
    m = SDDPModel(
        sense = :Max,
        stages = 2,
        objective_bound = 5.0,
        cut_oracle = LevelOneCutOracle(),
        solver = ClpSolver()
                                ) do sp, t

        @state(sp, x >= 0, x0==1.5)
        @variable(sp, 0 <= u <= 1)
        @constraint(sp, x == x0 - u)
        if t == 1
            @stageobjective(sp, 2 * u)
        else
            @stageobjective(sp, p=[1, 3], p * u)
            setnoiseprobability!(sp, noise_probability)
        end
    end
    solve(m,  iteration_limit = 5, print_level=0, cut_selection_frequency = 5)
    m
end

srand(111)
m = solve_model([0.5, 0.5])
@test isapprox(getbound(m), 3.0, atol=1e-3)

SDDP.writecuts!("cuts_valid.csv", m; onlyvalid=true)
SDDP.writecuts!("cuts_all.csv",   m; onlyvalid=false)
valid_cuts = readdlm("cuts_valid.csv", ',')
all_cuts   = readdlm("cuts_all.csv",   ',')
@test size(valid_cuts, 2) == size(all_cuts, 2)
@test size(valid_cuts, 1) < size(all_cuts, 1)
rm("cuts_valid.csv")
rm("cuts_all.csv")

m = solve_model([1.0, 0.0])
@test isapprox(getbound(m), 2.5, atol=1e-3)

m = solve_model([0.0, 1.0])
@test isapprox(getbound(m), 4.0, atol=1e-3)

m = solve_model([0.25, 0.75])
@test isapprox(getbound(m), 3.5, atol=1e-3)

@test_throws Exception solve_model([0.75, 0.75])
