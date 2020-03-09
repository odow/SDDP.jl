#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using JuMP, SDDP, LinearAlgebra, GLPK, Test

function generation_expansion(integrality_handler)
    build_cost = 1e4
    use_cost = 4
    num_units = 5
    capacities = ones(num_units)
    demand_vals =
        0.5 * [
            5 5 5 5 5 5 5 5
            4 3 1 3 0 9 8 17
            0 9 4 2 19 19 13 7
            25 11 4 14 4 6 15 12
            6 7 5 3 8 4 17 13
        ]
    # Cost of unmet demand
    penalty = 5e5
    # Discounting rate
    rho = 0.99
    model = SDDP.LinearPolicyGraph(
        stages = 5,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
        integrality_handler = integrality_handler,
    ) do sp, stage

        @variable(sp, 0 <= invested[1:num_units] <= 1, SDDP.State, Int, initial_value = 0)
        @variables(sp, begin
            generation >= 0
            unmet >= 0
            demand
        end)

        @constraints(
            sp,
            begin
                # Can't un-invest
                investment[i in 1:num_units], invested[i].out >= invested[i].in
                # Generation capacity
                sum(capacities[i] * invested[i].out for i = 1:num_units) >= generation
                # Meet demand or pay a penalty
                unmet >= demand - sum(generation)
                # For fewer iterations order the units to break symmetry, units are identical (tougher numerically)
                # [j in 1:(num_units - 1)], invested[j].out <= invested[j + 1].out
            end
        )

        # Demand is uncertain
        SDDP.parameterize(ω -> JuMP.fix(demand, ω), sp, demand_vals[stage, :])

        @expression(
            sp,
            investment_cost,
            build_cost * sum(invested[i].out - invested[i].in for i = 1:num_units)
        )
        @stageobjective(
            sp,
            (investment_cost + generation * use_cost) * rho^(stage - 1) + penalty * unmet
        )
    end
    SDDP.train(model, iteration_limit = 50, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 2.078860e6 atol = 1e3
end

for integrality_handler in [
# Solve a continuous relaxation only, tough for SDDiP
# SDDP.SDDiP(),
SDDP.ContinuousRelaxation()]
    generation_expansion(integrality_handler)
end
