#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#
#    Inspired by
#       R. McCardle, Farm management optimization. Masters thesis, University
#       of Louisville, Louisville, Kentucky, United States of America (2009)
#

using SDDP, JuMP, Clp
using Base.Test

function mccardle_farm_model()
    # cutting, stage
    S = [
        0 1 2;
        0 0 1;
        0 0 0
    ]

    # days in period
    t = [60,60,245]
    # demand
    D = [210,210,858]

    # selling price per bale
    q = [
        [4.5 4.5 4.5; 4.5 4.5 4.5; 4.5 4.5 4.5],
        [5.5 5.5 5.5; 5.5 5.5 5.5; 5.5 5.5 5.5],
        [6.5 6.5 6.5; 6.5 6.5 6.5; 6.5 6.5 6.5]
    ]

    # predicted yield (bales/acres) from cutting i in weather j
    b = [
        30 75 37.5;
        15 37.5 18.25;
        7.5 18.75 9.325
    ]
    # max storage
    w = 3000
    # cost to grow hay
    C = [50 50 50; 50 50 50; 50 50 50]
    # Cost per bale of hay from cutting i during weather condition j;
    r = [
        [5 5 5; 5 5 5; 5 5 5],
        [6 6 6; 6 6 6; 6 6 6],
        [7 7 7; 7 7 7; 7 7 7]
    ]

    # max acreage for planting
    M = 60.0
    # initial inventory
    H = 0.0
    # inventory cost
    V = [0.05, 0.05, 0.05]
    # max demand for hay
    L = 3000.0

    transition = Array{Float64, 2}[
        [1.0]',
        [0.14 0.69 0.17; 0.14 0.69 0.17; 0.14 0.69 0.17],
        [0.14 0.69 0.17; 0.14 0.69 0.17; 0.14 0.69 0.17],
        [0.14 0.69 0.17; 0.14 0.69 0.17; 0.14 0.69 0.17]
    ]

    m = SDDPModel(
                 stages = 4,
        objective_bound = [0.0, 0.0, 0.0, 0.0],
                  sense = :Min,
      markov_transition = transition,
                 solver = ClpSolver()
                            ) do sp, stage, weather
        @states(sp, begin
            # acres planted for each cutting
            0 <= acres <= M, acres0==M
            # bales from cutting i in storage
            bales[cutting=1:3] >= 0, bales0==[H,0,0][cutting]
        end)

        @variables(sp, begin
            # quantity of bales to buy from cutting
            buy[cutting=1:3] >= 0
            # quantity of bales to sell from cutting
            sell[cutting=1:3] >= 0
            # quantity of bales to eat from cutting
            eat[cutting=1:3] >= 0
            # penalties
            pen_p[cutting=1:3] >= 0
            pen_n[cutting=1:3] >= 0
        end)
        @expression(sp, total_penalties, sum(pen_p) + sum(pen_n))

        if stage == 1
            @constraint(sp, bales0 .== bales)
            @stageobjective(sp, 0.0)
        else
            @expression(sp, cut_ex[c=1:3], bales0[c] + buy[c]  - eat[c] - sell[c] + pen_p[c] - pen_n[c])
            @constraints(sp, begin
                # plan planting for next stage
                acres <= acres0
                # meet demand
                sum(eat) >= D[stage-1]

                bales[stage-1] == cut_ex[stage-1] + acres0 * b[stage-1, weather]
                # can buy and sell other cuttings
                [c=1:3;c!=stage-1], bales[c] == cut_ex[c]

                # max storage
                sum(bales) <= w

                # can only sell what is in storage
                [c=1:3], sell[c] <= bales0[c]

                sum(sell) <= L
            end)
            @stageobjective(sp,
                1000 * total_penalties +
                # cost of growing
                C[stage-1, weather] * acres0 +
                sum(
                    # inventory cost
                    V[stage-1] * bales0[cutting] * t[stage-1] +
                    # purchase cost
                    r[cutting][stage-1,weather] * buy[cutting] +
                    # feed cost
                    S[cutting,stage-1] * eat[cutting] -
                    # sell reward
                    q[cutting][stage-1,weather] * sell[cutting]
                for cutting in 1:3)
            )
        end
    end
end

srand(111)
m = mccardle_farm_model()
solution = solve(m, max_iterations=20)
@test isapprox(getbound(m), 4074.1391, atol=1e-5)

# results = simulate(m,  # Simulate the policy
#     100,               # number of monte carlo realisations
#     [:acres,:acres0,:buy,:sell,:eat,:pen_p,:pen_n,:bales,:bales0]       # variables to return
#     )
# plt = SDDP.newplot()
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->results[i][:stageobjective][t], title="Accumulated Profit", ylabel="Accumulated Profit (\$)", cumulative=true)
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->results[i][:stageobjective][t], title="Weekly Income",      ylabel="Week Profit (\$)")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->results[i][:acres][t],          title="acres")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->results[i][:acres0][t],         title="acres0")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:bales][t]),     title="bales")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:bales0][t]),    title="bales0")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:buy][t]),       title="buy")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:sell][t]),      title="sell")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:eat][t]),       title="eat")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:pen_p][t]),     title="pen_p")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->sum(results[i][:pen_n][t]),     title="pen_n")
# SDDP.addplot!(plt, 1:100, 1:4, (i,t)->results[i][:markov][t],         title="weather")
# SDDP.show(plt)
