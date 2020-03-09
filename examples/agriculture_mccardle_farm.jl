#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Test

"""
A farm planning problem. There are four stages. The first stage is a
deterministic  planning stage. The next three are wait-and-see operational
stages. The uncertainty in the three operational stages is a Markov chain for
weather. There are three Markov states: dry, normal, and wet.

Inspired by
  R. McCardle, Farm management optimization. Masters thesis, University
  of Louisville, Louisville, Kentucky, United States of America (2009)

All data, including short variable names, is taken from that thesis.
"""
function test_mccardle_farm_model()
    S = [  # cutting, stage
        0 1 2
        0 0 1
        0 0 0
    ]
    t = [60, 60, 245]  # days in period
    D = [210, 210, 858]  # demand
    q = [  # selling price per bale
        [4.5 4.5 4.5; 4.5 4.5 4.5; 4.5 4.5 4.5],
        [5.5 5.5 5.5; 5.5 5.5 5.5; 5.5 5.5 5.5],
        [6.5 6.5 6.5; 6.5 6.5 6.5; 6.5 6.5 6.5],
    ]
    b = [  # predicted yield (bales/acres) from cutting i in weather j.
        30 75 37.5
        15 37.5 18.25
        7.5 18.75 9.325
    ]
    w = 3000  # max storage
    C = [50 50 50; 50 50 50; 50 50 50]  # cost to grow hay
    r = [  # Cost per bale of hay from cutting i during weather condition j.
        [5 5 5; 5 5 5; 5 5 5],
        [6 6 6; 6 6 6; 6 6 6],
        [7 7 7; 7 7 7; 7 7 7],
    ]
    M = 60.0  # max acreage for planting
    H = 0.0  # initial inventory
    V = [0.05, 0.05, 0.05]  # inventory cost
    L = 3000.0  # max demand for hay

    graph = SDDP.MarkovianGraph([
        ones(Float64, 1, 1),
        [0.14 0.69 0.17],
        [0.14 0.69 0.17; 0.14 0.69 0.17; 0.14 0.69 0.17],
        [0.14 0.69 0.17; 0.14 0.69 0.17; 0.14 0.69 0.17],
    ])

    model = SDDP.PolicyGraph(
        graph,
        bellman_function = SDDP.BellmanFunction(lower_bound = 0.0),
        optimizer = GLPK.Optimizer,
    ) do subproblem, index
        stage, weather = index
        # ===================== State Variables =====================
        # Area planted.
        @variable(subproblem, 0 <= acres <= M, SDDP.State, initial_value = M)
        @variable(
            subproblem,
            bales[i = 1:3] >= 0,
            SDDP.State,
            initial_value = (i == 1 ? H : 0)
        )
        # ===================== Variables =====================
        @variables(subproblem, begin
            buy[1:3] >= 0  # Quantity of bales to buy from each cutting.
            sell[1:3] >= 0 # Quantity of bales to sell from each cutting.
            eat[1:3] >= 0  # Quantity of bales to eat from each cutting.
            pen_p[1:3] >= 0  # Penalties
            pen_n[1:3] >= 0  # Penalties
        end)
        # ===================== Constraints =====================
        if stage == 1
            @constraint(subproblem, acres.out <= acres.in)
            @constraint(subproblem, [i = 1:3], bales[i].in == bales[i].out)
        else
            @expression(
                subproblem,
                cut_ex[c = 1:3],
                bales[c].in + buy[c] - eat[c] - sell[c] + pen_p[c] - pen_n[c]
            )
            @constraints(
                subproblem,
                begin
                    # Cannot plant more land than previously cropped.
                    acres.out <= acres.in
                    # In each stage we need to meet demand.
                    sum(eat) >= D[stage-1]
                    # We can buy and sell other cuttings.
                    bales[stage-1].out == cut_ex[stage-1] + acres.in * b[stage-1, weather]
                    [c = 1:3; c != stage - 1], bales[c].out == cut_ex[c]
                    # There is some maximum storage.
                    sum(bales[i].out for i = 1:3) <= w
                    # We can only sell what is in storage.
                    [c = 1:3], sell[c] <= bales[c].in
                    # Maximum sales quantity.
                    sum(sell) <= L
                end
            )
        end
        # ===================== Stage objective =====================
        if stage == 1
            @stageobjective(subproblem, 0.0)
        else
            @stageobjective(
                subproblem,
                1000 * (sum(pen_p) + sum(pen_n)) +
                # cost of growing
                C[stage-1, weather] * acres.in +
                sum(
                    # inventory cost
                    V[stage-1] * bales[cutting].in * t[stage-1] +
                    # purchase cost
                    r[cutting][stage-1, weather] * buy[cutting] +
                    # feed cost
                    S[cutting, stage-1] * eat[cutting] -
                    # sell reward
                    q[cutting][stage-1, weather] * sell[cutting] for cutting = 1:3
                )
            )
        end
        return
    end
    SDDP.train(model, iteration_limit = 50, print_level = 0)
    @test SDDP.termination_status(model) == :iteration_limit
    @test SDDP.calculate_bound(model) ≈ 4074.1391 atol = 1e-5
end

test_mccardle_farm_model()
