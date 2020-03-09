#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Test, Distributions

function infinite_lin_HD()
    graph = SDDP.Graph(
        :root_node,
        [:week],
        [(:root_node => :week, 1.0), (:week => :week, 0.95)],
    )
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, node
        @variable(subproblem, -10 <= state <= 10, SDDP.State, initial_value = 0)
        @variables(subproblem, begin
            0 <= order_quantity
            0 <= lost_demand
            0 <= disposed_units
            0 <= backordered_units
            0 <= held_units
            demand
        end)
        @constraint(subproblem, backordered_units >= -state.out)
        @constraint(subproblem, held_units >= state.out)
        @constraint(
            subproblem,
            state.out == state.in + order_quantity - demand + lost_demand - disposed_units
        )
        # Truncated normal on [0, 10] with mean 5 and sd 2.
        Pg = rand(Distributions.TruncatedNormal(5, 2, 0, 10), 50)
        sort!(Pg)
        SDDP.parameterize(subproblem, Pg) do ω
            JuMP.fix(demand, ω)
        end
        @stageobjective(
            subproblem,
            20 * order_quantity +  # Ordering cost cp.
            2 * held_units +  # Holding cost ch.
            10 * backordered_units + # Backorder cost cb.
            10 * disposed_units +  # Disposal cost cd.
            100 * lost_demand  # Lost demand cost cl.
        )
    end
    return model
end

function infinite_lin_DH()
    graph = SDDP.Graph(
        :root_node,
        [:decision, :recourse],
        [
            (:root_node => :decision, 1.0),
            (:decision => :recourse, 1.0),
            (:recourse => :decision, 0.95),
        ],
    )
    model = SDDP.PolicyGraph(
        graph,
        bellman_function = SDDP.BellmanFunction(lower_bound = 0),
        optimizer = GLPK.Optimizer,
    ) do subproblem, node
        @variable(subproblem, -10 <= state <= 10, SDDP.State, initial_value = 0)
        @variable(subproblem, 0 <= order_quantity, SDDP.State, initial_value = 0)
        if node == :decision
            @constraint(subproblem, state.out == state.in)
            @stageobjective(subproblem, 20 * order_quantity.out)
        else
            @variables(subproblem, begin
                0 <= lost_demand
                0 <= disposed_units
                0 <= backordered_units
                0 <= held_units
                demand
            end)
            @constraints(
                subproblem,
                begin
                    state.out ==
                    state.in + order_quantity.in - demand + lost_demand - disposed_units
                    backordered_units >= -state.out
                    held_units >= state.out
                end
            )
            Pg = rand(Distributions.TruncatedNormal(5, 2, 0, 10), 50)
            sort!(Pg)
            SDDP.parameterize(subproblem, Pg) do ω
                JuMP.fix(demand, ω)
            end
            @stageobjective(
                subproblem,
                2 * held_units +  # Holding cost ch.
                10 * backordered_units + # Backorder cost cb.
                10 * disposed_units +  # Disposal cost cd.
                100 * lost_demand  # Lost demand cost cl.
            )
        end
    end
    return model
end

using Random
Random.seed!(1234)
begin
    model = infinite_lin_HD()
    SDDP.train(model, iteration_limit = 75, print_level = 1)
    results = SDDP.simulate(model, 500)
    objectives = [sum(s[:stage_objective] for s in simulation) for simulation in results]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    println("HD Confidence_interval = $(sample_mean) ± $(sample_ci)")
end
Random.seed!(1234)
begin
    model = infinite_lin_DH()
    SDDP.train(model, iteration_limit = 75, print_level = 1)
    results = SDDP.simulate(model, 500)
    objectives = [sum(s[:stage_objective] for s in simulation) for simulation in results]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    println("DH Confidence_interval = $(sample_mean) ± $(sample_ci)")
end
