#  Copyright 2017-20, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Random, Statistics, Test

function inventory_management_problem()
    demand_values = [1.0, 2.0]
    demand_prob = Dict(:Ah => [0.2, 0.8], :Bh => [0.8, 0.2])
    graph = SDDP.Graph(
        :root_node,
        [:Ad, :Ah, :Bd, :Bh],
        [
            (:root_node => :Ad, 0.5),
            (:root_node => :Bd, 0.5),
            (:Ad => :Ah, 1.0),
            (:Ah => :Ad, 0.8),
            (:Ah => :Bd, 0.1),
            (:Bd => :Bh, 1.0),
            (:Bh => :Bd, 0.8),
            (:Bh => :Ad, 0.1),
        ],
    )
    SDDP.add_ambiguity_set(graph, [:Ad, :Bd], 1e2)
    SDDP.add_ambiguity_set(graph, [:Ah, :Bh], 1e2)

    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, node
        @variables(subproblem, begin
                0 <= inventory <= 2, (SDDP.State, initial_value = 0.0)
                buy >= 0
                demand
            end)
        @constraint(subproblem, demand == inventory.in - inventory.out + buy)
        if node == :Ad || node == :Bd || node == :D
            JuMP.fix(demand, 0)
            @stageobjective(subproblem, buy)
        else
            SDDP.parameterize(subproblem, demand_values, demand_prob[node]) do ω
                JuMP.fix(demand, ω)
            end
            @stageobjective(subproblem, 2 * buy + inventory.out)
        end
    end

    # Train the policy.
    Random.seed!(123)
    SDDP.train(model; iteration_limit = 100, print_level = 0, cut_type = SDDP.SINGLE_CUT)
    results = SDDP.simulate(model, 500)
    objectives = [sum(s[:stage_objective] for s in simulation) for simulation in results]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    @test SDDP.calculate_bound(model) ≈ sample_mean atol = sample_ci
    return model
end

inventory_management_problem()
