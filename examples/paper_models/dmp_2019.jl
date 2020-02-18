#  Copyright 2017-20, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This file implements the example described in
#   Dowson, O., Morton, D.P., & Pagnoncelli, B. (2019). Partially observable
#   multistage stochastic programming.

using SDDP, GLPK, Random, Statistics, Test

const demand_values = [1.0, 2.0]

const demand_probs = Dict(:Ah => [0.2, 0.8], :Bh => [0.8, 0.2], :H => [0.5, 0.5])

"Create the policy graph for the problem."
function build_graph(model_name)
    if model_name == "hidden" || model_name == "visible"
        graph = SDDP.Graph(
            :root_node,
            [:Ad, :Ah, :Bd, :Bh],
            [
                (:root_node => :Ad, 0.5),
                (:root_node => :Bd, 0.5),
                (:Ad => :Ah, 1.0),
                (:Ah => :Ad, 0.9),
                (:Bd => :Bh, 1.0),
                (:Bh => :Bd, 0.9),
            ],
        )
        if model_name == "hidden"
            SDDP.add_ambiguity_set(graph, [:Ad, :Bd], 1e2)
            SDDP.add_ambiguity_set(graph, [:Ah, :Bh], 1e2)
        end
        return graph
    elseif model_name == "expected_value"
        graph = SDDP.Graph(
            :root_node,
            [:D, :H],
            [(:root_node => :D, 1.0), (:D => :H, 1.0), (:H => :D, 0.9)],
        )
        return graph
    else
        error("Invalid option: model_name=$(model_name)")
    end
end

function solve_inventory_management_problem(model_name, risk_measure)
    graph = build_graph(model_name)
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
            SDDP.parameterize(subproblem, demand_values, demand_probs[node]) do ω
                JuMP.fix(demand, ω)
            end
            @stageobjective(subproblem, 2 * buy + inventory.out)
        end
    end
    Random.seed!(123)
    SDDP.train(model; risk_measure = risk_measure, iteration_limit = 200)
    simulations =
        simulate_policy(model, model_name; terminate_on_leaf = false, discount = true)
    return (model = model, simulations = simulations)
end

function simulate_policy(model, model_name; terminate_on_leaf::Bool, discount::Bool)
    # Simulate policy using same realizations for each policy.
    Ad, Ah, Bd, Bh = :Ad, :Ah, :Bd, :Bh
    if model_name == "expected_value"
        Ad, Ah, Bd, Bh = :D, :H, :D, :H
    end
    Random.seed!(1234)
    scenarios = Any[]
    for replication = 1:2000
        actual_node = rand() < 0.5 ? :A : :B
        num_stages = if terminate_on_leaf
            ceil(Int, log(rand()) / log(0.9))
        else
            50
        end
        scenario = Tuple{Symbol,Union{Nothing,Float64}}[]
        for stage = 1:num_stages
            r = rand()
            if actual_node == :A
                push!(scenario, (Ad, nothing))
                push!(scenario, (Ah, r < demand_probs[:Ah][1] ? 1 : 2))
            else
                push!(scenario, (Bd, nothing))
                push!(scenario, (Bh, r < demand_probs[:Bh][1] ? 1 : 2))
            end
        end
        push!(scenarios, scenario)
    end
    simulations = Any[]
    for scenario in scenarios
        push!(
            simulations,
            SDDP.simulate(
                model,
                1,
                [:inventory, :buy],
                sampling_scheme = SDDP.Historical(scenario),
            )[1],
        )
    end
    function calculate_objective(simulation)
        if discount
            y = 0.0
            ρ = 1.0
            for (i, s) in enumerate(simulation)
                y += ρ * s[:stage_objective]
                if !isodd(i)
                    ρ *= 0.9
                end
            end
            return y
        else
            return sum(s[:stage_objective] for s in simulation)
        end
    end
    objectives = calculate_objective.(simulations)
    return (simulations = simulations, objectives = objectives)
end

function quantile_data(data...)
    return hcat([
        Statistics.quantile(data_i, [0.0, 0.25, 0.5, 0.75, 1.0]) for data_i in data
    ]...)
end

function run_paper_analysis()
    visible = solve_inventory_management_problem("visible", SDDP.Expectation())

    hidden = solve_inventory_management_problem("hidden", SDDP.Expectation())

    expected_value =
        solve_inventory_management_problem("expected_value", SDDP.Expectation())

    risk_averse_expected_value =
        solve_inventory_management_problem("expected_value", SDDP.ModifiedChiSquared(0.25))

    quantiles = quantile_data(
        visible.simulations.objectives,
        hidden.simulations.objectives,
        expected_value.simulations.objectives,
        risk_averse_expected_value.simulations.objectives,
    )
    open("quantiles.dat", "w") do io
        for i = 1:size(quantiles, 1)
            println(io, join(quantiles[i, :], " "))
        end
    end
end

if length(ARGS) > 0
    @assert ARGS[1] == "--run"
    run_paper_analysis()
end
