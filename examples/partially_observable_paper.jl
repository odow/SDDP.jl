#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Gurobi, Random, Statistics, Test

const demand_values = [1.0, 2.0]

const demand_probs = Dict(
    :Ah => [0.2, 0.8],
    :Bh => [0.8, 0.2],
    :H => [0.5, 0.5]
)

"Create the policy graph for the problem."
function build_graph(model_name)
    if model_name == "hidden" || model_name == "visible"
        graph = SDDP.Graph(
            :root_node,
            [:Ad, :Ah, :Bd, :Bh],
            [
                (:root_node => :Ad, 0.5), (:root_node => :Bd, 0.5),
                (:Ad => :Ah, 1.0), (:Ah => :Ad, 0.9),
                (:Bd => :Bh, 1.0), (:Bh => :Bd, 0.9)
            ]
        )
        if model_name == "hidden"
            SDDP.add_partition(graph, [:Ad, :Bd])
            SDDP.add_partition(graph, [:Ah, :Bh])
        end
        return graph
    elseif model_name == "expected_value"
        graph = SDDP.Graph(
            :root_node,
            [:D, :H],
            [(:root_node => :D, 1.0), (:D => :H, 1.0), (:H => :D, 0.9)]
        )
        return graph
    else
        error("Invalid option: model_name=$(model_name)")
    end
end

function solve_inventory_management_problem(model_name, risk_measure)
    graph = build_graph(model_name)
    model = SDDP.PolicyGraph(graph,
                lower_bound = 0.0,
                optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag=0),
                lipschitz_belief = Dict(:Ah => 1e2, :Bh => 1e2, :Ad => 1e2, :Bd => 1e2)
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
    SDDP.train(model; risk_measure=risk_measure, iteration_limit=200)
    simulations = simulate_policy(model, model_name; terminate_on_leaf = false, discount=false)
    expected_value = simulate_policy(model, model_name; terminate_on_leaf = false, discount=true)
    return (model=model, simulations=simulations, expected_value=expected_value)
end

function simulate_policy(model, model_name; terminate_on_leaf::Bool, discount::Bool)
    # Simulate policy using same realizations for each policy.
    Ad, Ah, Bd, Bh = :Ad, :Ah, :Bd, :Bh
    if model_name == "expected_value"
        Ad, Ah, Bd, Bh = :D, :H, :D, :H
    end
    Random.seed!(1234)
    scenarios = Any[]
    for replication in 1:2000
        actual_node = rand() < 0.5 ? :A : :B
        num_stages = if terminate_on_leaf
            ceil(Int, log(rand()) / log(0.9))
        else
            50
        end
        scenario = Tuple{Symbol, Union{Nothing, Float64}}[]
        for stage in 1:num_stages
            r = rand()
            if actual_node == :A
                push!(scenario, (Ad, nothing))
                push!(scenario, (Ah, r < demand_probs[:Ah][1] ? 1 : 2))
            else
                push!(scenario, (Bd, nothing))
                push!(scenario, (Ah, r < demand_probs[:Bh][1] ? 1 : 2))
            end
        end
        push!(scenarios, scenario)
    end
    simulations = Any[]
    for scenario in scenarios
        push!(simulations, SDDP.simulate(
            model, 1, [:inventory, :buy],
            sampling_scheme = SDDP.Historical(scenario))[1]
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
    return (simulations=simulations, objectives=objectives)
end

function get_hidden_value_function(hidden)
    model = hidden[:Ad].subproblem
    belief = hidden[:Ad].belief_state.belief
    JuMP.set_upper_bound(model[:buy], 0)
    B = 0:0.05:1
    X = 0:0.1:2
    if JuMP.has_lower_bound(model[:inventory].out)
        JuMP.delete_lower_bound(model[:inventory].out)
        JuMP.delete_upper_bound(model[:inventory].out)
    end
    Q = zeros(Float64, length(X), length(B))
    for (j, b) in enumerate(B)
        for (i, x) in enumerate(X)
            JuMP.fix(model[:inventory].in, x)
            JuMP.fix(model[:inventory].out, x)
            belief[:Ad] = b
            belief[:Bd] = 1 - b
            SDDP.set_objective(hidden, hidden[:Ad])
            JuMP.optimize!(model)
            Q[i, j] = JuMP.objective_value(model)
        end
    end
    open("value_function.dat", "w") do io
        for (i, x) in enumerate(X)
            for (j, b) in enumerate(B)
                println(io, "$(x) $(b) $(Q[i, j])")
            end
            println(io)
        end
    end
    return X, B, Q
end

function quantile_data(data...)
    return hcat([
        Statistics.quantile(
            data_i,
            [0.0, 0.25, 0.5, 0.75, 1.0]
        ) for data_i in data]...)
end

function run_paper_analysis()
    visible = solve_inventory_management_problem("visible", SDDP.Expectation())

    hidden = solve_inventory_management_problem("hidden", SDDP.Expectation())
    get_hidden_value_function(hidden.model)

    expected_value = solve_inventory_management_problem(
        "expected_value", SDDP.Expectation())

    risk_averse_expected_value = solve_inventory_management_problem(
        "expected_value", SDDP.ModifiedChiSquared(0.25))

    quantiles = quantile_data(
        visible.simulations.objectives,
        hidden.simulations.objectives,
        expected_value.simulations.objectives,
        risk_averse_expected_value.simulations.objectives
    )
    open("quantiles.dat", "w") do io
        for i in 1:size(quantiles, 1)
            println(io, join(quantiles[i, :], " "))
        end
    end

    quantiles = quantile_data(
        visible.expected_value.objectives,
        hidden.expected_value.objectives,
        expected_value.expected_value.objectives,
        risk_averse_expected_value.expected_value.objectives
    )
    open("expected_value_quantiles.dat", "w") do io
        for i in 1:size(quantiles, 1)
            println(io, join(quantiles[i, :], " "))
        end
    end

    open("expected_value_mean.dat", "w") do io
        println(io, Statistics.mean(visible.expected_value.objectives))
        println(io, Statistics.mean(hidden.expected_value.objectives))
        println(io, Statistics.mean(expected_value.expected_value.objectives))
        println(io, Statistics.mean(risk_averse_expected_value.expected_value.objectives))
    end
end

if length(ARGS) > 0
    @assert ARGS[1] == "--run"
    run_paper_analysis()
end

# simulation = expected_value.simulations.simulations[1]
# objectives = map(simulation -> begin
#     y = 0.0
#     for (i, d) in enumerate(simulation)
#         if isodd(i)
#             y += d[:buy]
#         else
#             y += 2 * d[:buy]
#             y += d[:inventory].out
#         end
#     end
#     return y
# end, expected_value.simulations.simulations)
#
# demand = map(simulation -> begin
#     y = 0.0
#     for (i, d) in enumerate(simulation)
#         if !isodd(i)
#             y += d[:noise_term]
#         end
#     end
#     return y
# end, expected_value.simulations.simulations)

objectives = map(simulation -> begin
    y = 0.0
    ρ = 1.0
    for (i, d) in enumerate(simulation)
        if isodd(i)
            y += ρ * d[:buy]
        else
            y += ρ * 2 * d[:buy]
            y += ρ * d[:inventory].out
            ρ *= 1.0  # 0.9
        end
    end
    return y
end, risk_averse_expected_value.expected_value.simulations)

X = map(
    simulation -> map(d -> d[:inventory].out, simulation),
    risk_averse_expected_value.expected_value.simulations
)
function foo(X)
    return extrema(map(simulation -> begin
        y = 0.0
        for (i, d) in enumerate(simulation)
            if !isodd(i)
                y += d[:noise_term]
            end
        end
        return y
    end, X))
end

#
plt =SDDP.SpaghettiPlot(risk_averse_expected_value.simulations.simulations[1:2])
SDDP.add_spaghetti(plt) do data
    data[:inventory].out
end
SDDP.add_spaghetti(plt) do data
    data[:noise_term] === nothing ? 0.0 : data[:noise_term]
end
SDDP.save(plt)
