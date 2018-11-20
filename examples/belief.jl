#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Kokako, GLPK, Random, Statistics, Test

function inventory_management_problem(graph, demand_values, demand_prob,
                                      optimizer=with_optimizer(GLPK.Optimizer))
    model = Kokako.PolicyGraph(graph,
                lower_bound = 0.0,
                optimizer = optimizer,
                lipschitz_belief = Dict(:Ah => 1e2, :Bh => 1e2, :Ad => 1e2, :Bd => 1e2)
                    ) do subproblem, node
        @variables(subproblem, begin
            0 <= inventory <= 2, (Kokako.State, initial_value = 0.0)
            buy >= 0
            demand
        end)
        @constraint(subproblem, demand == inventory.in - inventory.out + buy)
        if node == :Ad || node == :Bd || node == :D
            JuMP.fix(demand, 0)
            @stageobjective(subproblem, buy)
        else
            Kokako.parameterize(subproblem, demand_values, demand_prob[node]) do ω
                JuMP.fix(demand, ω)
            end
            @stageobjective(subproblem, 2 * buy + inventory.out)
        end
    end
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
            Kokako.set_objective(hidden, hidden[:Ad])
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


# ==============================================================================
# Run the script `belief.jl`.
#
# Arguments
#     --num_replications::Int - the number of simulations to conduct.
#     --iteration_limit::Int - the number of iterations to train the policy for.
#     --model_name::String - one of "hidden", "visible", or "expectedvalue"
#     --risk_measure::String - one of "Expectation", "ModifiedChiSquared(radius)"
#
# Example usage:
#     julia --project=. examples/belief.jl --num_replications=500 \
#         --time_limit=20 --model_name=expectedvalue \
#         --risk_measure=ModifiedChiSquared(0.5)
# ==============================================================================
if length(ARGS) > 0
    using Gurobi
    optimizer = with_optimizer(Gurobi.Optimizer, OutputFlag = 0)
    # Parse the input arguments.
    args = Dict{String, String}(
        "num_replications" => "2000",
        "time_limit" => "20",
        "model_name" => "expectedvalue",
        "risk_measure" => "Expectation"
    )
    for arg in ARGS
        if startswith(arg, "--") && occursin("=", arg)
            items = split(arg[3:end], "=")
            args[string(items[1])] = items[2]
        else
            @warn("Argument $(arg) not supported, ignoring it.")
        end
    end
    num_replications = parse(Int, args["num_replications"])
    time_limit = parse(Int, args["time_limit"])
    model_name = args["model_name"]
    risk_measure = if args["risk_measure"] == "Expectation"
        Kokako.Expectation()
    elseif args["risk_measure"] == "WorstCase"
        Kokako.WorstCase()
    elseif startswith(args["risk_measure"], "ModifiedChiSquared")
        m = match(r"\(([\d\.]+)\)", args["risk_measure"])
        if m === nothing
            error("Unable to parse radius of ModifiedChiSquared: $(args["risk_measure"])")
        end
        radius = parse(Float64, m[1])
        Kokako.ModifiedChiSquared(radius)
    end

    # Initialize noise terms.
    demand_values = [1.0, 2.0]
    demand_probs = Dict(
        :Ah => [0.2, 0.8],
        :Bh => [0.8, 0.2],
        :H => [0.5, 0.5]
    )
    # The default graph.
    graph = Kokako.Graph(
        :root_node,
        [:Ad, :Ah, :Bd, :Bh],
        [
            (:root_node => :Ad, 0.5), (:root_node => :Bd, 0.5),
            (:Ad => :Ah, 1.0), (:Ah => :Ad, 0.9),
            (:Bd => :Bh, 1.0), (:Bh => :Bd, 0.9)
        ]
    )
    Ad, Ah, Bd, Bh = :Ad, :Ah, :Bd, :Bh

    if model_name == "hidden"
        Kokako.add_partition(graph, [:Ad, :Bd])
        Kokako.add_partition(graph, [:Ah, :Bh])
    elseif model_name == "visible"
        # Do nothing :)
    elseif model_name == "expectedvalue"
        graph = Kokako.Graph(
            :root_node,
            [:D, :H],
            [(:root_node => :D, 1.0), (:D => :H, 1.0), (:H => :D, 0.9)]
        )
        Ad, Ah, Bd, Bh = :D, :H, :D, :H
    else
        error("The value of --model_name=$(model_name) is not supported.")
    end

    # Build the model.
    model = inventory_management_problem(
        graph, demand_values, demand_probs, optimizer
    )

    # Train the policy.
    Random.seed!(123)
    Kokako.train(model;
        risk_measure = (node) -> (node == :Ad || node == :Bd || node == :D) ? Kokako.Expectation() : risk_measure,
        time_limit = time_limit,
        print_level = 1
    )

    # Simulate policy using same realizations for each policy.
    Random.seed!(1234)
    simulations = Any[]
    for replication in 1:num_replications
        actual_node = rand() < 0.5 ? :A : :B
        num_stages = ceil(Int, log(rand()) / log(0.9))
        scenario = Tuple{Symbol, Union{Nothing, Float64}}[]
        for stage in 1:num_stages
            push!(scenario, ((actual_node == :A ? Ad : Bd), nothing))

            # Form second stage.
            node = (actual_node == :A ? Ah : Bh)
            r = rand()
            noise = last(demand_probs[node])
            for (idx, prob) in enumerate(demand_probs[node])
                r -= prob
                if r < 0.0
                    noise = demand_values[idx]
                    break
                end
            end
            push!(scenario, ((actual_node == :A ? Ah : Bh), noise))
        end
        push!(simulations,
            Kokako.simulate(model, 1, [:inventory, :buy], sampling_scheme=Kokako.Historical(scenario))[1]
        )
    end

    # Calculate some nice statistics.
    objectives = map(sim -> sum(s[:stage_objective] for s in sim), simulations)
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(
        1.96 * Statistics.std(objectives) / sqrt(num_replications); digits = 2)
    println("Normal assumption: μ = $(sample_mean) ± $(sample_ci)")

    λ = 1 / mean(objectives)
    half_width = 1.96 / sqrt(length(objectives))
    lower = round(1 / (λ * (1 + half_width)); digits = 2)
    upper = round(1 / (λ * (1 - half_width)); digits = 2)
    println("Exponential assumption: μ ∈ [$(lower), $(upper)]")
    # Save for later analysis.
    open("$(model_name).json", "w") do io
        write(io, Kokako.JSON.json(simulations))
    end

    if model_name == "hidden"
        get_hidden_value_function(model)
    end
else  # Test the example!
    demand_values = [1.0, 2.0]
    demand_probs = Dict(
        :Ah => [0.2, 0.8],
        :Bh => [0.8, 0.2],
        :H => [0.5, 0.5]
    )
    graph = Kokako.Graph(
        :root_node,
        [:Ad, :Ah, :Bd, :Bh],
        [
            (:root_node => :Ad, 0.5), (:root_node => :Bd, 0.5),
            (:Ad => :Ah, 1.0), (:Ah => :Ad, 0.9),
            (:Bd => :Bh, 1.0), (:Bh => :Bd, 0.9)
        ]
    )
    Kokako.add_partition(graph, [:Ad, :Bd])
    Kokako.add_partition(graph, [:Ah, :Bh])

    # Build the model.
    model = inventory_management_problem(graph, demand_values, demand_probs)
    # Train the policy.
    Random.seed!(123)
    Kokako.train(model; iteration_limit = 100, print_level = 1)
    results = Kokako.simulate(model, 500)
    objectives = [
        sum(s[:stage_objective] for s in simulation) for simulation in results
    ]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    println("Confidence_interval = $(sample_mean) ± $(sample_ci)")
    @test Kokako.calculate_bound(model) ≈ sample_mean atol = sample_ci
end
