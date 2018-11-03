using Kokako, GLPK, Test, Statistics

function infinite_hydro_thermal()
    Ω = [
        (inflow = 0.0, demand = 7.5),
        (inflow = 5.0, demand = 5),
        (inflow = 10.0, demand = 2.5)
    ]
    graph = Kokako.Graph(
        :root_node,
        [:week],
        [(:root_node => :week, 1.0), (:week => :week, 0.9)]
    )
    model = Kokako.PolicyGraph(graph,
                bellman_function = Kokako.AverageCut(lower_bound = 0),
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do subproblem, node
        @variable(subproblem,
            5.0 <= reservoir <= 15.0, Kokako.State, initial_value = 10.0)
        @variables(subproblem, begin
            thermal_generation >= 0
            hydro_generation >= 0
            spill >= 0
            inflow
            demand
        end)
        @constraints(subproblem, begin
            reservoir.out == reservoir.in - hydro_generation - spill + inflow
            hydro_generation + thermal_generation == demand
        end)
        @stageobjective(subproblem, 10 * spill + thermal_generation)
        Kokako.parameterize(subproblem, Ω) do ω
            JuMP.fix(inflow, ω.inflow)
            JuMP.fix(demand, ω.demand)
        end
    end
    Kokako.train(model, print_level = 1, time_limit = 5.0)
    @test Kokako.calculate_bound(model) ≈ 119.167 atol = 5.0

    results = Kokako.simulate(model, 500)
    objectives = [
        sum(s[:stage_objective] for s in simulation) for simulation in results
    ]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    println("Confidence_interval = $(sample_mean) ± $(sample_ci)")
    @test sample_mean ≈ 119.167 atol = 5.0
end

infinite_hydro_thermal()
