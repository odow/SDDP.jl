using Kokako, GLPK, Test

function infinite_hydro_thermal()
    Ω = [
        (inflow=0.0, demand=7.5),
        (inflow=5.0, demand=5),
        (inflow=10.0, demand=2.5)
    ]
    graph = Kokako.Graph(
        :root_node,
        [:week],
        [(:root_node => :week, 1.0), (:week => :week, 0.9)]
    )
    model = Kokako.PolicyGraph(graph,
                bellman_function = Kokako.AverageCut(lower_bound=0),
                optimizer = with_optimizer(GLPK.Optimizer)
                    ) do subproblem, node
        @variables(subproblem, begin
            incoming_reservoir
            5.0 <= outgoing_reservoir <= 15.0
            thermal_generation >= 0
            hydro_generation >= 0
            spill >= 0
            inflow
            demand
        end)
        Kokako.add_state_variable(subproblem, incoming_reservoir, outgoing_reservoir, 10.0)
        @constraints(subproblem, begin
            outgoing_reservoir == incoming_reservoir - hydro_generation - spill + inflow
            hydro_generation + thermal_generation == demand
        end)
        @stageobjective(subproblem, Min, 10 * spill + thermal_generation)
        Kokako.parameterize(subproblem, Ω) do ω
            JuMP.fix(inflow, ω.inflow)
            JuMP.fix(demand, ω.demand)
        end
    end
    Kokako.train(model, iteration_limit=100, print_level=1)
end

infinite_hydro_thermal()
