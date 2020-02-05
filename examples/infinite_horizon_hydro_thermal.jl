#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, GLPK, Test, Statistics

function infinite_hydro_thermal(; cut_type)
    Ω = [
        (inflow = 0.0, demand = 7.5),
        (inflow = 5.0, demand = 5),
        (inflow = 10.0, demand = 2.5),
    ]
    graph =
        SDDP.Graph(:root_node, [:week], [(:root_node => :week, 1.0), (:week => :week, 0.9)])
    model = SDDP.PolicyGraph(
        graph,
        bellman_function = SDDP.BellmanFunction(lower_bound = 0, cut_type = cut_type),
        optimizer = GLPK.Optimizer,
    ) do subproblem, node
        @variable(subproblem, 5.0 <= reservoir <= 15.0, SDDP.State, initial_value = 10.0)
        @variables(subproblem, begin
            thermal_generation >= 0
            hydro_generation >= 0
            spill >= 0
            inflow
            demand
        end)
        @constraints(
            subproblem,
            begin
                reservoir.out == reservoir.in - hydro_generation - spill + inflow
                hydro_generation + thermal_generation == demand
            end
        )
        @stageobjective(subproblem, 10 * spill + thermal_generation)
        SDDP.parameterize(subproblem, Ω) do ω
            JuMP.fix(inflow, ω.inflow)
            JuMP.fix(demand, ω.demand)
        end
    end
    SDDP.train(
        model;
        time_limit = 2.0,
        print_level = 0,
        sampling_scheme = SDDP.InSampleMonteCarlo(terminate_on_cycle = true),
        cycle_discretization_delta = 0.1,
    )
    @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1

    results = SDDP.simulate(model, 500)
    objectives = [sum(s[:stage_objective] for s in simulation) for simulation in results]
    sample_mean = round(Statistics.mean(objectives); digits = 2)
    sample_ci = round(1.96 * Statistics.std(objectives) / sqrt(500); digits = 2)
    println("Confidence_interval = $(sample_mean) ± $(sample_ci)")
    @test sample_mean - sample_ci <= 119.167 <= sample_mean + sample_ci
end

infinite_hydro_thermal(cut_type = SDDP.SINGLE_CUT)
infinite_hydro_thermal(cut_type = SDDP.MULTI_CUT)
