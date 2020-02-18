#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using GLPK
using SDDP
using Test

@testset "Read/write cuts to file" begin
    function create_model(graph)
        return SDDP.PolicyGraph(
            graph,
            lower_bound = 0.0,
            optimizer = GLPK.Optimizer,
        ) do subproblem, node
            @variable(
                subproblem,
                5.0 <= reservoir <= 15.0,
                SDDP.State,
                initial_value = 10.0
            )
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
            SDDP.parameterize(
                subproblem,
                [
                    (inflow = 0.0, demand = 7.5),
                    (inflow = 5.0, demand = 5),
                    (inflow = 10.0, demand = 2.5),
                ],
            ) do ω
                JuMP.fix(inflow, ω.inflow)
                JuMP.fix(demand, ω.demand)
            end
        end
    end

    graphs = [
        (
            "Symbol",
            SDDP.Graph(
                :root_node,
                [:week],
                [(:root_node => :week, 1.0), (:week => :week, 0.9)],
            ),
        ),
        ("Int", SDDP.Graph(0, [1], [(0 => 1, 1.0), (1 => 1, 0.9)])),
        (
            "NTuple",
            SDDP.Graph(
                (0, 1),
                [(1, 1)],
                [((0, 1) => (1, 1), 1.0), ((1, 1) => (1, 1), 0.9)],
            ),
        ),
    ]
    for (T, graph) in graphs
        @testset "$T" begin
            model = create_model(graph)
            @test SDDP.calculate_bound(model) ≈ 9.17 atol = 0.1
            SDDP.train(model; iteration_limit = 50, print_level = 0)
            @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1
            SDDP.write_cuts_to_file(model, "$(T).cuts.json")
            model_2 = create_model(graph)
            @test SDDP.calculate_bound(model_2) ≈ 9.17 atol = 0.1
            SDDP.read_cuts_from_file(model_2, "$(T).cuts.json")
            @test SDDP.calculate_bound(model_2) ≈ 119.167 atol = 0.1
            rm("$(T).cuts.json")
        end
    end
    @testset "String" begin
        graph = SDDP.Graph(
            "root_node",
            ["week"],
            [("root_node" => "week", 1.0), ("week" => "week", 0.9)],
        )
        model = create_model(graph)
        @test SDDP.calculate_bound(model) ≈ 9.17 atol = 0.1
        SDDP.train(model; iteration_limit = 50, print_level = 0, cut_type = SDDP.MULTI_CUT)
        @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1
        SDDP.write_cuts_to_file(model, "model.cuts.json")
        model_2 = create_model(graph)
        @test SDDP.calculate_bound(model_2) ≈ 9.17 atol = 0.1
        @test_throws Exception SDDP.read_cuts_from_file(model_2, "model.cuts.json")
        SDDP.read_cuts_from_file(
            model_2,
            "model.cuts.json",
            node_name_parser = (::Type{String}, x::String) -> x,
        )
        @test SDDP.calculate_bound(model_2) ≈ 119.167 atol = 0.1
        rm("model.cuts.json")
    end
end
