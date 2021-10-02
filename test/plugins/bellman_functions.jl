#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestBellmanFunctions

using SDDP
using Test
import GLPK

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function _create_model(graph)
    return SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do subproblem, _
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
                reservoir.out ==
                reservoir.in - hydro_generation - spill + inflow
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
            return JuMP.fix(demand, ω.demand)
        end
    end
end

function test_Read_write_cuts_to_file()
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
        model = _create_model(graph)
        @test SDDP.calculate_bound(model) ≈ 9.17 atol = 0.1
        SDDP.train(model; iteration_limit = 50, print_level = 0)
        @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1
        SDDP.write_cuts_to_file(model, "$(T).cuts.json")
        model_2 = _create_model(graph)
        @test SDDP.calculate_bound(model_2) ≈ 9.17 atol = 0.1
        SDDP.read_cuts_from_file(model_2, "$(T).cuts.json")
        @test SDDP.calculate_bound(model_2) ≈ 119.167 atol = 0.1
        rm("$(T).cuts.json")
    end
    return
end

function test_read_write_cuts_to_file_String()
    graph = SDDP.Graph(
        "root_node",
        ["week"],
        [("root_node" => "week", 1.0), ("week" => "week", 0.9)],
    )
    model = _create_model(graph)
    @test SDDP.calculate_bound(model) ≈ 9.17 atol = 0.1
    SDDP.train(
        model;
        iteration_limit = 50,
        print_level = 0,
        cut_type = SDDP.MULTI_CUT,
    )
    @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1
    SDDP.write_cuts_to_file(model, "model.cuts.json")
    model_2 = _create_model(graph)
    @test SDDP.calculate_bound(model_2) ≈ 9.17 atol = 0.1
    @test_throws Exception SDDP.read_cuts_from_file(model_2, "model.cuts.json")
    SDDP.read_cuts_from_file(
        model_2,
        "model.cuts.json",
        node_name_parser = (::Type{String}, x::String) -> x,
    )
    @test SDDP.calculate_bound(model_2) ≈ 119.167 atol = 0.1
    rm("model.cuts.json")
    return
end

function test_add_all_cuts_SINGLE_CUT()
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, 5 <= x <= 15, SDDP.State, initial_value = 10)
        @variable(sp, g >= 0)
        @variable(sp, h >= 0)
        @variable(sp, u >= 0)
        @constraint(sp, inflow, x.out == x.in - h - u)
        @constraint(sp, demand, h + g == 0)
        @stageobjective(sp, 10 * u + g)
        SDDP.parameterize(sp, [(0, 7.5), (5, 5.0), (10, 2.5)]) do ω
            set_normalized_rhs(inflow, ω[1])
            set_normalized_rhs(demand, ω[2])
            return
        end
    end
    SDDP.train(model; iteration_limit = 10)
    for (t, node) in model.nodes
        @test num_constraints(
            node.subproblem,
            AffExpr,
            MOI.GreaterThan{Float64},
        ) < 10
    end
    SDDP.add_all_cuts(model)
    for (t, node) in model.nodes
        n = num_constraints(node.subproblem, AffExpr, MOI.GreaterThan{Float64})
        @test t == 3 || n == 10
    end
    return
end

function test_add_all_cuts_MULTI_CUT()
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = GLPK.Optimizer,
    ) do sp, t
        @variable(sp, 5 <= x <= 15, SDDP.State, initial_value = 10)
        @variable(sp, g >= 0)
        @variable(sp, h >= 0)
        @variable(sp, u >= 0)
        @constraint(sp, inflow, x.out == x.in - h - u)
        @constraint(sp, demand, h + g == 0)
        @stageobjective(sp, 10 * u + g)
        SDDP.parameterize(sp, [(0, 7.5), (5, 5.0), (10, 2.5)]) do ω
            set_normalized_rhs(inflow, ω[1])
            set_normalized_rhs(demand, ω[2])
            return
        end
    end
    SDDP.train(model; iteration_limit = 10, cut_type = SDDP.MULTI_CUT)
    for (t, node) in model.nodes
        @test num_constraints(
            node.subproblem,
            AffExpr,
            MOI.GreaterThan{Float64},
        ) < 31
    end
    SDDP.add_all_cuts(model)
    for (t, node) in model.nodes
        n = num_constraints(node.subproblem, AffExpr, MOI.GreaterThan{Float64})
        @test t == 3 || n == 31
    end
    return
end

end  # module

TestBellmanFunctions.runtests()
