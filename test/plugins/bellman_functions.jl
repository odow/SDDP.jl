#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestBellmanFunctions

using SDDP
using Test
import HiGHS

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
        optimizer = HiGHS.Optimizer,
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

function test_read_write_cuts_to_file_ValueFunction()
    graph = SDDP.Graph(
        "root_node",
        ["week"],
        [("root_node" => "week", 1.0), ("week" => "week", 0.9)],
    )
    model = _create_model(graph)
    SDDP.train(model; iteration_limit = 50, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1
    V = SDDP.ValueFunction(model, node = "week")
    value_f = SDDP.evaluate(V, reservoir = 10)
    SDDP.write_cuts_to_file(model, "model.cuts.json")
    model_2 = _create_model(graph)
    @test SDDP.calculate_bound(model_2) ≈ 9.17 atol = 0.1
    SDDP.read_cuts_from_file(
        model_2,
        "model.cuts.json",
        node_name_parser = (::Type{String}, x::String) -> x,
    )
    V2 = SDDP.ValueFunction(model_2, node = "week")
    @test value_f == SDDP.evaluate(V2, reservoir = 10)
    rm("model.cuts.json")
    return
end

function test_read_read_cuts_from_file_nothing()
    graph = SDDP.Graph(
        "root_node",
        ["week"],
        [("root_node" => "week", 1.0), ("week" => "week", 0.9)],
    )
    model = _create_model(graph)
    SDDP.train(model; iteration_limit = 50, print_level = 0)
    @test SDDP.calculate_bound(model) ≈ 119.167 atol = 0.1
    V = SDDP.ValueFunction(model, node = "week")
    value_f = SDDP.evaluate(V, reservoir = 10)
    SDDP.write_cuts_to_file(
        model,
        "model.cuts.json";
        node_name_parser = s -> "myname_$s",
    )
    model_2 = _create_model(graph)
    @test SDDP.calculate_bound(model_2) ≈ 9.17 atol = 0.1
    function parser(::Type{String}, x::String)
        @test startswith(x, "myname_")
        return replace(x, "myname_" => "")
    end
    SDDP.read_cuts_from_file(
        model_2,
        "model.cuts.json",
        node_name_parser = parser,
    )
    N = num_constraints(
        model_2["week"].subproblem;
        count_variable_in_set_constraints = true,
    )
    SDDP.read_cuts_from_file(
        model_2,
        "model.cuts.json",
        node_name_parser = (::Any, s) -> nothing,
    )
    N2 = num_constraints(
        model_2["week"].subproblem;
        count_variable_in_set_constraints = true,
    )
    @test N == N2
    rm("model.cuts.json")
    return
end

function test_add_all_cuts_SINGLE_CUT()
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
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
        optimizer = HiGHS.Optimizer,
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

function test_belief_state_cut_selection()
    demand_values = [1.0, 2.0]
    demand_prob = Dict(:Ah => [0.2, 0.8], :Bh => [0.8, 0.2])
    graph = SDDP.Graph(
        :root_node,
        [:Ad, :Ah, :Bd, :Bh],
        [
            (:root_node => :Ad, 0.5),
            (:root_node => :Bd, 0.5),
            (:Ad => :Ah, 1.0),
            (:Bd => :Bh, 1.0),
        ],
    )
    SDDP.add_ambiguity_set(graph, [:Ad, :Bd], 1e2)
    SDDP.add_ambiguity_set(graph, [:Ah, :Bh], 1e2)
    model = SDDP.PolicyGraph(
        graph,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do subproblem, node
        @variables(
            subproblem,
            begin
                0 <= inventory <= 2, (SDDP.State, initial_value = 0.0)
                buy >= 0
                demand
            end
        )
        @constraint(subproblem, demand == inventory.in - inventory.out + buy)
        if node == :Ad || node == :Bd || node == :D
            JuMP.fix(demand, 0)
            @stageobjective(subproblem, buy)
        else
            SDDP.parameterize(subproblem, demand_values, demand_prob[node]) do ω
                return JuMP.fix(demand, ω)
            end
            @stageobjective(subproblem, 2 * buy + inventory.out)
        end
    end
    SDDP.train(
        model,
        iteration_limit = 20,
        cut_deletion_minimum = 20,
        print_level = 0,
    )
    n_cuts = count(model[:Ad].bellman_function.global_theta.cuts) do cut
        return cut.constraint_ref !== nothing
    end
    @test n_cuts == 20
    SDDP.train(
        model,
        iteration_limit = 1,
        add_to_existing_cuts = true,
        print_level = 0,
    )
    n_cuts_2 = count(model[:Ad].bellman_function.global_theta.cuts) do cut
        return cut.constraint_ref !== nothing
    end
    @test n_cuts_2 < n_cuts
    return
end

function test_biobjective_cut_selection()
    model = SDDP.LinearPolicyGraph(
        stages = 3,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    ) do subproblem, _
        @variable(subproblem, 0 <= v <= 200, SDDP.State, initial_value = 50)
        @variables(subproblem, begin
            0 <= g[i = 1:2] <= 100
            0 <= u <= 150
            s >= 0
            shortage_cost >= 0
        end)
        @expressions(subproblem, begin
            objective_1, g[1] + 10 * g[2]
            objective_2, shortage_cost
        end)
        @constraints(subproblem, begin
                inflow_constraint, v.out == v.in - u - s
                g[1] + g[2] + u == 150
                shortage_cost >= 40 - v.out
                shortage_cost >= 60 - 2 * v.out
                shortage_cost >= 80 - 4 * v.out
            end)
        ## You must call this for a biobjective problem!
        SDDP.initialize_biobjective_subproblem(subproblem)
        SDDP.parameterize(subproblem, 0.0:5:50.0) do ω
            JuMP.set_normalized_rhs(inflow_constraint, ω)
            ## You must call `set_biobjective_functions` from within
            ## `SDDP.parameterize`.
            return SDDP.set_biobjective_functions(
                subproblem,
                objective_1,
                objective_2,
            )
        end
    end
    SDDP.train_biobjective(
        model,
        solution_limit = 10,
        iteration_limit = 10,
        print_level = 0,
    )
    n_cuts = count(model[1].bellman_function.global_theta.cuts) do cut
        return cut.constraint_ref !== nothing
    end
    @test n_cuts < 100
    return
end

end  # module

TestBellmanFunctions.runtests()
