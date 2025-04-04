#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestInnerBellmanFunctions

using SDDP
using Test
import HiGHS

import SDDP: Inner

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

function build(subproblem, t)
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
        return JuMP.fix(demand, ω.demand)
    end
end

function _create_model(graph, bellman_function = nothing)
    return SDDP.PolicyGraph(
        build,
        graph;
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
        bellman_function,
    )
end

function test_Read_write_cuts_to_file()
    nstages = 4
    graphs = [
        (
            "Int",
            SDDP.Graph(
                0,
                [1, 2, 3, 4],
                [(0 => 1, 1.0), (1 => 2, 1.0), (2 => 3, 1.0), (3 => 4, 1.0)],
            ),
        ),
        (
            "NTuple",
            SDDP.Graph(
                (0, 1),
                [(1, 1), (1, 2), (1, 3), (1, 4)],
                [
                    ((0, 1) => (1, 1), 1.0),
                    ((1, 1) => (1, 2), 1.0),
                    ((1, 2) => (1, 3), 1.0),
                    ((1, 3) => (1, 4), 1.0),
                ],
            ),
        ),
        (
            "Symbol",
            SDDP.Graph(
                :root_node,
                [:w1, :w2, :w3, :w4],
                [
                    (:root_node => :w1, 1.0),
                    (:w1 => :w2, 1.0),
                    (:w2 => :w3, 1.0),
                    (:w3 => :w4, 1.0),
                ],
            ),
        ),
    ]
    # Upper bound via inner approximation
    base_Lip = (10.0 + 1.0)
    base_ub = (10.0 + 7.5) * base_Lip
    build_Lip(t) = base_Lip * (nstages - t)
    build_ub(t) = base_ub * (nstages - t)
    ibf = Inner.InnerBellmanFunction(
        build_Lip;
        upper_bound = build_ub,
        vertex_type = SDDP.SINGLE_CUT,
    )

    # TODO: generalize other graphs, esp. regarding types
    for (T, graph) in graphs[1:1]
        model = _create_model(graph)
        @test SDDP.calculate_bound(model) ≈ 9.17 atol = 0.1
        SDDP.train(model; iteration_limit = 50, print_level = 0)
        @test SDDP.calculate_bound(model) ≈ 45.833 atol = 0.1
        model_inner, ub = Inner.inner_dp(
            build,
            model;
            nstages,
            sense = :Min,
            optimizer = HiGHS.Optimizer,
            lower_bound = 0.0,
            bellman_function = ibf,
            risk_measures = SDDP.Expectation(),
            print_level = 1,
        )
        @test ub ≈ 45.833 atol = 0.1
        @test SDDP.calculate_bound(model_inner) ≈ ub
        SDDP.Inner.write_vertices_to_file(model_inner, "$(T).vertices.json")

        model_2 = _create_model(graph, ibf)
        for (k, node) in model_2.nodes
            SDDP.set_objective(node)
        end
        @test SDDP.calculate_bound(model_2) ≈ 801.666 atol = 0.1
        @test_throws JuMP.NoOptimizer SDDP.Inner.read_vertices_from_file(
            model_2,
            "$(T).vertices.json",
            vertex_selection = true,
        )
        SDDP.Inner.read_vertices_from_file(
            model_2,
            "$(T).vertices.json";
            vertex_selection = true,
            optimizer = HiGHS.Optimizer,
        )
        @test SDDP.calculate_bound(model_2) ≈ ub atol = 0.1
        rm("$(T).vertices.json")
    end
    return
end

function create_policy_graph_with_inner_approximation()
    nstages = 20
    graph = SDDP.LinearGraph(nstages)

    # for some reason, creating the model with the same
    # lipschitz estimates used in the test_Read_write_cuts_to_file
    # makes the test for calculating bound after training fails

    return SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        lower_bound = 0.0,
        upper_bound = 1000,
        optimizer = HiGHS.Optimizer,
        lipschitz_constant = 10.0,
    )
end

function test_InnerPolicyGraph_LinearGraph()
    nstages = 4
    graph = SDDP.LinearGraph(nstages)

    # missing lipschitz_constant
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        sense = :Min,
        lower_bound = 0.0,
        upper_bound = 1000.0,
        optimizer = HiGHS.Optimizer,
    )

    # missing lower_bound
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        sense = :Min,
        upper_bound = 1000.0,
        optimizer = HiGHS.Optimizer,
        lipschitz_constant = 10.0,
    )

    # missing upper_bound
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        sense = :Max,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
        lipschitz_constant = 10.0,
    )

    # missing upper_bound even when sense = :Min
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        sense = :Min,
        lower_bound = 0.0,
        optimizer = HiGHS.Optimizer,
        lipschitz_constant = 10.0,
    )

    # mising lower_bound even when sense = :Max
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        sense = :Max,
        upper_bound = 1000.0,
        optimizer = HiGHS.Optimizer,
        lipschitz_constant = 10.0,
    )
end

function test_InnerPolicyGraph_other_graphs()
    graph = SDDP.MarkovianGraph([
        ones(Float64, (1, 1)),
        [0.5 0.5],
        [0.5 0.5; 0.3 0.4],
        [0.5 0.5; 0.3 0.4],
        [0.5 0.5; 0.3 0.4],
    ])

    # rejects MarkovianGraph
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        lower_bound = 0.0,
        upper_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    )

    graph = SDDP.UnicyclicGraph(0.9; num_nodes = 4)

    # rejects UnicyclicGraph
    @test_throws Exception SDDP.Inner.InnerPolicyGraph(
        build,
        graph;
        lower_bound = 0.0,
        upper_bound = 0.0,
        optimizer = HiGHS.Optimizer,
    )
end

function test_InnerPolicyGraph_train()
    # calling SDDP.train with an "inner model" will probably result
    # in sub-optimal bounds, so this test checks for a bound with a
    # reasonable gap to the optimal (250.0 > 241.3833)
    #
    # however, the function call works and the inner bound might
    # decrease for some iterations
    model = create_policy_graph_with_inner_approximation()
    SDDP.train(model; iteration_limit = 200, print_level = 0)
    @test SDDP.calculate_bound(model) >= 250.0
end

function test_InnerPolicGraph_simulate()
    model = create_policy_graph_with_inner_approximation()
    SDDP.train(model; iteration_limit = 50, print_level = 0)
    results = SDDP.simulate(model, 50, [:vertex_coverage_distance])
    # TODO - make more meaningful tests with the simulation result
    @test typeof(results) === Vector{Vector{Dict{Symbol,Any}}}
end

function test_dp_vertices_from_visited_states()
    nstages = 20
    graph = SDDP.LinearGraph(nstages)
    cut_model = _create_model(graph)
    SDDP.train(cut_model; iteration_limit = 200, print_level = 0)
    vertex_model = create_policy_graph_with_inner_approximation()
    SDDP.Inner.dp_vertices_from_visited_states(
        vertex_model,
        cut_model;
        optimizer = HiGHS.Optimizer,
    )
    @test SDDP.calculate_bound(cut_model) ≈ 241.3833 atol = 0.1
    @test SDDP.calculate_bound(vertex_model) ≈ 241.3833 atol = 0.1
end

end  # module

TestInnerBellmanFunctions.runtests()
