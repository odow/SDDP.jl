#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Test, GLPK

@testset "Basic Graphs" begin
    @testset "LinearGraph" begin
        graph = SDDP.LinearGraph(5)
        @test graph.root_node == 0
        for stage in 0:4
            @test haskey(graph.nodes, stage)
            @test graph.nodes[stage] == [(stage + 1, 1.0)]
        end
        @test haskey(graph.nodes, 5)
        @test graph.nodes[5] == Tuple{Int, Float64}[]

        graph = SDDP.LinearGraph(3)
        @test sprint(show, graph) ==
            "Root\n" *
            " 0\n" *
            "Nodes\n" *
            " 1\n" *
            " 2\n" *
            " 3\n" *
            "Arcs\n" *
            " 0 => 1 w.p. 1.0\n" *
            " 1 => 2 w.p. 1.0\n" *
            " 2 => 3 w.p. 1.0\n"
    end

    @testset "MarkovianGraph" begin
        @testset "Error cases" begin
            # Not root transition matrix.
            @test_throws Exception SDDP.MarkovianGraph([[0.5 0.5; 0.5 0.5]])
            # Negative probability.
            @test_throws Exception SDDP.MarkovianGraph([[-0.5 0.75]])
            # Proability sums to greater than 1.
            @test_throws Exception SDDP.MarkovianGraph([[0.8 0.8]])
            # Mis-matched dimensions.
            @test_throws Exception SDDP.MarkovianGraph([
                [0.1 0.2 0.7], [0.5 0.5; 0.5 0.5]
            ])
        end
        @testset "keyword vs list" begin
            graph_1 = SDDP.MarkovianGraph(
                stages = 2,
                transition_matrix = [0.4 0.6; 0.25 0.75],
                root_node_transition = [0.7, 0.3]
            )
            graph_2 = SDDP.MarkovianGraph([
                [0.7 0.3], [0.4 0.6; 0.25 0.75]
            ])
            @test graph_1.root_node == graph_2.root_node
            @test graph_1.nodes == graph_2.nodes
        end
    end

    @testset "Graph" begin
        @testset "Construct Graph" begin
            graph = SDDP.Graph(:root)
            @test graph.root_node == :root
            @test collect(keys(graph.nodes)) == [:root]
        end
        @testset "Add node" begin
            graph = SDDP.Graph(:root)
            SDDP.add_node(graph, :x)
            @test collect(keys(graph.nodes)) == [:root, :x]
        end
        @testset "Add duplicate node" begin
            graph = SDDP.Graph(:root)
            SDDP.add_node(graph, :x)
            @test_throws Exception SDDP.add_node(graph, :x)
        end
        @testset "Add edge" begin
            graph = SDDP.Graph(:root)
            SDDP.add_node(graph, :x)
            SDDP.add_edge(graph, :root => :x, 1.0)
            @test haskey(graph.nodes, :root)
            @test graph.nodes[:root] == [(:x, 1.0)]
        end
        @testset "Add edge of wrong type" begin
            graph = SDDP.Graph(:root)
            @test_throws Exception SDDP.add_node(graph, 1)
        end
        @testset "Add edge to missing node" begin
            graph = SDDP.Graph(:root)
            SDDP.add_node(graph, :x)
            @test_throws Exception SDDP.add_edge(graph, :x => :y, 1.0)
            @test_throws Exception SDDP.add_edge(graph, :y => :x, 1.0)
        end
        @testset "Add edge to root" begin
            graph = SDDP.Graph(:root)
            SDDP.add_node(graph, :x)
            @test_throws Exception SDDP.add_edge(graph, :x => :root, 1.0)
        end
    end
end

@testset "PolicyGraph constructor" begin
    @testset "LinearGraph" begin
        model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                                   lower_bound = 0.0,
                                   direct_mode=false) do node, stage
        end

        @test_throws Exception SDDP.PolicyGraph(SDDP.LinearGraph(2),
                                                  lower_bound = 0.0
                                                      ) do node, stage
        end
        nodes = Set{Int}()
        model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                                   lower_bound = 0.0,
                                   optimizer = with_optimizer(GLPK.Optimizer)
                                       ) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([1, 2])
        @test sprint(show, model) ==
            "A policy graph with 2 nodes.\n Node indices: 1, 2\n"
    end

    @testset "MarkovianGraph" begin
        graph = SDDP.MarkovianGraph([
                ones(Float64, (1, 1)),
                [0.5 0.5],
                [0.5 0.5; 0.3 0.4],
                [0.5 0.5; 0.3 0.4],
                [0.5 0.5; 0.3 0.4]
            ]
        )
        nodes = Set{Tuple{Int, Int}}()
        model = SDDP.PolicyGraph(graph, lower_bound = 0.0,
                                   direct_mode = false) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (5, 2)])
    end

    @testset "MarkovianPolicyGraph" begin
        nodes = Set{Tuple{Int, Int}}()
        model = SDDP.MarkovianPolicyGraph(
                transition_matrices = [
                    ones(Float64, (1, 1)),
                    [0.5 0.5],
                    [0.5 0.5; 0.3 0.4],
                    [0.5 0.5; 0.3 0.4],
                    [0.5 0.5; 0.3 0.4]],
                lower_bound = 0.0, direct_mode = false) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (5, 2)])
    end

    @testset "General" begin
        graph = SDDP.Graph(
            :root,
            [:stage_1, :stage_2, :stage_3],
            [
                (:root => :stage_1, 1.0),
                (:stage_1 => :stage_2, 1.0),
                (:stage_2 => :stage_3, 1.0),
                (:stage_3 => :stage_1, 0.9)
            ]
        )
        nodes = Set{Symbol}()
        model = SDDP.PolicyGraph(graph, lower_bound = 0.0,
                                   direct_mode = false) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([:stage_1, :stage_2, :stage_3])
    end
end

@testset "SDDP.State" begin
    model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                               lower_bound = 0.0,
                               direct_mode = false) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0)
    end
    for stage in 1:2
        node = model[stage]
        @test haskey(node.states, :x)
        @test length(keys(node.states)) == 1
        @test node.states[:x] == node.subproblem[:x]
    end
end

@testset "SDDP.parameterize" begin
    model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                               lower_bound = 0.0,
                               direct_mode = false) do node, stage
        @variable(node, 0 <= x <= 1)
        SDDP.parameterize(node, [1, 2, 3], [0.4, 0.5, 0.1]) do ω
            JuMP.set_upper_bound(x, ω)
        end
    end
    node = model[2]
    @test length(node.noise_terms) == 3
    @test JuMP.upper_bound(node.subproblem[:x]) == 1
    node.parameterize(node.noise_terms[2].term)
    @test JuMP.upper_bound(node.subproblem[:x]) == 2
    node.parameterize(3)
    @test JuMP.upper_bound(node.subproblem[:x]) == 3
end

@testset "SDDP.set_stage_objective" begin
    @testset ":Min" begin
        model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                                   lower_bound = 0.0,
                                   direct_mode = false) do node, stage
            @variable(node, 0 <= x <= 1)
            @stageobjective(node, 2x)
        end
        node = model[2]
        @test node.stage_objective == 2 * node.subproblem[:x]
        @test model.objective_sense == SDDP.MOI.MIN_SENSE
    end

    @testset ":Max" begin
        model = SDDP.PolicyGraph(SDDP.LinearGraph(2),
                                   lower_bound = 0.0,
                                   sense = :Max,
                                   direct_mode = false) do node, stage
            @variable(node, 0 <= x <= 1)
            @stageobjective(node, 2x)
        end
        node = model[2]
        @test node.stage_objective == 2 * node.subproblem[:x]
        @test model.objective_sense == SDDP.MOI.MAX_SENSE
    end
end

@testset "Errors" begin
    @testset "<=0 stages" begin
        exception = ErrorException(
            "You must create a LinearPolicyGraph with `stages >= 1`.")
        @test_throws exception SDDP.LinearPolicyGraph(stages = 0) do sp, t
        end
    end
    @testset "missing bounds" begin
        exception = ErrorException(
            "You must specify a bound on the objective value, through " *
            "`lower_bound` if minimizing, or `upper_bound` if maximizing.")
        @test_throws exception SDDP.LinearPolicyGraph(stages = 1) do sp, t
        end
    end
    @testset "parameterize!" begin
        exception = ErrorException(
            "Duplicate calls to SDDP.parameterize detected.")
        @test_throws exception SDDP.LinearPolicyGraph(
                stages = 2, lower_bound = 0.0, sense = :Max, direct_mode = false
                ) do node, stage
            @variable(node, 0 <= x <= 1)
            SDDP.parameterize(node, [1, 2]) do ω
                @stageobjective(node, ω * x)
            end
            SDDP.parameterize(node, [3, 4]) do ω
                @stageobjective(node, ω * x)
            end
        end
    end
    @testset "no initial_value" begin
        exception = ErrorException(
            "In `@variable(node, x, SDDP.State)`: When creating a state " *
            "variable, you must set the `initial_value` keyword to the value " *
            "of the state variable at the root node.")
        @test_throws exception SDDP.LinearPolicyGraph(
                stages = 2, lower_bound = 0.0, sense = :Max, direct_mode = false
                ) do node, stage
            @variable(node, x, SDDP.State)
            @stageobjective(node, x.out)
        end
    end
end

@testset "Numerical stability report" begin
    model = SDDP.LinearPolicyGraph(
            stages = 2, lower_bound = -1e10, direct_mode=false) do subproblem, t
        @variable(subproblem, x >= -1e7, SDDP.State, initial_value=1e-5)
        @constraint(subproblem, 1e9 * x.out >= 1e-6 * x.in + 1e-8)
        @stageobjective(subproblem, 1e9 * x.out)
    end
    report = sprint(SDDP.numerical_stability_report, model)
    @test occursin("WARNING", report)
    report_2 = sprint(
        io -> SDDP.numerical_stability_report(io, model, by_node=true))
    @test occursin("Numerical stability report for node: 1", report_2)
    @test occursin("Numerical stability report for node: 2", report_2)
end
