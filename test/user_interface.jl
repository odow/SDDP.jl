#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Kokako, Test, GLPK

@testset "Basic Graphs" begin
    @testset "LinearGraph" begin
        graph = Kokako.LinearGraph(5)
        @test graph.root_node == 0
        for stage in 0:4
            @test haskey(graph.nodes, stage)
            @test graph.nodes[stage] == [(stage + 1, 1.0)]
        end
        @test haskey(graph.nodes, 5)
        @test graph.nodes[5] == Tuple{Int, Float64}[]
    end

    @testset "MarkovianGraph" begin
        @testset "Error cases" begin
            # Not root transition matrix.
            @test_throws Exception Kokako.MarkovianGraph([[0.5 0.5; 0.5 0.5]])
            # Negative probability.
            @test_throws Exception Kokako.MarkovianGraph([[-0.5 0.75]])
            # Proability sums to greater than 1.
            @test_throws Exception Kokako.MarkovianGraph([[0.8 0.8]])
            # Mis-matched dimensions.
            @test_throws Exception Kokako.MarkovianGraph([
                [0.1 0.2 0.7], [0.5 0.5; 0.5 0.5]
            ])
        end
        @testset "keyword vs list" begin
            graph_1 = Kokako.MarkovianGraph(
                stages = 2,
                transition_matrix = [0.4 0.6; 0.25 0.75],
                root_node_transition = [0.7, 0.3]
            )
            graph_2 = Kokako.MarkovianGraph([
                [0.7 0.3], [0.4 0.6; 0.25 0.75]
            ])
            @test graph_1.root_node == graph_2.root_node
            @test graph_1.nodes == graph_2.nodes
        end
    end

    @testset "Graph" begin
        @testset "Construct Graph" begin
            graph = Kokako.Graph(:root)
            @test graph.root_node == :root
            @test collect(keys(graph.nodes)) == [:root]
        end
        @testset "Add node" begin
            graph = Kokako.Graph(:root)
            Kokako.add_node(graph, :x)
            @test collect(keys(graph.nodes)) == [:root, :x]
        end
        @testset "Add duplicate node" begin
            graph = Kokako.Graph(:root)
            Kokako.add_node(graph, :x)
            @test_throws Exception Kokako.add_node(graph, :x)
        end
        @testset "Add edge" begin
            graph = Kokako.Graph(:root)
            Kokako.add_node(graph, :x)
            Kokako.add_edge(graph, :root => :x, 1.0)
            @test haskey(graph.nodes, :root)
            @test graph.nodes[:root] == [(:x, 1.0)]
        end
        @testset "Add edge of wrong type" begin
            graph = Kokako.Graph(:root)
            @test_throws Exception Kokako.add_node(graph, 1)
        end
        @testset "Add edge to missing node" begin
            graph = Kokako.Graph(:root)
            Kokako.add_node(graph, :x)
            @test_throws Exception Kokako.add_edge(graph, :x => :y, 1.0)
            @test_throws Exception Kokako.add_edge(graph, :y => :x, 1.0)
        end
        @testset "Add edge to root" begin
            graph = Kokako.Graph(:root)
            Kokako.add_node(graph, :x)
            @test_throws Exception Kokako.add_edge(graph, :x => :root, 1.0)
        end
        @testset "Invalid probability" begin
            graph = Kokako.Graph(:root)
            Kokako.add_node(graph, :x)
            Kokako.add_edge(graph, :root => :x, 0.5)
            Kokako.add_edge(graph, :root => :x, 0.75)
            @test_throws Exception Kokako.validate_graph(graph)
        end
    end
end

@testset "PolicyGraph constructor" begin
    @testset "LinearGraph" begin
        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   lower_bound = 0.0,
                                   direct_mode=false) do node, stage
        end

        @test_throws Exception Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                                  lower_bound = 0.0
                                                      ) do node, stage
        end
        nodes = Set{Int}()
        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   lower_bound = 0.0,
                                   optimizer = with_optimizer(GLPK.Optimizer)
                                       ) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([1, 2])
    end

    @testset "MarkovianGraph" begin
        graph = Kokako.MarkovianGraph([
                ones(Float64, (1, 1)),
                [0.5 0.5],
                [0.5 0.5; 0.3 0.4],
                [0.5 0.5; 0.3 0.4],
                [0.5 0.5; 0.3 0.4]
            ]
        )
        nodes = Set{Tuple{Int, Int}}()
        model = Kokako.PolicyGraph(graph, lower_bound = 0.0,
                                   direct_mode = false) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (5, 2)])
    end

    @testset "General" begin
        graph = Kokako.Graph(
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
        model = Kokako.PolicyGraph(graph, lower_bound = 0.0,
                                   direct_mode = false) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([:stage_1, :stage_2, :stage_3])
    end
end

@testset "Kokako.State" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               lower_bound = 0.0,
                               direct_mode = false) do node, stage
        @variable(node, x, Kokako.State, initial_value = 0)
    end
    for stage in 1:2
        node = model[stage]
        @test haskey(node.states, :x)
        @test length(keys(node.states)) == 1
        @test node.states[:x] == node.subproblem[:x]
    end
end

@testset "Kokako.parameterize" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               lower_bound = 0.0,
                               direct_mode = false) do node, stage
        @variable(node, 0 <= x <= 1)
        Kokako.parameterize(node, [1, 2, 3], [0.4, 0.5, 0.1]) do ω
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

@testset "Kokako.set_stage_objective" begin
    @testset ":Min" begin
        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   lower_bound = 0.0,
                                   direct_mode = false) do node, stage
            @variable(node, 0 <= x <= 1)
            @stageobjective(node, 2x)
        end
        node = model[2]
        @test node.stage_objective == 2 * node.subproblem[:x]
        @test model.objective_sense == Kokako.MOI.MIN_SENSE
    end

    @testset ":Max" begin
        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   lower_bound = 0.0,
                                   sense = :Max,
                                   direct_mode = false) do node, stage
            @variable(node, 0 <= x <= 1)
            @stageobjective(node, 2x)
        end
        node = model[2]
        @test node.stage_objective == 2 * node.subproblem[:x]
        @test model.objective_sense == Kokako.MOI.MAX_SENSE
    end
end
