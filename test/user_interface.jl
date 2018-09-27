using Kokako, Test, GLPK

@testset "Basic Graphs" begin
    @testset "LinearGraph" begin
        graph = Kokako.LinearGraph(5)
        @test graph.root_node == 0
        for stage in 0:4
            @test haskey(graph.nodes, stage)
            @test graph.nodes[stage] == [(stage+1, 1.0)]
        end
        @test haskey(graph.nodes, 5)
        @test graph.nodes[5] == Tuple{Int, Float64}[]
    end

    @testset "MarkovianGraph" begin
        @testset "Error cases" begin
            # Not root transition matrix.
            @test_throws Exception Kokako.MarkovianGraph([ [0.5 0.5; 0.5 0.5] ])
            # Negative probability.
            @test_throws Exception Kokako.MarkovianGraph([ [-0.5 0.75] ])
            # Proability sums to greater than 1.
            @test_throws Exception Kokako.MarkovianGraph([ [0.8 0.8] ])
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
                                   direct_mode=false) do node, stage
        end

        @test_throws Exception Kokako.PolicyGraph(Kokako.LinearGraph(2)
                                                      ) do node, stage
        end

        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   optimizer = with_optimizer(GLPK.Optimizer)
                                       ) do node, stage
        end
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
        model = Kokako.PolicyGraph(graph, direct_mode=false) do node, stage
        end
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
        model = Kokako.PolicyGraph(graph, direct_mode=false) do node, stage
        end
    end
end


@testset "Kokako.add_state_variable" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode=false) do node, stage
        @variable(node, x)
        @variable(node, x′)
        Kokako.add_state_variable(node, :x, x, x′)
    end
    for stage in 1:2
        node = model[stage]
        ext = Kokako.extension(node)
        @test haskey(ext.states, :x)
        @test length(keys(ext.states)) == 1
        @test ext.states[:x] == Kokako.State(node[:x], node[:x′])
    end
end

@testset "Kokako.parameterize" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode=false) do node, stage
        @variable(node, 0 <= x <= 1)
        Kokako.parameterize(node, [1, 2, 3], [0.4, 0.5, 0.1]) do ω
            JuMP.set_upper_bound(x, ω)
        end
    end
    node = model[2]
    ext = Kokako.extension(node)
    @test length(ext.noise_terms) == 3
    @test JuMP.upper_bound(node[:x]) == 1
    ext.parameterize(ext.noise_terms[2].term)
    @test JuMP.upper_bound(node[:x]) == 2
    ext.parameterize(3)
    @test JuMP.upper_bound(node[:x]) == 3
end

@testset "Kokako.set_stage_objective" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode=false) do node, stage
        @variable(node, 0 <= x <= 1)
        Kokako.set_stage_objective(node, :Min, 2x)
    end
    node = model[2]
    @test JuMP.objective_function(node, typeof(2 * node[:x])) ==
        2 * node[:x] + Kokako.extension(node).cost_to_go
end
