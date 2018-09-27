using Kokako, Test

@testset "Linear Graphs" begin
    graph = Kokako.LinearGraph(5)
    @test graph.root_node == 0
    for stage in 0:4
        @test haskey(graph.nodes, stage)
        @test graph.nodes[stage] == [(stage+1, 1.0)]
    end
    @test haskey(graph.nodes, 5)
    @test graph.nodes[5] == Tuple{Int, Float64}[]
end

@testset "Markovian Graph" begin
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

@testset "General Graph" begin
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
