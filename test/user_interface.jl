#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using GLPK
using SDDP
using Test

@testset "Basic Graphs" begin
    @testset "LinearGraph" begin
        graph = SDDP.LinearGraph(5)
        @test graph.root_node == 0
        for stage = 0:4
            @test haskey(graph.nodes, stage)
            @test graph.nodes[stage] == [(stage + 1, 1.0)]
        end
        @test haskey(graph.nodes, 5)
        @test graph.nodes[5] == Tuple{Int,Float64}[]

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
        @test length(graph.belief_partition) == 0
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
            @test_throws Exception SDDP.MarkovianGraph([[0.1 0.2 0.7], [0.5 0.5; 0.5 0.5]])
        end
        @testset "keyword vs list" begin
            graph_1 = SDDP.MarkovianGraph(
                stages = 2,
                transition_matrix = [0.4 0.6; 0.25 0.75],
                root_node_transition = [0.7, 0.3],
            )
            graph_2 = SDDP.MarkovianGraph([[0.7 0.3], [0.4 0.6; 0.25 0.75]])
            @test graph_1.root_node == graph_2.root_node
            @test graph_1.nodes == graph_2.nodes
            @test length(graph_1.belief_partition) == 0
            @test length(graph_2.belief_partition) == 0
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
        @testset "Belief partition" begin
            graph = SDDP.Graph(:root)
            SDDP.add_node(graph, :x)
            SDDP.add_node(graph, :y)
            @test_throws ErrorException(
                "You must provide on Lipschitz contsant for every element in " *
                "the ambiguity set.",
            ) SDDP.add_ambiguity_set(graph, [:x], Float64[])
            @test_throws ErrorException("Cannot provide negative Lipschitz constant: [-1.0]") SDDP.add_ambiguity_set(
                graph,
                [:x],
                -1.0,
            )
            SDDP.add_ambiguity_set(graph, [:x])
            SDDP.add_ambiguity_set(graph, [:y])
            @test graph.belief_partition == [[:x], [:y]]

            graph = SDDP.Graph(
                :root,
                [:x, :y],
                [(:root => :x, 0.5), (:root => :y, 0.5)],
                belief_partition = [[:x, :y]],
                belief_lipschitz = [[1.0, 1.0]],
            )
            @test graph.belief_partition == [[:x, :y]]
            @test sprint(show, graph) == join(
                [
                    "Root",
                    " root",
                    "Nodes",
                    " x",
                    " y",
                    "Arcs",
                    " root => x w.p. 0.5",
                    " root => y w.p. 0.5",
                    "Partition",
                    " {",
                    "    x",
                    "    y",
                    " }\n",
                ],
                "\n",
            )



            graph = SDDP.Graph(:root, [:x, :y], [(:root => :x, 0.5), (:root => :y, 0.5)])
            @test length(graph.belief_partition) == 0

        end
    end
end

@testset "PolicyGraph constructor" begin
    @testset "LinearGraph" begin
        model = SDDP.PolicyGraph(
            SDDP.LinearGraph(2),
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
        end

        @test_throws Exception SDDP.PolicyGraph(
            SDDP.LinearGraph(2),
            lower_bound = 0.0,
            direct_mode = true,
        ) do node, stage
        end
        nodes = Set{Int}()
        model = SDDP.PolicyGraph(
            SDDP.LinearGraph(2),
            lower_bound = 0.0,
            optimizer = GLPK.Optimizer,
        ) do node, stage
            push!(nodes, stage)
        end
        @test nodes == Set([1, 2])
        @test sprint(show, model) == "A policy graph with 2 nodes.\n Node indices: 1, 2\n"
    end

    @testset "MarkovianGraph" begin
        graph = SDDP.MarkovianGraph([
            ones(Float64, (1, 1)),
            [0.5 0.5],
            [0.5 0.5; 0.3 0.4],
            [0.5 0.5; 0.3 0.4],
            [0.5 0.5; 0.3 0.4],
        ])
        nodes = Set{Tuple{Int,Int}}()
        model =
            SDDP.PolicyGraph(graph, lower_bound = 0.0, direct_mode = false) do node, stage
                push!(nodes, stage)
            end
        @test nodes ==
              Set([(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (5, 2)])
    end

    @testset "MarkovianPolicyGraph" begin
        nodes = Set{Tuple{Int,Int}}()
        model = SDDP.MarkovianPolicyGraph(
            transition_matrices = [
                ones(Float64, (1, 1)),
                [0.5 0.5],
                [0.5 0.5; 0.3 0.4],
                [0.5 0.5; 0.3 0.4],
                [0.5 0.5; 0.3 0.4],
            ],
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
            push!(nodes, stage)
        end
        @test nodes ==
              Set([(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (4, 1), (4, 2), (5, 1), (5, 2)])
    end

    @testset "General" begin
        graph = SDDP.Graph(
            :root,
            [:stage_1, :stage_2, :stage_3],
            [
                (:root => :stage_1, 1.0),
                (:stage_1 => :stage_2, 1.0),
                (:stage_2 => :stage_3, 1.0),
                (:stage_3 => :stage_1, 0.9),
            ],
        )
        nodes = Set{Symbol}()
        model =
            SDDP.PolicyGraph(graph, lower_bound = 0.0, direct_mode = false) do node, stage
                push!(nodes, stage)
            end
        @test nodes == Set([:stage_1, :stage_2, :stage_3])
    end
end

@testset "SDDP.State" begin
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
        @variable(node, x, SDDP.State, initial_value = 0)
    end
    for stage = 1:2
        node = model[stage]
        @test haskey(node.states, :x)
        @test length(keys(node.states)) == 1
        @test node.states[:x] == node.subproblem[:x]
    end
end

@testset "SDDP.parameterize" begin
    model = SDDP.PolicyGraph(
        SDDP.LinearGraph(2),
        lower_bound = 0.0,
        direct_mode = false,
    ) do node, stage
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
        model = SDDP.PolicyGraph(
            SDDP.LinearGraph(2),
            sense = :Min,
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
            @variable(node, 0 <= x <= 1)
            @stageobjective(node, 2x)
        end
        node = model[2]
        @test node.stage_objective == 2 * node.subproblem[:x]
        @test model.objective_sense == SDDP.MOI.MIN_SENSE

        @test_throws Exception SDDP.LinearPolicyGraph(
            stages = 2,
            sense = :Min,
            upper_bound = 0.0,
            direct_mode = false,
        ) do node, stage
        end
    end

    @testset ":Max" begin
        model = SDDP.PolicyGraph(
            SDDP.LinearGraph(2),
            upper_bound = 0.0,
            sense = :Max,
            direct_mode = false,
        ) do node, stage
            @variable(node, 0 <= x <= 1)
            @stageobjective(node, 2x)
        end
        node = model[2]
        @test node.stage_objective == 2 * node.subproblem[:x]
        @test model.objective_sense == SDDP.MOI.MAX_SENSE

        @test_throws Exception SDDP.LinearPolicyGraph(
            stages = 2,
            sense = :Max,
            lower_bound = 0.0,
            direct_mode = false,
        ) do node, stage
        end
    end
end

@testset "Errors" begin
    @testset "<=0 stages" begin
        exception =
            ErrorException("You must create a LinearPolicyGraph with `stages >= 1`.")
        @test_throws exception SDDP.LinearPolicyGraph(stages = 0) do sp, t
        end
    end
    @testset "missing bounds" begin
        exception = ErrorException(
            "You must specify a finite lower bound on the objective value" *
            " using the `lower_bound = value` keyword argument.",
        )
        @test_throws exception SDDP.LinearPolicyGraph(stages = 1, sense = :Min) do sp, t
        end

        exception = ErrorException(
            "You must specify a finite upper bound on the objective value" *
            " using the `upper_bound = value` keyword argument.",
        )
        @test_throws exception SDDP.LinearPolicyGraph(stages = 1, sense = :Max) do sp, t
        end
    end
    @testset "parameterize!" begin
        exception = ErrorException("Duplicate calls to SDDP.parameterize detected.")
        @test_throws exception SDDP.LinearPolicyGraph(
            stages = 2,
            upper_bound = 0.0,
            sense = :Max,
            direct_mode = false,
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
            "of the state variable at the root node.",
        )
        @test_throws exception SDDP.LinearPolicyGraph(
            stages = 2,
            upper_bound = 0.0,
            sense = :Max,
            direct_mode = false,
        ) do node, stage
            @variable(node, x, SDDP.State)
            @stageobjective(node, x.out)
        end
    end

    @testset "termination_status" begin
        model = SDDP.LinearPolicyGraph(
            stages = 2,
            upper_bound = 0.0,
            sense = :Max,
            direct_mode = false,
        ) do node, stage
            @variable(node, x, SDDP.State, initial_value = 0.0)
            @stageobjective(node, x.out)
        end
        @test SDDP.termination_status(model) == :model_not_solved
    end
end

@testset "Numerical stability report" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = -1e10,
        direct_mode = false,
    ) do subproblem, t
        @variable(subproblem, x >= -1e7, SDDP.State, initial_value = 1e-5)
        @variable(subproblem, 1 <= y <= 5, Int)  # Note: this is just to test range fallback
        @constraint(subproblem, 1e9 * x.out >= 1e-6 * x.in + 1e-8)
        @stageobjective(subproblem, 1e9 * x.out)
    end
    report = sprint(SDDP.numerical_stability_report, model)
    @test occursin("WARNING", report)
    report_2 = sprint(io -> SDDP.numerical_stability_report(io, model, by_node = true))
    @test occursin("Numerical stability report for node: 1", report_2)
    @test occursin("Numerical stability report for node: 2", report_2)
end


@testset "objective_state" begin
    model = SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0,
        direct_mode = false,
    ) do subproblem, t
        @variable(subproblem, x, SDDP.State, initial_value = 0)
        SDDP.parameterize(subproblem, [1, 2]) do ω
            price = SDDP.objective_state(subproblem)
            @stageobjective(subproblem, price * x.out)
        end
    end
    @test_throws ErrorException("No objective state defined.") SDDP.simulate(model, 1)

    @test_throws ErrorException("add_objective_state can only be called once.") SDDP.LinearPolicyGraph(
        stages = 2,
        lower_bound = 0,
        direct_mode = false,
    ) do subproblem, t
        @variable(subproblem, x, SDDP.State, initial_value = 0)
        SDDP.add_objective_state(
            subproblem,
            initial_value = 1.5,
            lower_bound = 0.75,
            upper_bound = 2.25,
            lipschitz = 100.0,
        ) do y, ω
            return y + ω
        end
        SDDP.add_objective_state(
            subproblem,
            initial_value = 1.5,
            lower_bound = 0.75,
            upper_bound = 2.25,
            lipschitz = 100.0,
        ) do y, ω
            return y + ω
        end
        SDDP.parameterize(subproblem, [1, 2]) do ω
            price = SDDP.objective_state(subproblem)
            @stageobjective(subproblem, price * x.out)
        end
    end
end

@testset "Belief Updater" begin
    graph = SDDP.LinearGraph(2)
    SDDP.add_edge(graph, 2 => 1, 0.9)
    model =
        SDDP.PolicyGraph(graph, lower_bound = 0.0, direct_mode = false) do subproblem, node
            beliefs = [[0.2, 0.8], [0.7, 0.3]]
            SDDP.parameterize(subproblem, [:A, :B], beliefs[node]) do ω
                return nothing
            end
        end
    belief_updater = SDDP.construct_belief_update(model, [Set([1]), Set([2])])
    belief = Dict(1 => 1.0, 2 => 0.0)
    belief′ = copy(belief)
    @test belief_updater(belief′, belief, 2, :A) == Dict(1 => 0.0, 2 => 1.0)
    @test belief′ == Dict(1 => 0.0, 2 => 1.0)
    belief = Dict(1 => 0.0, 2 => 1.0)
    @test belief_updater(belief′, belief, 1, :B) == Dict(1 => 1.0, 2 => 0.0)
    @test belief′ == Dict(1 => 1.0, 2 => 0.0)


    belief_updater = SDDP.construct_belief_update(model, [Set([1, 2])])

    belief = Dict(1 => 1.0, 2 => 0.0)
    @test belief_updater(belief′, belief, 1, :A) == Dict(1 => 0.0, 2 => 1.0)
    belief = Dict(1 => 0.0, 2 => 1.0)
    @test belief_updater(belief′, belief, 1, :B) == Dict(1 => 1.0, 2 => 0.0)

    function is_approx(x::Dict{T,Float64}, y::Dict{T,Float64}) where {T}
        if length(x) != length(y)
            return false
        end
        for (key, value) in x
            if !(value ≈ y[key])
                return false
            end
        end
        return true
    end
    belief = Dict(1 => 0.6, 2 => 0.4)
    @test is_approx(belief_updater(belief′, belief, 1, :A), Dict(1 => 6 / 41, 2 => 35 / 41))
end

@testset "Ensure root printed first" begin
    g = SDDP.Graph(:root, [:a], [(:root => :a, 1.0)])
    @test sprint(show, g) == """
    Root
     root
    Nodes
     a
    Arcs
     root => a w.p. 1.0
    """
end

@testset "Tuple{Int,Float64} nodes sorted" begin
    @test SDDP.sort_nodes([(1, 1.0), (2, 0.1), (1, 0.5)]) == [(1, 0.5), (1, 1.0), (2, 0.1)]
    g = SDDP.Graph((0, 0.0), [(1, 1.0), (2, 0.1), (1, 0.5)], [((0, 0.0) => (2, 0.1), 1.0)])
    @test sprint(show, g) == """
    Root
     (0, 0.0)
    Nodes
     (1, 0.5)
     (1, 1.0)
     (2, 0.1)
    Arcs
     (0, 0.0) => (2, 0.1) w.p. 1.0
    """
end

@testset "String nodes unsorted" begin
    @test SDDP.sort_nodes(["c", "b"]) == ["c", "b"]
    g = SDDP.Graph("a", ["c", "b"], [("a" => "b", 1.0), ("b" => "c", 1.0)])
    @test sprint(show, g) == """
    Root
     a
    Nodes
     c
     b
    Arcs
     a => b w.p. 1.0
     b => c w.p. 1.0
    """
end
