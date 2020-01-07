#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

@testset "find_min" begin
    @test SDDP.find_min([1.0, 2.0, 3.0], 2.1) == (abs(2.0 - 2.1), 2)
    @test SDDP.find_min([1.0, 2.0, 3.0], 0.0) == (1.0, 1)
    @test SDDP.find_min([1.0, 2.0, 3.0], 5.0) == (2.0, 3)
end

@testset "allocate_support_budget" begin
    f() = rand(10)
    @inferred SDDP.allocate_support_budget(f, 20, 100)
    states = SDDP.allocate_support_budget(f, 20, 100)
    @test isa(states, Vector{Int})
    @test sum(states) == 20
    @test all(states .> 0)

    f() = [1, 2, 3 + rand()]
    @inferred SDDP.allocate_support_budget(f, 17, 31)
    states = SDDP.allocate_support_budget(f, 17, 31)
    @test sum(states) == 17
    @test all(states .> 0)

    f() = [1.0, 2.0, 3.0]
    @inferred SDDP.allocate_support_budget(f, 5, 10)
    states = SDDP.allocate_support_budget(f, 5, 10)
    @test states == [1, 1, 1]
    @test all(states .> 0)

    states = [1, 3, 5]
    @test states === SDDP.allocate_support_budget(f, states, 19)

    @test SDDP.allocate_support_budget(() -> rand(3), 2, 10) == [1, 1, 1]
end

@testset "lattice_approximation" begin
    support, probability = SDDP.lattice_approximation(() -> rand(5), [1, 2, 3, 4, 5], 100)
    for (t, s) in enumerate(support)
        @test length(s) == t
        @test all(x -> 0 <= x <= 1, s)
        @test !any(isnan, probability[t])
        @test all(isapprox.(sum(probability[t], dims = 2), 1.0))
    end
end

@testset "MarkovianGraph" begin
    g = SDDP.MarkovianGraph(() -> rand(5), budget = 10, scenarios = 100)
    @test g.root_node == (0, 0.0)
    @test length(g.nodes) == 11
    for (k, node) in g.nodes
        if length(node) > 0
            @test sum(arc[2] for arc in node) ≈ 1.0
        else
            @test length(node) == 0
        end
    end
end
