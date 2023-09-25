#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module TestModelingAids

using SDDP
using Test

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

function test_find_min()
    @test SDDP.find_min([1.0, 2.0, 3.0], 2.1) == (abs(2.0 - 2.1), 2)
    @test SDDP.find_min([1.0, 2.0, 3.0], 0.0) == (1.0, 1)
    @test SDDP.find_min([1.0, 2.0, 3.0], 5.0) == (2.0, 3)
    return
end

function test__allocate_support_budget()
    @inferred SDDP._allocate_support_budget(() -> rand(10), 20, 100)
    states = SDDP._allocate_support_budget(() -> rand(10), 20, 100)
    @test isa(states, Vector{Int})
    @test sum(states) == 20
    @test all(states .> 0)
    @inferred SDDP._allocate_support_budget(() -> [1, 2, 3 + rand()], 17, 31)
    states = SDDP._allocate_support_budget(() -> [1, 2, 3 + rand()], 17, 31)
    @test sum(states) == 17
    @test all(states .> 0)
    @inferred SDDP._allocate_support_budget(() -> [1.0, 2.0, 3.0], 5, 10)
    states = SDDP._allocate_support_budget(() -> [1.0, 2.0, 3.0], 5, 10)
    @test states == [1, 1, 1]
    @test all(states .> 0)
    states = [1, 3, 5]
    new_states =
        SDDP._allocate_support_budget(() -> [1.0, 2.0, 3.0], states, 19)
    @test states == new_states
    @test SDDP._allocate_support_budget(() -> rand(3), 2, 10) == [1, 1, 1]
    return
end

function test__lattice_approximation()
    support, probability =
        SDDP._lattice_approximation(() -> rand(5), [1, 2, 3, 4, 5], 100)
    for (t, s) in enumerate(support)
        @test length(s) == t
        @test all(x -> 0 <= x <= 1, s)
        @test !any(isnan, probability[t])
        @test all(isapprox.(sum(probability[t], dims = 2), 1.0))
    end
    return
end

function test_MarkovianGraph()
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
    return
end

function test_duplicate_nodes()
    function simulator()
        inflow = zeros(3)
        current = 50.0
        Ω = [-10.0, 0.1, 9.6]
        for t in 1:3
            current += rand(Ω)
            inflow[t] = current
        end
        return inflow
    end
    num_nodes = Int[]
    for _ in 1:100
        g = SDDP.MarkovianGraph(simulator; budget = 8, scenarios = 30)
        push!(num_nodes, length(g.nodes))
    end
    @test minimum(num_nodes) < 9
    @test maximum(num_nodes) == 9
    return
end

end  # module

TestModelingAids.runtests()
