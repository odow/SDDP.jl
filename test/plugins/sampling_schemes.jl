#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Kokako, Test

@testset "InSampleMonteCarlo" begin
    @testset "Acyclic" begin
        model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                                   bellman_function = Kokako.AverageCut(),
                                   direct_mode=false) do node, stage
            @variable(node, 0 <= x <= 1)
            Kokako.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        scenario, terminated_due_to_cycle = Kokako.sample_scenario(
            model, Kokako.InSampleMonteCarlo()
        )
        @test length(scenario) == 2
        @test !terminated_due_to_cycle
        for (stage, (node, noise)) in enumerate(scenario)
            @test stage == node
            @test noise in stage * [1, 3]
        end
    end

    @testset "Cyclic" begin
        graph = Kokako.LinearGraph(2)
        Kokako.add_edge(graph, 2=>1, 0.9)
        model = Kokako.PolicyGraph(
                graph, bellman_function = Kokako.AverageCut(),
                direct_mode=false) do node, stage
            @variable(node, 0 <= x <= 1)
            Kokako.parameterize(node, stage * [1, 3], [0.5, 0.5]) do ω
                JuMP.set_upper_bound(x, ω)
            end
        end
        scenario, terminated_due_to_cycle = Kokako.sample_scenario(
            model, Kokako.InSampleMonteCarlo(
                terminate_on_dummy_leaf = false,
                max_depth = 4
            )
        )
        @test length(scenario) == 4
        @test !terminated_due_to_cycle  # Terminated due to max depth.
        for (index, (node, noise)) in enumerate(scenario)
            stage = (index - 1) % 2 + 1
            @test stage == node
            @test noise in stage * [1, 3]
        end
    end
end
