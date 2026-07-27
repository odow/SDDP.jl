# Copyright (c) 2017-26: Oscar Dowson and SDDP.jl contributors.
#
# This Source Code Form is subject to the terms of the Mozilla Public License,
# v2.0. If a copy of the MPL was not distributed with this file, You can obtain
# one at http://mozilla.org/MPL/2.0/.

module TestVisualization

using SDDP
using Test

import Plots

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

function test_SpaghettiPlot()
    simulations = [
        [
            Dict(:x => 1.0, :y => 4.0),
            Dict(:x => 2.0, :y => 5.0),
            Dict(:x => 3.0, :y => 6.0),
        ],
        [
            Dict(:x => 1.5, :y => 4.5),
            Dict(:x => 2.5, :y => 5.5),
            Dict(:x => 3.5, :y => 6.5),
        ],
    ]
    plt = SDDP.SpaghettiPlot(simulations)
    SDDP.add_spaghetti(plt; cumulative = true) do data
        return data[:x]
    end
    SDDP.add_spaghetti(plt; title = "y") do data
        return 2 * data[:y]
    end
    SDDP.plot(plt, "test.html"; open = false)
    @test sprint(show, plt) == "A spaghetti plot with 2 scenarios and 3 stages."
    control = joinpath(@__DIR__, "control.html")
    @test read("test.html", String) == read(control, String)
    SDDP.save(plt, "test.html"; open = false)
    @test sprint(show, plt) == "A spaghetti plot with 2 scenarios and 3 stages."
    @test read("test.html", String) == read(control, String)
    rm("test.html")
    @test SDDP.launch_file("test.html", identity) isa Cmd
    return
end

function test_PublicationPlot()
    simulations = [
        [Dict{Symbol,Any}(:x => 1), Dict{Symbol,Any}(:x => 5)],
        [Dict{Symbol,Any}(:x => 2), Dict{Symbol,Any}(:x => 6)],
        [Dict{Symbol,Any}(:x => 3), Dict{Symbol,Any}(:x => 4)],
    ]
    plot = SDDP.publication_plot(simulations) do data
        return data[:x]
    end
    @test plot isa Plots.Plot
    data = SDDP.publication_data(simulations, [0.0, 0.25, 0.5, 1.0], d -> d[:x])
    @test data == [1 4; 1.5 4.5; 2 5; 3 6]
    for val in (-Inf, Inf, NaN)
        simulations[2][2] = Dict{Symbol,Any}(:x => val)
        @test_throws(
            ErrorException(
                "Unable to plot `publication_plot` because stage 2 of " *
                "replication 2 contains data that is not finite. The data " *
                "function must return a finite real-valued scalar. Got: $val",
            ),
            SDDP.publication_data(simulations, [0.5], d -> d[:x]),
        )
    end
    return
end

function test_PublicationPlot_different_lengths()
    simulations = [
        [Dict{Symbol,Any}(:x => 1), Dict{Symbol,Any}(:x => 5)],
        [Dict{Symbol,Any}(:x => 2)],
        [Dict{Symbol,Any}(:x => 3), Dict{Symbol,Any}(:x => 4)],
    ]
    data = SDDP.publication_data(simulations, [0.0, 0.25, 0.5, 1.0], d -> d[:x])
    @test data == [1.0 4.0; 1.5 4.25; 2.0 4.5; 3.0 5.0]
    return
end

function test_plot_graph()
    model = SDDP.LinearPolicyGraph(; stages = 3, lower_bound = 0.0) do sp, t
        @variable(sp, 0 <= x <= 100, SDDP.State, initial_value = 0)
        @variable(sp, 0 <= u_production <= 200)
        @variable(sp, u_overtime >= 0)
        @constraint(sp, demand, x.in - x.out + u_production + u_overtime == 0)
        Ω = [[100.0], [100.0, 300.0], [100.0, 300.0]]
        SDDP.parameterize(ω -> set_normalized_rhs(demand, ω), sp, Ω[t])
        @stageobjective(sp, 100 * u_production + 300 * u_overtime + 50 * x.out)
        return
    end
    dir = mktempdir()
    filename = joinpath(dir, "plot_graph.html")
    SDDP.plot(model, filename; open = false)
    contents = read(filename, String)
    for line in [
        "{data: {id: '0', shape: 'ellipse', meta: 'x = 0.0'}},",
        "{data: {id: '1', meta: 'Node: 1\\nNoise terms: 1'}},",
        "{data: {id: '2', meta: 'Node: 2\\nNoise terms: 2', has_noise: true}},",
        "{data: {id: '3', meta: 'Node: 3\\nNoise terms: 2', has_noise: true}},",
        "{data: {id: 'edge_1', source: '0', target: '1', meta: 'From: 0\\nTo: 1\\nProbablity: 1.0'}},",
        "{data: {id: 'edge_2', source: '1', target: '2', meta: 'From: 1\\nTo: 2\\nProbablity: 1.0'}},",
        "{data: {id: 'edge_3', source: '2', target: '3', meta: 'From: 2\\nTo: 3\\nProbablity: 1.0'}},",
    ]
        @test occursin(line, contents)
    end
    return
end

end  # module

TestVisualization.runtests()
