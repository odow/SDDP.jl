#  Copyright (c) 2017-24, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
    if Sys.WORD_SIZE == 64
        # This fails on 32-bit machines.
        @test read("test.html", String) == read(control, String)
    end
    SDDP.save(plt, "test.html"; open = false)
    @test sprint(show, plt) == "A spaghetti plot with 2 scenarios and 3 stages."
    if Sys.WORD_SIZE == 64
        # This fails on 32-bit machines.
        @test read("test.html", String) == read(control, String)
    end
    rm("test.html")
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

end  # module

TestVisualization.runtests()
