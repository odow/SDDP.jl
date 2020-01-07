#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP
using Test

@testset "SpaghettiPlot" begin
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
    SDDP.add_spaghetti(plt, cumulative = true) do data
        return data[:x]
    end
    SDDP.add_spaghetti(plt, title = "y") do data
        return 2 * data[:y]
    end
    SDDP.plot(plt, "test.html", open = false)
    @test sprint(show, plt) == "A spaghetti plot with 2 scenarios and 3 stages."
    @test read("test.html", String) == read("control.html", String)
    SDDP.save(plt, "test.html", open = false)
    @test sprint(show, plt) == "A spaghetti plot with 2 scenarios and 3 stages."
    @test read("test.html", String) == read("control.html", String)
    rm("test.html")
end

@testset "PublicationPlot" begin
    data = SDDP.publication_data(
        [
            [Dict{Symbol,Any}(:x => 1), Dict{Symbol,Any}(:x => 5)],
            [Dict{Symbol,Any}(:x => 2), Dict{Symbol,Any}(:x => 6)],
            [Dict{Symbol,Any}(:x => 3), Dict{Symbol,Any}(:x => 4)],
        ],
        [0.0, 1.0],
        (d) -> d[:x],
    )
    @test data == [1 4; 3 6]
end
