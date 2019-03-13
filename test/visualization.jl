#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Kokako, Test

@testset "SpaghettiPlot" begin
    simulations = [
        [
            Dict(:x => 1.0, :y => 4.0),
            Dict(:x => 2.0, :y => 5.0),
            Dict(:x => 3.0, :y => 6.0)
        ], [
            Dict(:x => 1.5, :y => 4.5),
            Dict(:x => 2.5, :y => 5.5),
            Dict(:x => 3.5, :y => 6.5)
        ]
    ]
    plt = Kokako.SpaghettiPlot(simulations)
    Kokako.add_spaghetti(plt, cumulative = true) do data
        return data[:x]
    end
    Kokako.add_spaghetti(plt, title = "y") do data
        return 2 * data[:y]
    end
    Kokako.prep_html(plt, "test.html")
    @test read("test.html", String) == read("control.html", String)
    rm("test.html")
end
