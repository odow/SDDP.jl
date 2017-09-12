#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

@testset "Visualization" begin
    @testset "Simulation" begin
        p = SDDP.newplot()
        SDDP.addplot!(p, 1:2, 1:2, (i,t)->i+t;
            xlabel="x", ylabel="y", cumulative=true, title="Title",
            interpolate="step-before", ymin="0", ymax="1")
        @test length(p.data) == 1
        d = p.data[1]
        @test d["xlabel"] == "x"
        @test d["ylabel"] == "y"
        @test d["title"] == "Title"
        @test d["cumulative"] == true
        @test d["interpolate"] == "step-before"
        @test d["ymin"] == "0"
        @test d["ymax"] == "1"
        @test d["data"] == [[2.0,5.0], [3.0,7.0]]

        SDDP.addplot!(p, 1:3, 2:3, (i,t)->i-t;
            xlabel="x2", ylabel="y2", cumulative=false, title="Title2",
            interpolate="step-after", ymin="20", ymax="12")
        @test length(p.data) == 2
        d = p.data[2]
        @test d["xlabel"] == "x2"
        @test d["ylabel"] == "y2"
        @test d["title"] == "Title2"
        @test d["cumulative"] == false
        @test d["interpolate"] == "step-after"
        @test d["ymin"] == "20"
        @test d["ymax"] == "12"
        @test d["data"] == [[-1.0, -2.0], [0.0, -1.0], [1.0, 0.0]]
    end

    @testset "ValueFunction" begin
        cuts = [
            SDDP.Cut(1.0, [-1.0, -1.0]),
            SDDP.Cut(0.5, [-0.25, 0.0])
        ]
        x, y = SDDP.processvaluefunctiondata(cuts, true, [0.0, 1.0], [0.0, 1.0])
        @test x == [0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0]
        @test y == [1.0, 0.5, 0.25, 0.25]

        x, y = SDDP.processvaluefunctiondata(cuts, false, [0.0, 1.0], 0.5)
        @test x == [0.0 1.0;]
        @test y == [0.5, -0.5]
    end
end
