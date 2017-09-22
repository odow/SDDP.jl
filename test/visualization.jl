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
        @testset "Two prices" begin
            cuts = [
                SDDP.Cut(1.0, [-1.0, -1.0]),
                SDDP.Cut(0.5, [-0.25, 0.0])
            ]
            x, y = SDDP.processvaluefunctiondata(cuts, true, [0.0, 1.0], [0.0, 1.0])
            @test x == [0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0]
            @test y == [1.0, 0.5, 0.25, 0.25]

            (plotly_data, scene_text) = SDDP.getplotlydata(x, y, "lab1", "lab2")
            @test plotly_data["x"] == x[1,:]
            @test plotly_data["y"] == x[2,:]
            @test plotly_data["z"] == y
            @test plotly_data["type"] == "scatter3d"
        end
        @testset "One price" begin
            cuts = [
                SDDP.Cut(1.0, [-1.0, -1.0]),
                SDDP.Cut(0.5, [-0.25, 0.0])
            ]
            x, y = SDDP.processvaluefunctiondata(cuts, false, [0.0, 1.0], 0.5)
            @test x == [0.0 1.0;]
            @test y == [0.5, -0.5]

            (plotly_data, scene_text) = SDDP.getplotlydata(x, y, "lab1", "lab2")
            @test plotly_data["x"] == [0.0, 1.0]
            @test plotly_data["y"] == y
            @test haskey(plotly_data, "z") == false
            @test plotly_data["type"] == "scatter"
            @test scene_text == "xaxis: {title: \"lab1\"}, yaxis: {title: \"Future Cost\"}"
        end
        @testset "fallback" begin
            cuts = [
                SDDP.Cut(1.0, [-1.0, -1.0]),
                SDDP.Cut(0.5, [-0.25, 0.0])
            ]
            vf = SDDP.DefaultValueFunction()
            push!(vf.cutmanager.cuts, cuts[1])
            push!(vf.cutmanager.cuts, cuts[2])

            x, y = SDDP.processvaluefunctiondata(vf, false, [0.0, 1.0], 0.5)
            @test x == [0.0 1.0;]
            @test y == [0.5, -0.5]
        end

        @testset "Simulation Example" begin
            include(joinpath(examples_dir, "HydroValleys", "hydro_valley.jl"))
            # stagewise inflows and markov prices
            srand(12345)
            
            markov_stagewise_model = hydrovalleymodel(hasstagewiseinflows=true, hasmarkovprice=true)
            SDDP.solve(markov_stagewise_model, max_iterations=10, print_level=0)
            @test isapprox(getbound(markov_stagewise_model), 855.0, atol=1e-3)

            results = simulate(markov_stagewise_model, 5, [:reservoir])
            @test length(results) == 5

            p = SDDP.newplot()
            SDDP.addplot!(p, 1:5, 1:3, (i, t)->round(results[i][:stageobjective][t], 2),
                title="Accumulated Profit", ylabel="Accumulated Profit (\$)", cumulative=true)
            SDDP.addplot!(p, 1:5, 1:3, (i, t)->round(results[i][:stageobjective][t], 2),
                title="Weekly Income", ylabel="Week Profit (\$)")
            SDDP.addplot!(p, 1:5, 1:3, (i, t)->round(results[i][:reservoir][t][1], 2),
                title="Upper Reservoir", ylabel="Level")
            SDDP.addplot!(p, 1:5, 1:3, (i, t)->round(results[i][:reservoir][t][2], 2),
                title="Lower Reservoir")
            prices = [1 2 0; 2 1 0; 3 4 0]
            SDDP.addplot!(p, 1:5, 1:3, (i, t)->round(prices[t, results[i][:markov][t]], 2),
                ylabel="Price", interpolate="step-after")
            # for real-world display, use
            # SDDP.show(p)

            html = SDDP.prephtml(p)
            @test html == readstring(joinpath(dirname(@__FILE__), "visualize_examples", "simulation.html"))
        end

        @testset "Value Function Example" begin
            include(joinpath(examples_dir, "visualize_value_function.jl"))
            html = SDDP.prepvaluefunctionplot(m, 3, 1, "Bonds", "", 50.0, 0.0:2.0:100)
            @test html == readstring(joinpath(dirname(@__FILE__), "visualize_examples", "value_function.html"))
        end


    end
end
