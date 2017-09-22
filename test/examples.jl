#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

@testset "Examples" begin
    @testset "HydroValleys" begin
        for example in [
                "hydro_valley_tests.jl",
                "simplified_hydrothermal_dispatch.jl"
            ]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, "HydroValleys", example))
            end
        end
    end

    @testset "Binary Problems" begin
        for example in [
                "booking_management.jl",
                "vehicle_location.jl"
            ]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, "BinaryProblems", example))
            end
        end
    end

    @testset "SDDP.jl" begin
        for example in [
                "newsvendor.jl",
                "newsvendor_historic_simulation.jl",
                "async_newsvendor.jl",
                "load_cuts.jl",
                "asset_management_stagewise.jl",
                "asset_management.jl",
                "simple_objective_noise.jl"
            ]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, example))
            end
        end
    end

    @testset "StochDynamicProgramming.jl" begin
        for example in [
                "stock-example.jl",
                "multistock-example.jl"
            ]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, "StochDynamicProgramming.jl", example))
            end
        end
    end

    @testset "StructDualDynProg.jl" begin
        for example in [
                "prob5.2_2stages.jl",
                "prob5.2_3stages.jl"
            ]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, "StructDualDynProg.jl", example))
            end
        end
    end

    @testset "FAST" begin
        for example in [
                "hydro_thermal.jl",
                "quickstart.jl",
                "production_management_multiple_stages.jl"
            ]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, "FAST", example))
            end
        end
    end
end
