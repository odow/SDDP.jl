#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const examples_dir = joinpath(dirname(dirname(@__FILE__)), "examples")

const Examples = Dict(
    "HydroValleys" => [
            "hydro_valley_tests.jl",
            "simplified_hydrothermal_dispatch.jl"
        ],
    "BinaryProblems" => [
            "booking_management.jl",
            "vehicle_location.jl",
            "airconditioning.jl"
        ],
    "Newsvendor" => [
            "async_newsvendor.jl",
            "newsvendor.jl"
        ],
    "StochDynamicProgramming.jl" => [
            "stock-example.jl",
            "multistock-example.jl"
        ],
    "StructDualDynProg.jl" => [
            "prob5.2_2stages.jl",
            "prob5.2_3stages.jl"
        ],
    "FAST" => [
            "hydro_thermal.jl",
            "quickstart.jl",
            "production_management_multiple_stages.jl"
        ],
    "AssetManagement" => [
            "asset_management_stagewise.jl",
            "asset_management.jl"
        ],
    "Agriculture" => [
            "mccardle_farm.jl"
        ],
    "Nonlinear" => [
            "logs.jl",
            "inventory-control.jl",
            "log_auto_regressive.jl"
        ],
    "PriceInterpolation" => [
            "simple_contracting.jl",
            "newsvendor.jl",
            "widget_producer.jl",
            "riverchain.jl",
            "simple_markov.jl"
        ]

)

@testset "Examples" begin
    for (key, examples) in Examples
        @testset "$(key)" begin
            for example in examples
                @testset "$example" begin
                    println("Running $(key):$(example)")
                    include(joinpath(examples_dir, key, example))
                end
            end
        end
    end

    @testset "SDDP.jl" begin
        for example in ["simple_objective_noise.jl"]
            @testset "$example" begin
                println("Running $(example)")
                include(joinpath(examples_dir, example))
            end
        end
    end
end
