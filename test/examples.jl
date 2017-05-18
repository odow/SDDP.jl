#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

const examples_dir = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "Examples" begin
    for example in [
            "hydro_valley_deterministic.jl",
            "hydro_valley_stagewise.jl",
            "hydro_valley_markov.jl",
            "hydro_valley_stagewise_markov.jl",
            "newsvendor.jl",
            "newsvendor_historic_simulation.jl",
            "stock-example.jl",
            "multistock-example.jl",
            "risk_aversion.jl",
            "dematos_cutselection.jl",
            "async_newsvendor.jl"
        ]
        @testset "$example" begin
            include(joinpath(examples_dir, example))
        end
    end
end
