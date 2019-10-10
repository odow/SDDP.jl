#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Random
using SDDP
using Test

function read_dir(dir, exclude = String[])
    return filter(s->!(s in exclude) && endswith(s, ".jl"), readdir(dir))
end

const EXCLUDED_EXAMPLES = [
    "inventory_management.jl",
    "msppy_hydro_thermal.jl",
    "tiger_problem.jl",
]

const EXAMPLES_DIR = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "SDDP.jl" begin
    @testset "Unit Tests" begin
        @testset "plugins/$(file)" for file in read_dir("plugins")
            include(joinpath("plugins", file))
        end
        @testset "$(file)" for file in read_dir(".", ["runtests.jl"])
            include(file)
        end
    end

    @testset "Examples" begin
        @testset "$example" for example in read_dir(EXAMPLES_DIR, EXCLUDED_EXAMPLES)
            Random.seed!(12345)
            include(joinpath(EXAMPLES_DIR, example))
        end
    end
end
