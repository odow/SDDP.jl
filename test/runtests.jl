#  Copyright 2017-21, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Distributed
using Random
using SDDP
using Test

function read_dir(dir, exclude = String[])
    return filter(s -> !(s in exclude) && endswith(s, ".jl"), readdir(dir))
end

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "docs", "src", "examples")
const PLUGINS_DIR = joinpath(@__DIR__, "plugins")
const VISUALIZATION_DIR = joinpath(@__DIR__, "visualization")

@testset "SDDP.jl" begin
    @testset "Unit Tests" begin
        @testset "plugins/$(file)" for file in read_dir(
            PLUGINS_DIR,
            ["parallel_schemes.jl"],
        )
            @info file
            include(joinpath(PLUGINS_DIR, file))
        end
        @testset "visualization/$(file)" for file in read_dir(VISUALIZATION_DIR)
            @info file
            include(joinpath(VISUALIZATION_DIR, file))
        end
        @testset "$(file)" for file in read_dir(".", ["runtests.jl"])
            @info file
            include(file)
        end
    end

    @testset "Examples" begin
        @testset "$example" for example in read_dir(EXAMPLES_DIR)
            @info example
            Random.seed!(12345)
            include(joinpath(EXAMPLES_DIR, example))
        end
    end

    @testset "Parallel" begin
        procs = Distributed.addprocs(4)
        include(joinpath(@__DIR__, "plugins", "parallel_schemes.jl"))
        Distributed.rmprocs(procs)
    end
end
