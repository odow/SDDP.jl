#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Random
using Test

function util_test_directory(dir, exclude = String[])
    for (root, _, files) in walkdir(dir)
        for file in files
            if endswith(file, ".jl") && !(file in exclude)
                @testset "$(file)" begin
                    @info file
                    Random.seed!(12345)
                    include(joinpath(root, file))
                end
            end
        end
    end
    return
end

@testset "SDDP.jl" begin
    util_test_directory(".", ["runtests.jl"])
    util_test_directory(joinpath(dirname(@__DIR__), "docs", "src", "examples"))
end
