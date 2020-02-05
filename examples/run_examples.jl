#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Random
using Test

const EXCLUDED_EXAMPLES = [
    "run_examples.jl",
    "inventory_management.jl",
    "msppy_hydro_thermal.jl",
    "tiger_problem.jl",
]

@testset "run_examples.jl" begin
    @testset "$(example)" for example in filter(
        ex -> endswith(ex, ".jl") && !(ex in EXCLUDED_EXAMPLES),
        readdir(@__DIR__),
    )
        Random.seed!(12345)
        include(example)
    end
end
