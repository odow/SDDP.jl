#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Kokako, Test, Random

@testset "Unit Tests" begin
    for file in ["plugins/risk_measures.jl",
                 "plugins/sampling_schemes.jl",
                 "plugins/stopping_rules.jl",
                 "sddp.jl",
                 "user_interface.jl",
                 "visualization.jl"]
        @testset "$(file)" begin
            include(file)
        end
    end
end

const EXCLUDED_EXAMPLES = [
    "inventory_management.jl",
    "daniel_hydro_complete.jl"
]

const EXAMPLES_DIR = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "Examples" begin
    @testset "$(example)" for example in filter(
            s-> !(s in EXCLUDED_EXAMPLES) && endswith(s, ".jl"), readdir(EXAMPLES_DIR))
        Random.seed!(12345)
        include(joinpath(EXAMPLES_DIR, example))
    end
end
