using Kokako, Test, Random

@testset "Unit Tests" begin
    for file in ["plugins/risk_measures.jl",
                 "plugins/sampling_schemes.jl",
                 "plugins/stopping_rules.jl",
                 "sddp.jl",
                 "user_interface.jl"]
        @testset "$(file)" begin
            include(file)
        end
    end
end

const EXAMPLES = [
    "agriculture_mccardle_farm.jl",
    "asset_management_simple.jl",
    "FAST_hydro_thermal.jl",
    "FAST_production_management.jl",
    "FAST_quickstart.jl",
    "infinite_horizon_trivial.jl",
    "infinite_horizon_hydro_thermal.jl",
    "StructDualDynProg.jl_prob5.2_2stages.jl",
    "StructDualDynProg.jl_prob5.2_3stages.jl",
    "StochDynamicProgramming.jl_stock.jl",
    "StochDynamicProgramming.jl_multistock.jl"
]

const EXAMPLES_DIR = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "Examples" begin
    for example in EXAMPLES
        @testset "$(example)" begin
            Random.seed!(12345)
            include(joinpath(EXAMPLES_DIR, example))
        end
    end
end
