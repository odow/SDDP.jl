using Kokako, Test, Random

@testset "Unit Tests" begin
    for file in ["risk_measures.jl",
                 "sampling_schemes.jl",
                 "sddp.jl",
                 "user_interface.jl"]
        @testset "$(file)" begin
            include(file)
        end
    end
end

const EXAMPLES = Dict(
    "agriculture" => [
        "mccardle_farm.jl"],
    "FAST" => [
        "hydro_thermal.jl",
        "production_management.jl",
        "quickstart.jl"],
    "infinite_horizon" => [
        "trivial_example.jl",
        "hydro_thermal.jl"],
    "StructDualDynProg.jl" => [
        "prob5.2_2stages.jl",
        "prob5.2_3stages.jl"],
    "StochDynamicProgramming.jl" => [
        "stock_example.jl",
        "multistock_example.jl"]
)

const EXAMPLES_DIR = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "Examples" begin
    for (directory, example_list) in EXAMPLES
        @testset "$(directory)" begin
            for example in example_list
                @testset "$(example)" begin
                    Random.seed!(12345)
                    evalfile(joinpath(EXAMPLES_DIR, directory, example))
                end
            end
        end
    end
end
