using Kokako, Test

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
    "Agriculture" => ["mccardle_farm.jl"],
    "StructDualDynProg.jl" => ["prob5.2_2stages.jl"]
)

const EXAMPLES_DIR = joinpath(dirname(dirname(@__FILE__)), "examples")

@testset "Examples" begin
    for (directory, example_list) in EXAMPLES
        @testset "$(directory)" begin
            for example in example_list
                @testset "$(example)" begin
                    evalfile(joinpath(EXAMPLES_DIR, directory, example))
                end
            end
        end
    end
end
