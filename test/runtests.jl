using Kokako, Test

for file in ["graphs.jl", "policy_graphs.jl", "states.jl"]
    @testset "$(file)" begin
        include(file)
    end
end
