using Kokako, Test

for file in ["graphs.jl", "policy_graphs.jl"]
    @testset "$(file)" begin
        include(file)
    end
end
