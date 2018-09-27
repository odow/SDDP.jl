using Kokako, Test

for file in ["user_interface.jl"]
    @testset "$(file)" begin
        include(file)
    end
end
