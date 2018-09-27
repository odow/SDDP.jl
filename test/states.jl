using Kokako, Test

@testset "Linear" begin
    model = Kokako.PolicyGraph(Kokako.LinearGraph(2),
                               direct_mode=false) do subproblem, stage
        @variables(subproblem, begin
            x
            x′
        end)
        Kokako.add_state_variable(subproblem, :x, x, x′)
    end
    for stage in 1:2
        node = model[stage]
        ext = Kokako.extension(node)
        @test haskey(ext.states, :x)
        @test length(keys(ext.states)) == 1
        @test ext.states[:x] == Kokako.State(node[:x], node[:x′])
    end
end
