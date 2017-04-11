struct MyRiskMeasure <: SDDP.AbstractRiskMeasure end

@testset "Risk Measures" begin
    @testset "Expectation" begin
        ex = Expectation()
        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        SDDP.modifyprobability!(ex, y, x, obj)
        @test y == x
    end

    @testset "User-defined MyRiskMeasure" begin
        measure = MyRiskMeasure()
        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        @test_throws Exception SDDP.modifyprobability!(measure, y, x, obj)
    end
end
