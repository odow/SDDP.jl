using Kokako, Test

@testset "Expectation" begin
    risk_adjusted_probability = Vector{Float64}(undef, 5)
    Kokako.adjust_probability(
        Kokako.Expectation(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.3, 0.4, 0.5],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        true)
    @test risk_adjusted_probability == [0.1, 0.2, 0.3, 0.4, 0.5]

    risk_adjusted_probability = Vector{Float64}(undef, 5)
    Kokako.adjust_probability(
        Kokako.Expectation(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.3, 0.3, 0.1],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        false)
    @test risk_adjusted_probability == [0.1, 0.2, 0.3, 0.3, 0.1]
end

@testset "WorstCase" begin
    risk_adjusted_probability = Vector{Float64}(undef, 5)
    Kokako.adjust_probability(
        Kokako.WorstCase(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.0, 0.4, 0.5],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        true)
    @test risk_adjusted_probability == [1.0, 0.0, 0.0, 0.0, 0.0]

    risk_adjusted_probability = Vector{Float64}(undef, 5)
    Kokako.adjust_probability(
        Kokako.WorstCase(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.3, 0.4, 0.5],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        false)
    @test risk_adjusted_probability == [0.0, 0.0, 0.0, 0.0, 1.0]
end
