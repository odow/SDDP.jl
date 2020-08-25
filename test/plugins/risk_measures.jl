#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using GLPK
using SDDP
using Test

@testset "Expectation" begin
    @test sprint(show, SDDP.Expectation()) == "SDDP.Expectation()"
    risk_adjusted_probability = Vector{Float64}(undef, 5)
    SDDP.adjust_probability(
        SDDP.Expectation(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.3, 0.4, 0.5],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        true,
    )
    @test risk_adjusted_probability == [0.1, 0.2, 0.3, 0.4, 0.5]

    risk_adjusted_probability = Vector{Float64}(undef, 5)
    SDDP.adjust_probability(
        SDDP.Expectation(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.3, 0.3, 0.1],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        false,
    )
    @test risk_adjusted_probability == [0.1, 0.2, 0.3, 0.3, 0.1]
end

@testset "WorstCase" begin
    risk_adjusted_probability = Vector{Float64}(undef, 5)
    SDDP.adjust_probability(
        SDDP.WorstCase(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.0, 0.4, 0.5],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        true,
    )
    @test risk_adjusted_probability == [1.0, 0.0, 0.0, 0.0, 0.0]

    risk_adjusted_probability = Vector{Float64}(undef, 5)
    SDDP.adjust_probability(
        SDDP.WorstCase(),
        risk_adjusted_probability,
        [0.1, 0.2, 0.3, 0.4, 0.5],
        [:a, :b, :c, :d, :e],
        [5.0, 4.0, 6.0, 2.0, 1.0],
        false,
    )
    @test risk_adjusted_probability == [0.0, 0.0, 0.0, 0.0, 1.0]
end

@testset "Constructors" begin
    a = SDDP.Expectation()
    b = SDDP.AVaR(0.5)
    c = SDDP.WorstCase()
    d = 0.5a + 0.3b + 0.2c
    @test d.measures[1] == (0.5, a)
    @test d.measures[2] == (0.3, b)
    @test d.measures[3] == (0.2, c)

    aa = SDDP.EAVaR(lambda = 0.5, beta = 0.25)
    @test aa.measures[1] == (0.5, SDDP.Expectation())
    @test aa.measures[2] == (0.5, SDDP.AVaR(0.25))
end

@testset "AV@R" begin
    @test_throws Exception SDDP.AVaR(-0.1)
    @test_throws Exception SDDP.AVaR(1.1)
    @testset "beta=0.2" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.AVaR(0.2),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            false,
        )
        @test risk_adjusted_probability == [0.5, 0.5, 0.0, 0.0]
    end
    @testset "beta=0" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.AVaR(0.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            false,
        )
        @test risk_adjusted_probability == [1.0, 0.0, 0.0, 0.0]
    end
    @testset "beta=1" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.AVaR(1.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            false,
        )
        @test risk_adjusted_probability == [0.1, 0.2, 0.3, 0.4]
    end
end

@testset "EAV@R" begin
    @test sprint(show, SDDP.EAVaR(lambda = 0.2, beta = 0.3)) ==
          "A convex combination of 0.2 * SDDP.Expectation() + 0.8 * SDDP.AVaR(0.3)"
    @test_throws Exception SDDP.EAVaR(lambda = 1.1)
    @test_throws Exception SDDP.EAVaR(lambda = -0.1)
    @test_throws Exception SDDP.EAVaR(beta = 1.1)
    @test_throws Exception SDDP.EAVaR(beta = -0.1)
    @testset "Max - (0.25, 0.2)" begin
        nominal_probability = [0.1, 0.2, 0.3, 0.4]
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.EAVaR(lambda = 0.25, beta = 0.2),
            risk_adjusted_probability,
            nominal_probability,
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            false,
        )
        @test risk_adjusted_probability ≈
              0.25 * nominal_probability + 0.75 * [1 / 2, 1 / 2, 0, 0]
    end
    @testset "Min - (0.25, 0.2)" begin
        nominal_probability = [0.1, 0.2, 0.3, 0.4]
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.EAVaR(lambda = 0.25, beta = 0.2),
            risk_adjusted_probability,
            nominal_probability,
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            true,
        )
        @test risk_adjusted_probability ≈ 0.25 * nominal_probability + 0.75 * [0, 0, 0, 1.0]
    end
    @testset "Max - (0.5, 0.0)" begin
        nominal_probability = [0.1, 0.2, 0.3, 0.4]
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.EAVaR(lambda = 0.5, beta = 0.0),
            risk_adjusted_probability,
            nominal_probability,
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            false,
        )
        @test risk_adjusted_probability ≈ 0.5 * nominal_probability + 0.5 * [1.0, 0, 0, 0]
    end
    @testset "Max - (0.5, 0.0) 2" begin
        nominal_probability = [0.0, 0.2, 0.4, 0.4]
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            SDDP.EAVaR(lambda = 0.5, beta = 0.0),
            risk_adjusted_probability,
            nominal_probability,
            [:a, :b, :c, :d],
            [1.0, 2.0, 3.0, 4.0],
            false,
        )
        @test risk_adjusted_probability ≈ 0.5 * nominal_probability + 0.5 * [0.0, 1.0, 0, 0]
    end
end

@testset "ModifiedChiSquared" begin
    @test sprint(show, SDDP.ModifiedChiSquared(0.1)) == "ModifiedChiSquared with radius=0.1"
    @testset "Min - R=0.0" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.0),
            risk_adjusted_probability,
            fill(0.2, 5),
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            true,
        )
        @test risk_adjusted_probability ≈ [0.2, 0.2, 0.2, 0.2, 0.2] atol = 1e-6
    end
    @testset "Non-uniform Expectation" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.2, 0.2],
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            true,
        )
        @test risk_adjusted_probability ≈ [0.1, 0.2, 0.3, 0.2, 0.2] atol = 1e-6
    end
    @testset "Non-uniform Expectation Max" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.2, 0.2],
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            false,
        )
        @test risk_adjusted_probability ≈ [0.1, 0.2, 0.3, 0.2, 0.2] atol = 1e-6
    end

    @testset "Non-uniform Worst Case" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(6.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.3, 0.1],
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            true,
        )
        @test risk_adjusted_probability ≈ [0.0, 1.0, 0.0, 0.0, 0.0] atol = 1e-6
    end

    @testset "Non-uniform Worst Case Max" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(6.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.3, 0.1],
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            false,
        )
        @test risk_adjusted_probability ≈ [0.0, 0.0, 0.0, 0.0, 1.0] atol = 1e-6
    end

    @testset "Non-uniform" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.45),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.3, 0.1],
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -0.5],
            true,
        )
        @test risk_adjusted_probability ≈ [0.115714, 0.372861, 0.158568, 0.001421,0.351435] atol = 1e-6
    end

    @testset "Non-uniform Max" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.45),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.3, 0.1],
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -0.5],
            false,
        )
        @test risk_adjusted_probability ≈ [0.0, 0.0, 0.323223, 0.676777, 0.0] atol = 1e-6
    end
    @testset "Min - R=0.25" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.25),
            risk_adjusted_probability,
            fill(0.2, 5),
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            true,
        )
        @test risk_adjusted_probability ≈ [0.279057, 0.358114, 0.2, 0.120943, 0.0418861] atol =
            1e-6
    end
    @testset "Max - R=0.25" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.25),
            risk_adjusted_probability,
            fill(0.2, 5),
            [:a, :b, :c, :d, :e],
            [2.0, 1.0, 3.0, 4.0, 5.0],
            false,
        )
        @test risk_adjusted_probability ≈ [0.279057, 0.358114, 0.2, 0.120943, 0.0418861] atol =
            1e-6
    end
    @testset "Min - R=0.4" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.4),
            risk_adjusted_probability,
            fill(0.2, 5),
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            true,
        )
        @test risk_adjusted_probability ≈ [0.324162, 0.472486, 0.175838, 0.027514, 0.0] atol =
            1e-6
    end
    @testset "Max - R=0.4" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(0.4),
            risk_adjusted_probability,
            fill(0.2, 5),
            [:a, :b, :c, :d, :e],
            [2.0, 1.0, 3.0, 4.0, 5.0],
            false,
        )
        @test risk_adjusted_probability ≈ [0.324162, 0.472486, 0.175838, 0.027514, 0.0] atol =
            1e-6
    end
    @testset "Min - R=√0.8" begin
        risk_adjusted_probability = Vector{Float64}(undef, 5)
        SDDP.adjust_probability(
            SDDP.ModifiedChiSquared(sqrt(0.8)),
            risk_adjusted_probability,
            fill(0.2, 5),
            [:a, :b, :c, :d, :e],
            [-2.0, -1.0, -3.0, -4.0, -5.0],
            true,
        )
        @test risk_adjusted_probability ≈ [0, 1.0, 0, 0, 0]
    end
end

@testset "Wasserstein" begin
    function default_wasserstein(alpha)
        return SDDP.Wasserstein(GLPK.Optimizer; alpha = alpha) do x, y
            return abs(x - y)
        end
    end
    @test sprint(show, default_wasserstein(0.1)) == "SDDP.Wasserstein"
    @test_throws Exception default_wasserstein(-1.0)
    @testset ":Max Worst case" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            default_wasserstein(10.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.3, 0.6, 0.4],
            [1.1, 1.2, 0.6, 1.3],
            false,
        )
        @test risk_adjusted_probability ≈ [0, 0, 1.0, 0]
    end
    @testset ":Min Worst case" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            default_wasserstein(10.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.3, 0.6, 0.4],
            [1.1, 1.2, 0.6, 1.3],
            true,
        )
        @test risk_adjusted_probability ≈ [0, 0, 0, 1.0]
    end
    @testset ":Max Expectation" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            default_wasserstein(0.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.3, 0.6, 0.4],
            [1.1, 1.2, 0.6, 1.3],
            false,
        )
        @test risk_adjusted_probability ≈ [0.1, 0.2, 0.3, 0.4]
    end
    @testset ":Min Expectation" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            default_wasserstein(0.0),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.3, 0.6, 0.4],
            [1.1, 1.2, 0.6, 1.3],
            true,
        )
        @test risk_adjusted_probability ≈ [0.1, 0.2, 0.3, 0.4]
    end
    @testset ":Max Intermediate" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            default_wasserstein(0.1),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.3, 0.6, 0.4],
            [1.1, 1.2, 0.6, 1.3],
            false,
        )
        @test risk_adjusted_probability ≈ [0.0, 1 / 6, 5 / 6, 0.0]
    end
    @testset ":Min Intermediate" begin
        risk_adjusted_probability = Vector{Float64}(undef, 4)
        SDDP.adjust_probability(
            default_wasserstein(0.1),
            risk_adjusted_probability,
            [0.1, 0.2, 0.3, 0.4],
            [0.5, 0.3, 0.6, 0.4],
            -[1.1, 1.2, 0.6, 1.3],
            true,
        )
        @test risk_adjusted_probability ≈ [0.0, 1 / 6, 5 / 6, 0.0]
    end
end

@testset "Entropic" begin
    @test sprint(show, SDDP.Entropic(0.1)) == "Entropic risk measure with γ = 0.1"
    @testset "Minimization" begin
        # Test that increasing values of θ lead to larger values for F[X].
        X = [1.0, 2.0, 3.0]
        P = [0.5, 0.5, 0.0]
        Q = [NaN, NaN, NaN]
        last, last_q2 = -Inf, 0.0
        for i = -4:4
            θ = 10.0^i
            α = SDDP.adjust_probability(SDDP.Entropic(θ), Q, P, [], X, true)
            current = Q' * X + α
            @test current >= last
            @test Q[2] >= last_q2
            last, last_q2 = current, Q[2]
        end
    end
    @testset "Maximization" begin
        # Test that increasing values of θ lead to smaller values for F[X].
        X = [1.0, 2.0, 3.0]
        P = [0.5, 0.5, 0.0]
        Q = [NaN, NaN, NaN]
        last, last_q1 = Inf, 0.0
        for i = -4:4
            θ = 10.0^i
            α = SDDP.adjust_probability(SDDP.Entropic(θ), Q, P, [], X, false)
            current = Q' * X + α
            if Q[1] < 1 - 1e-6
                # If Q[1] ≈ 1, this test can fail due to numerical error.
                @test current <= last
            end
            @test Q[1] >= last_q1
            last, last_q1 = current, Q[1]
        end
    end
end
