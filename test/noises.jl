#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

@testset "Noises" begin
    @testset "Simple RHS" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        c = @noise(m, i=1:2, x <= i)
        noises = SDDP.ext(m).noises
        @test length(noises) == 2
        @test noises[1].constraints == [c]
        @test noises[1].values == [1]
        @test noises[2].constraints == [c]
        @test noises[2].values == [2]
    end

    @testset "Function RHS" begin
        Ω = [1.1, 1.2, 1.3]
        f = (w) -> 2w
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        c = @noise(m, i=Ω, x <= f(i))
        noises = SDDP.ext(m).noises
        @test length(noises) == length(Ω)
        for (i, ω) in enumerate(Ω)
            @test noises[i].constraints == [c]
            @test noises[i].values == [f(ω)]
        end
    end

    @testset "Incompatible Noise lengths" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        c = @noise(m, i=1:2, x <= i)
        @test_throws AssertionError @noise(m, i=1:3, x >= -i)
    end

    @testset "Noises" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        @noises(m, i=1:2, begin
            x <= i
            x >= -i
        end)
        noises = SDDP.ext(m).noises
        @test length(noises) == 2
        @test noises[1].values == [1, -1]
        @test noises[2].values == [2, -2]
    end

    @testset "Noise probability" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        setnoiseprobability!(m, [0.1, 0.2, 0.7])
        @test SDDP.ext(m).noiseprobability == [0.1, 0.2, 0.7]
        @test_throws AssertionError setnoiseprobability!(m, [0.1, 0.2])
    end
end
