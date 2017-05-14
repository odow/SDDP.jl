#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

immutable NewCutOracle <: SDDP.AbstractCutOracle end
immutable NotACutOracle end

@testset "Cut Oracles" begin
    @testset "New Cut Oracle" begin
        newcutoracle = NewCutOracle()
        cut = SDDP.Cut(1, [1])
        @test_throws Exception SDDP.storecut!(newcutoracle, cut)
        @test_throws Exception SDDP.validcuts(newcutoracle)
    end

    @testset "Not A Cut Oracle" begin
        notacutoracle = NotACutOracle()
        cut = SDDP.Cut(1, [1])
        @test_throws Exception SDDP.storecut!(notacutoracle, cut)
        @test_throws Exception SDDP.validcuts(notacutoracle)
    end

    @testset "Default Cut Oracle" begin
        oracle = SDDP.DefaultCutOracle()
        cut = SDDP.Cut(1, [1])
        SDDP.storecut!(oracle, cut)
        @test SDDP.validcuts(oracle) == [cut]
    end
end
