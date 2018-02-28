#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

struct NewCutOracle <: SDDP.AbstractCutOracle end
struct NotACutOracle end

@testset "Cut Oracles" begin
    @testset "New Cut Oracle" begin
        newcutoracle = NewCutOracle()
        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            cut_oracle = newcutoracle
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end
        cut = SDDP.Cut(1, [1])
        @test_throws Exception SDDP.storecut!(newcutoracle, m, SDDP.getsubproblem(m, 1, 1), cut)
        @test_throws Exception SDDP.validcuts(newcutoracle)
    end

    @testset "Not A Cut Oracle" begin
        notacutoracle = NotACutOracle()
        @test_throws Exception SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            cut_oracle = notacutoracle
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end
    end

    @testset "Default Cut Oracle" begin
        oracle = SDDP.DefaultCutOracle()
        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            cut_oracle = oracle
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end
        cut = SDDP.Cut(1, [1])
        SDDP.storecut!(oracle, m, SDDP.getsubproblem(m, 1, 1), cut)
        @test SDDP.validcuts(oracle) == [cut]
    end

    @testset "Level One Oracle" begin
        oracle = SDDP.LevelOneCutOracle()
    end
end
