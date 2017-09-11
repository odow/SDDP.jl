#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

@testset "SDDPModel" begin
    @testset "Test kwargs" begin
        # bad sense
        @test_throws Exception SDDPModel(sense=:Minimization, stages=3, objective_bound=10) do sp, t
        end
        # no bound
        @test_throws Exception SDDPModel(sense=:Min, stages=3) do sp, t
        end
        # too few args
        @test_throws Exception SDDPModel(sense=:Min, stages=3, objective_bound=10) do sp
        end
        # too many args
        @test_throws Exception SDDPModel(sense=:Min, stages=3, objective_bound=10) do sp, t, i, j
        end
        # markov but no argument
        @test_throws Exception SDDPModel(sense=:Min, stages=3, objective_bound=10, markov_transition=[0.5 0.5; 0.5 0.5]) do sp, t
        end
        # test sp is subproblem
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10) do sp, t
            @test SDDP.isext(sp)
        end
        @test length(m.stages) == 3
        # can't get bound of unsolved problem
        @test_throws Exception getbound(m)
    end

    @testset "Related Utilities" begin
        # test sp is subproblem
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10) do sp, t
            @state(sp, x>=0, x0==0.0)
        end
        sp = SDDP.getsubproblem(m, 1, 1)
        @test SDDP.cuttoaffexpr(sp, SDDP.Cut(1.0, [0.5])) == 1.0 + 0.5*sp[:x]
        # @test SDDP.cutoracle(sp) == SDDP.valueoracle(sp).cutmanager
    end
end
