#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using Clp

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
        # no solver
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10) do sp, t
            @state(sp, x>=0, x0==0.0)
            @rhsnoise(sp, i=1:2, x <= i)
            @stageobjective(sp, i=1:2, i * x)
        end
        @test_throws Exception solve(m, print_level=0)

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
        @test SDDP.cutoracle(sp) == SDDP.valueoracle(sp).cutmanager
    end

    @testset "save/load" begin
        # test sp is subproblem
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10) do sp, t
            @state(sp, x>=0, x0==0.0)
            @rhsnoise(sp, i=1:2, x <= i)
            @stageobjective(sp, i=1:2, i * x)
        end
        SDDP.savemodel!("test.model", m)
        m2 = SDDP.loadmodel("test.model")
        @test length(m.stages) == length(m2.stages)
        for t in 1:3
            @test m.stages[t].transitionprobabilities == m2.stages[t].transitionprobabilities
            @test m.stages[t].ext == m2.stages[t].ext
            @test m.stages[t].ext == m2.stages[t].ext
            ex1 = m.stages[t].subproblems[1].ext[:SDDP]
            ex2 = m2.stages[t].subproblems[1].ext[:SDDP]
            @test length(ex1.noises) == length(ex2.noises)
            @test length(ex1.states) == length(ex2.states)
            @test length(ex1.problembound) == length(ex2.problembound)
            @test ex1.markovstate == ex2.markovstate
            @test ex1.stage == ex2.stage
            @test ex1.noiseprobability == ex2.noiseprobability
        end
        rm("test.model")
    end

    @testset "time limit" begin
        # test sp is subproblem
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10, solver=ClpSolver()) do sp, t
            @state(sp, x>=0, x0==0.0)
            @rhsnoise(sp, i=1:2, x <= i)
            @stageobjective(sp, i=1:2, i * x)
        end
        status = solve(m, time_limit=0.0, print_level=0)
        @test status == :time_limit
    end

    @testset "memory footprint" begin
        # test sp is subproblem
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10, solver=ClpSolver()) do sp, t
            @state(sp, x>=0, x0==0.0)
            @rhsnoise(sp, i=1:2, x <= i)
            @stageobjective(sp, i=1:2, i * x)
        end
        status = solve(m, max_iterations=1, print_level=0, reduce_memory_footprint=true)
        @test status == :max_iterations
        sp = SDDP.getsubproblem(m, 1, 1)
        for con in sp.linconstr
            @test con.terms == JuMP.AffExpr(0.0)
        end
    end
end
