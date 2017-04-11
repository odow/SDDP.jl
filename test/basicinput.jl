#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################
using MathProgBase, Clp

@testset "Deterministic 2 stage problem" begin

    m = SDDPModel(stages=2, sense=:Min, solver=ClpSolver()) do subproblem, stage
        @state(subproblem, x >= 0, x0 == 5)
        @variable(subproblem, 0 <= u <= 2)
        @constraint(subproblem, x == x0 - u)
        @stageobjective(subproblem, stage * u)
    end
    @testset "SDDPModel fields" begin
        @test length(stages(m)) == 2
    end
    @testset "Stages" begin
        for (i, stage) in enumerate(stages(m))
            @test length(stage.subproblems) == 1
            @test length(stages.transitionprobabilities) == 1
            @test_approx_eq stages.transitionprobabilities[1] 1.0

            @testset "Subproblems" begin
                for subproblem in subproblems(stage)
                    @test subproblem.stage == i
                    @test subproblem.markovstate == 1
                    @test length(subproblem.states) == 1
                    @test length(subproblem.scenarios) == 0
                    @test subproblem.riskmeasure == Expectation()
                    @test getsense(subproblem.m) == :Min
                end
            end
        end
    end
end

@testset "Scenario constraints" begin

    m = SDDPModel(stages=2, sense=:Min) do subproblem, stage
        @state(subproblem, x >= 0, x0 == 5)
        @variable(subproblem, 0 <= u <= 2)
        @scenarioconstraint(subproblem, i=[-0.5, 1.0], x == x0 - u + i)
        @stageobjective(subproblem, stage * u)
    end
    @testset "SDDPModel fields" begin
        @test length(stages(m)) == 2
    end
    @testset "Stages" begin
        for (i, stage) in enumerate(stages(m))
            @test length(stage.subproblems) == 1
            @test length(stages.transitionprobabilities) == 1
            @test_approx_eq stages.transitionprobabilities[1] 1.0

            @testset "Subproblems" begin
                for (j, subproblem) in enumerate(subproblems(stage))
                    @test subproblem.stage == i
                    @test subproblem.markovstate == 1
                    @test length(subproblem.states) == 1
                    @test length(subproblem.scenarios) == 2
                    @testset "Stage $i, subproblem $j: scenarios" begin
                        @test length(subproblem.scenarios.scenarios) == 2
                        @test length(subproblem.scenarios.probability) == 2
                        @test subproblem.scenarios.probability == [0.5, 0.5]
                        @test subproblem.scenarios.scenarios[1].row == [2]
                        @test subproblem.scenarios.scenarios[1].value == [-0.5]
                        @test subproblem.scenarios.scenarios[2].row == [2]
                        @test subproblem.scenarios.scenarios[1].value == [1.0]
                    end
                    @test subproblem.riskmeasure == Expectation()
                    @test getsense(subproblem.m) == :Min
                end
            end
        end
    end
end
