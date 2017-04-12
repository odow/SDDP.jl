#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

@testset "Scenarios" begin
    @testset "Simple RHS" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        c = @scenario(m, i=1:2, x <= i)
        scenarios = SDDP.ext(m).scenarios
        @test length(scenarios) == 2
        @test scenarios[1].constraints == [c]
        @test scenarios[1].values == [1]
        @test scenarios[2].constraints == [c]
        @test scenarios[2].values == [2]
    end

    @testset "Function RHS" begin
        Ω = [1.1, 1.2, 1.3]
        f = (w) -> 2w
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        c = @scenario(m, i=Ω, x <= f(i))
        scenarios = SDDP.ext(m).scenarios
        @test length(scenarios) == length(Ω)
        for (i, ω) in enumerate(Ω)
            @test scenarios[i].constraints == [c]
            @test scenarios[i].values == [f(ω)]
        end
    end

    @testset "Incompatible Scenario lengths" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        c = @scenario(m, i=1:2, x <= i)
        @test_throws AssertionError @scenario(m, i=1:3, x >= -i)
    end

    @testset "Scenarios" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        @variable(m, x)
        @scenarios(m, i=1:2, begin
            x <= i
            x >= -i
        end)
        scenarios = SDDP.ext(m).scenarios
        @test length(scenarios) == 2
        @test scenarios[1].values == [1, -1]
        @test scenarios[2].values == [2, -2]
    end

    @testset "Scenario probability" begin
        m = SDDP.Subproblem()
        subproblem = SDDP.ext(m)
        setscenarioprobability!(m, [0.1, 0.2, 0.7])
        @test SDDP.ext(m).scenarioprobability == [0.1, 0.2, 0.7]
        @test_throws AssertionError setscenarioprobability!(m, [0.1, 0.2])
    end
end
