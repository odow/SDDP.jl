#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

immutable MyRiskMeasure <: SDDP.AbstractRiskMeasure end

@testset "Risk Measures" begin
    @testset "Expectation" begin
        measure = Expectation()
        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            risk_measure    = measure
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            stageobjective!(sp, x)
        end

        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
        @test y == x
    end

    @testset "AV@R" begin
        @test_throws Exception NestedAVaR(lambda=1.1)
        @test_throws Exception NestedAVaR(lambda=-0.1)
        @test_throws Exception NestedAVaR(beta=1.1)
        @test_throws Exception NestedAVaR(beta=-0.1)
        measure = NestedAVaR(lambda=0.25, beta=0.2)
        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            risk_measure    = measure
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end

        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = [1.0,2.0,3.0,4.0]
        SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
        @test y == 0.25 * x + 0.75 * [1/2, 1/2, 0, 0]

        measure = NestedAVaR(lambda=0.25, beta=0.2)
        m = SDDPModel(
            sense           = :Min,
            stages          = 2,
            objective_bound = 10,
            risk_measure    = measure
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end

        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = [1.0,2.0,3.0,4.0]
        SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
        @test y == 0.25 * x + 0.75 * [0, 0, 0, 1.0]

        measure = NestedAVaR(lambda=0.5, beta=0.0)
        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            risk_measure    = measure
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end
        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = [1.0,2.0,3.0,4.0]
        SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
        @test y == 0.5 * x + 0.5 * [1.0, 0, 0, 0]

        measure = NestedAVaR(lambda=0.5, beta=0.0)
        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            risk_measure    = measure
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end
        y = zeros(4)
        x = [0.0, 0.2, 0.4, 0.4]
        obj = [1.0,2.0,3.0,4.0]
        SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
        @test y == 0.5 * x + 0.5 * [0.0, 1.0, 0, 0]
    end

    @testset "User-defined MyRiskMeasure" begin
        measure = MyRiskMeasure()

        m = SDDPModel(
            sense           = :Max,
            stages          = 2,
            objective_bound = 10,
            risk_measure    = measure
            ) do sp, t
            @state(sp, x>=0, x0==0)
            @rhsnoise(sp, w=1:2, x <= w)
            @stageobjective(sp, x)
        end

        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        @test_throws Exception SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
    end
end
