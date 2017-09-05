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
            @noise(sp, w=1:2, x <= w)
            stageobjective!(sp, x)
        end

        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
        @test y == x
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
            @noise(sp, w=1:2, x <= w)
            stageobjective!(sp, x)
        end

        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        @test_throws Exception SDDP.modifyprobability!(measure, y, x, obj, m, SDDP.getsubproblem(m, 1, 1))
    end
end
