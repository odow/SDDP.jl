#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

struct MyRiskMeasure <: SDDP.AbstractRiskMeasure end

@testset "Risk Measures" begin
    @testset "Expectation" begin
        ex = Expectation()
        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        SDDP.modifyprobability!(ex, y, x, obj)
        @test y == x
    end

    @testset "User-defined MyRiskMeasure" begin
        measure = MyRiskMeasure()
        y = zeros(4)
        x = [0.1, 0.2, 0.3, 0.4]
        obj = ones(4)
        @test_throws Exception SDDP.modifyprobability!(measure, y, x, obj)
    end
end
