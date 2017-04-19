#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

@testset "Stage Objectives" begin
    @testset "Default Price Oracle" begin
        m = SDDP.Subproblem()
        @variable(m, x)
        stageobjective!(m, x)
        @test getobjective(m) == convert(QuadExpr, x + SDDP.ext(m).valueoracle.theta)
    end
end
