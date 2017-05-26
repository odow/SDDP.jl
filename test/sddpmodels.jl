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
        @test_throws Exception SDDPModel(sense=:Minimisation) do sp, t
        end
        m = SDDPModel(sense=:Max, stages=3, objective_bound=10) do sp, t
        end
        @test length(m.stages) == 3
    end
end
