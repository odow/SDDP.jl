#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# SDDP
# Stochastic Dual Dynamic Programming in Julia
# See http://github.com/odow/SDDP.jl
#############################################################################

@testset "Utilities" begin
    @test SDDP.getsense(SDDP.Max) == :Max
    @test SDDP.getsense(SDDP.Min) == :Min
    m=SDDP.Subproblem()
    @test SDDP.getsense(m) == :Min
    @test_throws Exception SDDP.dominates(:Minimize, 1.0, 2.0)
    @test_throws Exception SDDP.sample([0.0, 0.0])
    @test SDDP.rtol(1.0, 0.0) == 1.0
    @test SDDP.worstcase(:Min) == Inf
end


@testset "Print" begin
    @test SDDP.humanize(2000) == "  2.0K"
    @test SDDP.humanize(2000, "5.2f") == " 2.00K"
    @test_throws Exception SDDP.humanize(2000, "5.3f")
end
