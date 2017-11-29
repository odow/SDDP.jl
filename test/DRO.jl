#  Copyright 2017,  Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

@testset "PopStd" begin
    x = rand(10)
    @test isapprox(SDDP.popstd(x), std(x, corrected=false), atol=1e-9)
end
@testset "UpdateProbabilities" begin
    S = 5
    r = 0.25
    theta = -float([2; 1; 3; 4; 5])
    newprobabilities = zeros(S)
    SDDP.modifyprobability!(newprobabilities, DRO(r), :Min, theta, S)
    @test isapprox(newprobabilities, [0.279057,0.358114,0.2,0.120943,0.0418861], atol=1e-6)
    SDDP.modifyprobability!(newprobabilities, DRO(r), :Max, -theta, S)
    @test isapprox(newprobabilities, [0.279057,0.358114,0.2,0.120943,0.0418861], atol=1e-6)
    r = 0.4
    SDDP.modifyprobability!(newprobabilities, DRO(r), :Min, theta, S)
    @test isapprox(newprobabilities, [0.324162,0.472486,0.175838,0.027514,0.0], atol=1e-6)
    SDDP.modifyprobability!(newprobabilities, DRO(r), :Max, -theta, S)
    @test isapprox(newprobabilities, [0.324162,0.472486,0.175838,0.027514,0.0], atol=1e-6)
    r = sqrt(0.8)
    SDDP.modifyprobability!(newprobabilities, DRO(r), :Min, theta, S)
    @test isapprox(newprobabilities, [0.0,1.0,0.0,0.0,0.0], atol=1e-6)
end
