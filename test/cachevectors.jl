#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

@testset "Cached Vectors" begin
    y = SDDP.CachedVector(Int)
    for i in [2, 1, 4]
        push!(y, i)
    end
    @test y[2] == 1
    SDDP.reset!(y)
    @test length(y) == 0
    @test_throws BoundsError y[1]
    @assert length(y.data) == 3
    @test_throws BoundsError y[1] = 2
    push!(y, 5)
    @test y[1] == 5
end
