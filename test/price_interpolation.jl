#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

@testset "DiscreteDistribution" begin
    @test_throws Exception DiscreteDistribution([1,2], [0.2, 0.9])
end

# @testset "StaticPriceInterpolation" begin
#     @test_throws Exception StaticPriceInterpolation(dynamics=(p,w,t,i)->0.0, initial_price=0.0, rib_locations=[0.0, 1.0])
#     @test_throws Exception StaticPriceInterpolation(dynamics=(p,w,t,i)->0.0, initial_price=0.0, noise=DiscreteDistribution([1.0]))
#     @test_throws Exception StaticPriceInterpolation(dynamics=(p,w,t,i)->0.0, rib_locations=[0.0, 1.0], noise=DiscreteDistribution([1.0]))
#     @test_throws Exception StaticPriceInterpolation(initial_price=0.0, rib_locations=[0.0, 1.0], noise=DiscreteDistribution([1.0]))
#     @test_throws Exception StaticPriceInterpolation(dynamics=(p,w,t,i)->0.0, initial_price=(1,1), rib_locations=[0.0, 1.0], noise=DiscreteDistribution([1.0]))
# end
