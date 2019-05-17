#  Copyright 2017-19, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP, Test

@testset "_assign_to_centroid" begin
    centroids = [0.0, 1.0]
    samples = [-0.5, 1.1, -0.1, 0.9, 0.0, 0.49, 0.51]
    assignments = zeros(Int, 7)
    SDDP._assign_to_centroid(samples, centroids, assignments)
    @test assignments == [1, 2, 1, 2, 1, 1, 2]
end

@testset "_k_means" begin
    samples = [1.0, 2.0]
    centroids, score = SDDP._k_means(samples, 1)
    @test centroids == [1.5]
    centroids, score = SDDP._k_means(samples, 2)
    if score > 0.1
        @test centroids[1] == 1.5
        @test isnan(centroids[2])
        @test score == 0.25
    else
        @test centroids == [1.0, 2.0]
        @test score == 0.0
    end
end

@testset "Fitted Markov" begin
    graph = SDDP.MarkovianGraph(
            stages = 3, n_simulations = 1_000, n_markov_states=2) do t, x
        if t == 1
            return rand([-0.5, 0.5])
        end
        return 0.5 * x + rand([-0.25, 0.25])
    end
end
