#  Copyright 2017-20, Oscar Dowson, Lea Kapelevich.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using SDDP: binexpand, bincontract
using Test

@testset "Binary Expansion" begin
    int_len = floor(Int, log(typemax(Int)) / log(2))
    @test_throws Exception binexpand(0)
    @test binexpand(1, 1) == [1]
    @test binexpand(2, 2) == [0, 1]
    @test binexpand(3, 3) == [1, 1]
    @test binexpand(4, 4) == [0, 0, 1]
    @test binexpand(5, 5) == [1, 0, 1]
    @test binexpand(6, 6) == [0, 1, 1]
    @test binexpand(7, 7) == [1, 1, 1]
    @test_throws Exception binexpand(8, 7)
    @test binexpand(typemax(Int), typemax(Int)) == ones(Int, int_len)
    @test binexpand(0.5, 0.5) == binexpand(5, 5)
    @test binexpand(0.54, 0.54) == binexpand(5, 5)
    @test binexpand(0.56, 0.56, 0.1) == binexpand(6, 6)
    @test binexpand(0.5, 0.5, 0.01) == binexpand(50, 50)
    @test binexpand(0.54, 0.54, 0.01) == binexpand(54, 54)
    @test binexpand(0.56, 0.56, 0.01) == binexpand(56, 56)

    @test 0 == bincontract([0])
    @test 1 == bincontract([1])
    @test 0 == bincontract([0, 0])
    @test 1 == bincontract([1, 0])
    @test 2 == bincontract([0, 1])
    @test 3 == bincontract([1, 1])
    @test 2 == bincontract([0, 1, 0])
    @test 3 == bincontract([1, 1, 0])
    @test 4 == bincontract([0, 0, 1])
    @test 5 == bincontract([1, 0, 1])
    @test 6 == bincontract([0, 1, 1])
    @test 7 == bincontract([1, 1, 1])
    @test typemax(Int) == bincontract(ones(Int, int_len))
    @test bincontract([0], 0.1) ≈ 0.0
    @test bincontract([1], 0.1) ≈ 0.1
    @test bincontract([0, 1], 0.1) ≈ 0.2
    @test bincontract([1, 1], 0.1) ≈ 0.3
    @test bincontract([0, 1, 0], 0.1) ≈ 0.2
    @test bincontract([1, 1, 0], 0.1) ≈ 0.3
    @test bincontract([1, 0, 1], 0.1) ≈ 0.5
end
