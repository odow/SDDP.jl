# test the triangular dispatch
@testset "Triangular Dispatch" begin
    @test SDDP.getel(Int, 1, 1, 1) == 1
    @test_throws Exception SDDP.getel(Int, 1.0, 1, 1)

    @test SDDP.getel(Int, [2,1,3], 2, 1) == 1
    @test SDDP.getel(Int, [2,1,3], 2, 3) == 1
    @test_throws Exception SDDP.getel(Int, [2,1,3], 4, 3)

    @test SDDP.getel(Integer, Int16[2,1,3], 2, 1) == 1
    @test SDDP.getel(Integer, Int16[2,1,3], 3, 2) == 3
    @test_throws Exception SDDP.getel(Integer, Int16[2,1,3], 4, 3)

    @test SDDP.getel(Integer, [ [2], [1, 5], [3, 4, 2]], 2, 1) == 1
    @test_throws Exception SDDP.getel(Integer, [ [2], [1, 5], [3, 4, 2]], 1, 2)
    @test_throws Exception SDDP.getel(Integer, [ [2], [1, 5], [3, 4, 2]], 4, 1)
    @test SDDP.getel(Integer, [ [2], [1, 5], [3, 4, 2]], 3, 3) == 2
    @test SDDP.getel(Integer, [ Int16[2], [1, 5], Int32[3, 4, 2]], 2, 2) == 5
end

# issue 64
@testset "Stage dependent inputs via functions" begin
    @test SDDP.getel(Int, (t, i) -> t + i, 1, 1) == 2
    @test SDDP.getel(Int, (t, i) -> t + i, 2) == 3
    @test_throws Exception SDDP.getel(Float64, (t, i) -> t + i, 2)
    @test_throws Exception SDDP.getel(Float64, (t, i) -> t + i)
    @test_throws Exception SDDP.getel(Float64, (t) -> 2t, 1.0)
end
