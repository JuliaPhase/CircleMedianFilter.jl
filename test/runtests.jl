using CircleMedianFilter
using CircleMedianFilter: d, mirror_index, arc_distance__median_filter!

using Test

@testset "CircleMedianFilter.jl" begin
    # Write your tests here.
    @testset "Distance Function" begin
        # Test the distance function with some known values
        @test d(0.0, 0.0) == 0.0
        @test d(π / 2, -π / 2) ≈ π
        @test isapprox(d(π, -π), 0.0, atol=2eps(Float64))
        @test d(π / 4, 3π / 4) ≈ π / 2
        @test d(3π / 4, -3π / 4) ≈ π / 2
        @test isapprox(d(2π, 0), 0.0, atol=2eps(Float64))
    end

    @testset "Mirror Index Function" begin
        # Test the mirror index function with a sample array
        a = [1 2 3; 4 5 6; 7 8 9]

        @test mirror_index(a, 1, 1) == 1
        @test mirror_index(a, 0, -1) == 2
        @test mirror_index(a, 4, 2) == 8
        @test mirror_index(a, -2, 3) == 9
        @test mirror_index(a, 3, -2) == 9
    end
end
