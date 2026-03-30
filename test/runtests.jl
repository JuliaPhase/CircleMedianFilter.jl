using CircleMedianFilter
using CircleMedianFilter: arc_distance_median_filter!, d, mirror_index
using Test


@testitem "CircleMedianFilter.jl" begin
    using CircleMedianFilter
    using CircleMedianFilter: arc_distance_median_filter!, d, mirror_index

    using Test

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
        a = [1 2; 3 4]

        # Normal access
        @test mirror_index(a, 1, 1) == 1
        @test mirror_index(a, 2, 2) == 4

        # Mirroring at boundaries (virtual concatenation of flipped array)
        # For rows
        @test mirror_index(a, 0, 1) == 1   # row 0 -> row 1
        @test mirror_index(a, -1, 2) == 4  # row -1 -> row 1, col 2
        @test mirror_index(a, 4, 1) == 1   # row 3 -> row 1

        # For columns
        @test mirror_index(a, 1, 0) == 1   # col 0 -> col 1
        @test mirror_index(a, 2, -1) == 4  # col -1 -> col 2
    end

    @testset "Arc Distance Median Filter" begin
        @testset "Input Validation" begin
            # Test size mismatch
            y = rand(3, 3) * 2π
            u_wrong = rand(2, 2)
            @test_throws ErrorException arc_distance_median_filter!(u_wrong, y, 1, 1)

            # Test negative radius
            u = similar(y)
            @test_throws ErrorException arc_distance_median_filter!(u, y, -1, 1)
            @test_throws ErrorException arc_distance_median_filter!(u, y, 1, -1)

            # Test non-2D arrays
            y_1d = rand(10)
            u_1d = similar(y_1d)
            @test_throws ErrorException arc_distance_median_filter!(u_1d, y_1d, 1, 1)
        end

        @testset "Basic Functionality" begin
            # Test with small uniform array (should preserve values)
            y_uniform = fill(π / 4, 5, 5)
            u = similar(y_uniform)
            arc_distance_median_filter!(u, y_uniform, 1, 1)
            @test all(u .≈ π / 4)

            # Test with zero radius (should copy input)
            y = rand(5, 5) * 2π
            u = similar(y)
            arc_distance_median_filter!(u, y, 0, 0)
            @test u ≈ y
        end

        @testset "Noise Filtering" begin
            # Simple test: uniform data with one outlier should be corrected
            y = fill(π / 4, 5, 5)
            y[3, 3] = π  # Single outlier

            u = similar(y)
            arc_distance_median_filter!(u, y, 1, 1)

            # The outlier should be replaced with something closer to neighbors
            @test u[3, 3] == π / 4
        end

        @testset "Edge Cases" begin
            # Test with 1x1 array
            y_small = [π / 3;;]
            u_small = similar(y_small)
            arc_distance_median_filter!(u_small, y_small, 1, 1)
            @test u_small[1, 1] ≈ π / 3
        end

        @testset "Circular Data Properties" begin
            # Test that filter respects circular nature of data
            # Create data near 0/2π boundary
            y = fill(0.1, 5, 5)
            y[3, 3] = 2π - 0.1  # Value near 2π, should be close to others in circular space

            u = similar(y)
            arc_distance_median_filter!(u, y, 1, 1)

            # The filtered center value should be close to neighbors in circular space
            # Either close to 0.1 or close to 2π-0.1, but not somewhere in between like π
            center_dist_to_01 = d(u[3, 3], 0.1)
            center_dist_to_2pi_01 = d(u[3, 3], 2π - 0.1)
            center_dist_to_pi = d(u[3, 3], π)

            @test min(center_dist_to_01, center_dist_to_2pi_01) < center_dist_to_pi
        end

        @testset "In-place Modification" begin
            # Test that function modifies input array u in place
            y = rand(5, 5) * 2π
            u = rand(5, 5) * 2π  # Initialize with different values
            u_original = copy(u)

            result = arc_distance_median_filter!(u, y, 1, 1)

            # Check that u was modified
            @test u != u_original
            # Check that returned value is the same object as u
            @test result === u
        end
    end
end
