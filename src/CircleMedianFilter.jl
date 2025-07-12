"""
    CircleMedianFilter

Implementation of Algorithm 1 from the research paper:
"Fast Median Filtering for Phase or Orientation Data"
by Martin Storath and Andreas Weinmann
IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 40, No. 3, March 2018
DOI: 10.1109/TPAMI.2017.2692779

This module provides functionality for computing median filters on circular data,
particularly useful for processing phase data and orientation data where the
circular nature must be preserved. The implementation focuses on the arc distance
median, which provides robust, edge-preserving and value-preserving smoothing
for circle-valued data.

Applications include: phase images from interferometric synthetic aperture radar,
planar flow fields from optical flow, wind direction time series, and other
data living on the unit circle.
"""
module CircleMedianFilter

"""
    d(a, b)

Calculate the angular distance between two angles on a circle.

This function computes the shortest angular distance between two angles `a` and `b`
(in radians) on a unit circle. The result is always between 0 and π radians,
representing the smallest arc length between the two angles.

# Arguments
- `a`: First angle in radians
- `b`: Second angle in radians

# Returns
- The shortest angular distance between `a` and `b` in radians (0 ≤ distance ≤ π)

# Examples
```julia
d(0, π/2)     # π/2 (90 degrees)
d(0, 3π/2)    # π/2 (90 degrees, shorter path)
d(π/4, 5π/4)  # π (180 degrees)
```
"""
d(a, b) = abs(rem2pi(a - b, RoundNearest)) # distance between two angles in radians


"""
    mirror_index(a, i, j)

Access array elements with mirrored boundary conditions for out-of-bounds indices.

This function provides symmetric mirroring boundary conditions when accessing array elements
with indices that fall outside the array bounds. For negative indices or indices exceeding
the array dimensions, it returns values from the array using mirror reflection at the boundaries.

## Arguments
- `a`: Input 2D array
- `i`: Row index (can be negative or exceed array bounds)
- `j`: Column index (can be negative or exceed array bounds)

## Returns
- Array element at the mirrored position

## Boundary Behavior
- For negative indices: `a[0]` → `a[1]`, `a[-1]` → `a[2]`, etc.
- For indices beyond bounds: `a[n+1]` → `a[n]`, `a[n+2]` → `a[n-1]`, etc.
- The mirroring is symmetric and preserves the structure of the data

## Examples
```julia
a = [1 2 3; 4 5 6; 7 8 9]  ## 3×3 array

## Normal access
mirror_index(a, 2, 2)  ## Returns 5 (same as a[2, 2])

## Negative indices (mirror at beginning)
mirror_index(a, 0, 2)   ## Returns 2 (mirrors to a[1, 2])
mirror_index(a, -1, 2)  ## Returns 5 (mirrors to a[2, 2])

## Indices beyond bounds (mirror at end)
mirror_index(a, 4, 2)   ## Returns 5 (mirrors to a[2, 2])
mirror_index(a, 5, 2)   ## Returns 2 (mirrors to a[1, 2])
```

## Notes
This function is particularly useful for implementing filters with symmetric boundary
conditions, ensuring that edge effects are minimized while maintaining the circular
or periodic nature of the data being processed.

See also: [`arc_distance_median_filter!`](@ref)
"""
function mirror_index(a, i, j)
    ## Get the size of the array
    n, m = size(a)
    ## Handle negative indices
    if i < 1
        i = 1 - i
    elseif i > n
        i = 2n + 1 - i
    end
    ## Handle indices larger than the size of the array
    if j < 1
        j = 1 - j
    elseif j > m
        j = 2m + 1 - j
    end
    return a[i, j]
end




"""
    arc_distance__median_filter!(u, y, r, t)

Compute the arc distance median filter on circular data.

# Arguments
- `u`: Filtered output array, must be the same size as `y`
- `y`: Input array containing circular data (angles in radians)
 - `r`: half-width of the filter, defining the neighborhood size R = 2r+1 in the first coordinate
 - `t`: half-width of the filter, defining the neighborhood size T = 2t+1 in the second coordinate



"""
function arc_distance__median_filter!(u, y, r::Integer, t::Integer)
    # The algorithm calcultes (using iterative relations) Gₘₙ -- a table of distances from the center pixel yₘₙ to all other pixel in the window R×T
    # mirroring of the data is used for the pixels outside the image
    #
    # For current values of m and n, after calculating Gₘₙ it is saved in the n-th  position of the Gbuf array.

    #check the sizes of the arrays
    size(u) == size(y) || error("Input and output arrays must have the same size")
    ndims(u) == 2 || error("Input and output arrays must be 2-dimensional")
    #check the radius
    r < 0 && error("Radius r must be non-negative")
    t < 0 && error("Radius t must be non-negative")
    #check the type of the arrays
    (eltype(u) <: Real) || error("Input and output arrays must be of a real type")
    (eltype(y) <: Real) || error("Input array must be of a real type")

    #calculate R and T
    R = 2r + 1
    T = 2t + 1




    #allocate buffer arrays for G and Z
    # each buffer has the length of the the second dimension of the input array
    # and its element is array RxT
    Gbuf = Array{Float64}(undef, size(y, 2), R, T)
    Zbuf = Array{Float64}(undef, size(y, 2), R, T)
    Gcurrent = Array{Float64}(undef, R, T)
    Zcurrent = Array{Float64}(undef, R, T)



    # Process the first column
    m = 1
    # the first G is calculated by the definition
    n = 1
    for i in (-r):r
        for j in (-t):t
            Gcurrent[i + r + 1, j + t + 1] = d(y[m, n], y[mirror_index(y, m + i, n + j)])
        end
    end
    #  save to the buffer:
    Gbuf[n, :, :] .= Gcurrent[:, :]

    # Process the rest of the first column:
    for n in 2:size(y, 2)
        # Fill the buffer Gbuf for the first column
        for i in (-r):r
            for j in (-t):t
                Gcurrent[i + r + 1, j + t + 1] =
                    Gbuf[n - 1, i + r + 1, j + 1 + t + 1] +
                    d(y[m, n], y[mirror_index(y, m + i, n + j)])
                d(y[m, n], y[mirror_index(y, m + i, n + j)])
            end
        end
        Gbuf[n, :, :] .= Gcurrent[:, :]
    end

    for n in 1:size(y, 2)
        # Fill the buffer Gbuf and Zbuf for the first column
        for i in 1:R
            for j in 1:T
                Gbuf[n, i, j] = d(y[m, n], y[m + i - 1, n])
                Zbuf[n, i, j] = y[m + i - 1, n]
            end
        end
    end


    # Process the other rows



end

end
