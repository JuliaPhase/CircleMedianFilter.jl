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
export arc_distance_median_filter!

using OffsetArrays

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

Compute the arc distance median filter on circular data using Algorithm 1 from
Storath & Weinmann (2018).

This function implements an efficient algorithm for computing median filters on circular
data that preserves the circular nature of the input. The algorithm uses iterative
relations to compute distance tables efficiently, making it suitable for large images
with circular/angular data such as phase images or orientation fields.

## Algorithm Details

The function computes for each pixel the circular median within a (2r+1) × (2t+1)
rectangular neighborhood. The circular median minimizes the sum of arc distances
to all pixels in the neighborhood, providing robust edge-preserving smoothing.

The implementation uses:
- Mirrored boundary conditions via `mirror_index`
- Iterative computation of distance tables G_{m,n} for efficiency
- OffsetArrays for natural indexing with negative indices

## Arguments
- `u`: Pre-allocated output array for filtered results (modified in-place)
- `y`: Input 2D array containing circular data (angles in radians)
- `r::Integer`: Half-width of filter in first dimension (rows), so neighborhood size is 2r+1
- `t::Integer`: Half-width of filter in second dimension (columns), so neighborhood size is 2t+1

## Returns
- `u`: The filtered output array (same object as input `u`)

## Examples
```julia
using CircleMedianFilter
using OffsetArrays

## Create test phase data with some noise
y = rand(50, 50) * 2π  ## Random phase data
u = similar(y)         ## Pre-allocate output

## Apply 3×3 median filter (r=1, t=1)
arc_distance_median_filter!(u, y, 1, 1)

## Apply 5×7 median filter (r=2, t=3)
arc_distance_median_filter!(u, y, 2, 3)
```

## Notes
- Input angles should be in radians
- The function modifies the output array `u` in-place for memory efficiency
- Mirrored boundary conditions are used to handle edge pixels
- The algorithm has O(MN·R·T) complexity where M×N is image size and R×T is filter size

## References
- Storath, M., & Weinmann, A. (2018). Fast median filtering for phase or orientation data.
  IEEE Transactions on Pattern Analysis and Machine Intelligence, 40(3), 639-652.

See also: [`d`](@ref), [`mirror_index`](@ref)
"""
function arc_distance_median_filter!(u, y, r::Integer, t::Integer)
    ## Input validation
    size(u) == size(y) || error("Input and output arrays must have the same size")
    ndims(u) == 2 || error("Input and output arrays must be 2-dimensional")
    r >= 0 || error("Radius r must be non-negative")
    t >= 0 || error("Radius t must be non-negative")
    eltype(u) <: Real || error("Output array must be of a real type")
    eltype(y) <: Real || error("Input array must be of a real type")

    ## Calculate filter dimensions
    R = 2r + 1
    T = 2t + 1
    N = maximum(size(y))

    ## Allocate buffer arrays with offset indexing for natural -r:r, -t:t access
    Gbuf = [OffsetArray(zeros(Float64, R, T), (-r):r, (-t):t) for _ in 1:N]
    Zbuf = [OffsetArray(zeros(Float64, R, T), (-r):r, (-t):t) for _ in 1:N]
    Gcurrent = OffsetArray(zeros(Float64, R, T), (-r):r, (-t):t)
    Zcurrent = OffsetArray(zeros(Float64, R, T), (-r):r, (-t):t)

    ## Delegate to core computation function
    return _arc_distance_median_filter_core!(u, y, r, t, Gbuf, Zbuf, Gcurrent, Zcurrent)
end

"""
    _arc_distance_median_filter_core!(u, y, r, t, Gbuf, Zbuf, Gcurrent, Zcurrent)

Core computation function for arc distance median filtering.

This function performs the actual filtering computation using pre-allocated buffers.
It implements the iterative algorithm from Storath & Weinmann (2018) without
input validation or buffer allocation overhead.

## Arguments
- `u`: Pre-allocated output array for filtered results (modified in-place)
- `y`: Input 2D array containing circular data (angles in radians)
- `r::Integer`: Half-width of filter in first dimension (rows)
- `t::Integer`: Half-width of filter in second dimension (columns)
- `Gbuf`: Vector of OffsetArrays for storing G matrices for each column
- `Zbuf`: Vector of OffsetArrays for storing Z matrices for each column
- `Gcurrent`: OffsetArray for current G computation
- `Zcurrent`: OffsetArray for current Z computation

## Returns
- `u`: The filtered output array (same object as input `u`)

## Notes
This is an internal function. Use `arc_distance_median_filter!` for the public interface.
The buffers must be properly sized: each should have dimensions (2r+1) × (2t+1) with
offset indexing from (-r):r and (-t):t.

See also: [`arc_distance_median_filter!`](@ref)
"""
function _arc_distance_median_filter_core!(
    u, y, r::Integer, t::Integer, Gbuf, Zbuf, Gcurrent, Zcurrent
)
    M, N = size(y)

    ## Distance function with boundary handling
    dist(m, n, i, j) = d(mirror_index(y, m, n), mirror_index(y, i, j))

    ## Process first column (m = 1)
    m = 1

    ## Calculate first element G[1,1] by definition
    n = 1
    for k in (-r):r
        for l in (-t):t
            Gcurrent[k, l] = sum(
                dist(m + k, n + l, m + i, n + j) for i in (-r):r, j in (-t):t
            )
        end
    end
    Gbuf[n] .= Gcurrent
    ## Find median: pixel that minimizes sum of distances
    min_idx = argmin(Gcurrent)
    u[m, n] = mirror_index(y, m + min_idx[1], n + min_idx[2])

    ## Process rest of first column using iterative updates
    for n in 2:N
        for k in (-r):r
            for l in (-t):(t - 1)
                Gcurrent[k, l] =
                    Gbuf[n - 1][k, l + 1] + sum(
                        (
                            dist(m + k, n + l, m + i, n + t) -
                            dist(m + k, n + l, m + i, n - t - 1)
                        ) for i in (-r):r
                    )
            end
            ## Handle last column separately
            Gcurrent[k, t] = sum(
                dist(m + k, n + t, m + i, n + j) for i in (-r):r, j in (-t):t
            )
        end
        Gbuf[n] .= Gcurrent
        min_idx = argmin(Gcurrent)
        u[m, n] = mirror_index(y, m + min_idx[1], n + min_idx[2])
    end

    ## Process remaining rows (m >= 2)
    for m in 2:M
        ## Initialize auxiliary array Z for first element of row
        n = 1
        for k in (-r):r
            for l in (-t):t
                Zbuf[n][k, l] = sum(
                    (
                        dist(m + k, n + l, m + r, n + j) -
                        dist(m + k, n + l, m - r - 1, n + j)
                    ) for j in (-t):t
                )
            end
        end

        ## Compute first element G[m,1] using previous row
        for k in (-r):(r - 1)
            for l in (-t):t
                Gcurrent[k, l] = Gbuf[n][k + 1, l] + Zbuf[n][k, l]
            end
        end
        ## Handle last row separately
        for l in (-t):t
            Gcurrent[r, l] = sum(
                dist(m + r, n + l, m + i, n + j) for i in (-r):r, j in (-t):t
            )
        end

        Gbuf[n] .= Gcurrent
        min_idx = argmin(Gcurrent)
        u[m, n] = mirror_index(y, m + min_idx[1], n + min_idx[2])

        ## Process remaining elements in current row
        for n in 2:N
            ## Update Z array
            for k in (-r):r
                for l in (-t):(t - 1)
                    Zcurrent[k, l] =
                        Zbuf[n - 1][k, l + 1] + dist(m + k, n + l, m + r, n + t) -
                        dist(m + k, n + l, m + r, n - t - 1) -
                        dist(m + k, n + l, m - r - 1, n + t) +
                        dist(m + k, n + l, m - r - 1, n - t - 1)
                end
                Zcurrent[k, t] = sum(
                    (
                        dist(m + k, n + t, m + r, n + j) -
                        dist(m + k, n + t, m - r - 1, n + j)
                    ) for j in (-t):t
                )
            end
            Zbuf[n - 1] .= Zcurrent

            ## Compute G[m,n] using iterative relations
            for k in (-r):(r - 1)
                for l in (-t):t
                    Gcurrent[k, l] = Gbuf[n][k + 1, l] + Zcurrent[k, l]
                end
            end

            ## Handle last row
            for l in (-t):(t - 1)
                Gcurrent[r, l] =
                    Gbuf[n - 1][r, l + 1] + sum(
                        (
                            dist(m + r, n + l, m + i, n + t) -
                            dist(m + r, n + l, m + i, n - t - 1)
                        ) for i in (-r):r
                    )
            end
            Gcurrent[r, t] = sum(
                dist(m + r, n + t, m + i, n + j) for i in (-r):r, j in (-t):t
            )

            Gbuf[n] .= Gcurrent
            min_idx = argmin(Gcurrent)
            u[m, n] = mirror_index(y, m + min_idx[1], n + min_idx[2])
        end
    end

    return u
end

end
