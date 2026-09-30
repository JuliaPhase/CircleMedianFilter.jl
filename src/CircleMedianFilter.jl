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
function mirror_index(a::AbstractMatrix, i::Int, j::Int)
    ## Get the size of the array
    n, m = size(a)

    return mirror_index(a, i, j, n, m)
end

function mirror_index(a::AbstractMatrix, i::Int, j::Int, n::Int, m::Int)

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
    arc_distance_median_filter!(u, y, r, t)

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
    Tbuf = float(promote_type(eltype(u), eltype(y)))

    ## Two allocations total, matching the C reference layout:
    ## G is a flat (R × T × N) array storing distance sums for every column.
    ## Z is a single (R × T) array updated in-place, eliminating the N-copy Zbuf.
    G = zeros(Tbuf, R, T, N)
    Z = zeros(Tbuf, R, T)

    return _arc_distance_median_filter_core!(u, y, r, t, G, Z)
end

"""
    _arc_distance_median_filter_core!(u, y, r, t, G, Z)

Core computation function for arc distance median filtering.

This function performs the actual filtering computation using pre-allocated buffers.
It implements the iterative algorithm from Storath & Weinmann (2018) without
input validation or buffer allocation overhead.

## Arguments
- `u`: Pre-allocated output array for filtered results (modified in-place)
- `y`: Input 2D array containing circular data (angles in radians)
- `r::Integer`: Half-width of filter in first dimension (rows)
- `t::Integer`: Half-width of filter in second dimension (columns)
- `G`: Pre-allocated `(2r+1) × (2t+1) × N` array for distance sums (one slice per column)
- `Z`: Pre-allocated `(2r+1) × (2t+1)` auxiliary array, updated in-place each column

## Returns
- `u`: The filtered output array (same object as input `u`)

## Notes
This is an internal function. Use `arc_distance_median_filter!` for the public interface.
The layout matches the C reference implementation: G stores all N column slices in a
single contiguous array; Z is a single slice updated in-place (l ascending), so reading
`Z[k, l+1]` always sees the previous column's value before `Z[k, l]` is overwritten.

See also: [`arc_distance_median_filter!`](@ref)
"""
function _arc_distance_median_filter_core!(u, y, r::Integer, t::Integer, G, Z)
    M, N = size(y)
    R = 2r + 1
    T = 2t + 1
    Tbuf = eltype(G)

    ## Distance function with boundary handling
    ## dist(m::Int, n::Int, i::Int, j::Int)::Tbuf =
    ##     Tbuf(d(mirror_index(y, m, n, M, N), mirror_index(y, i, j, M, N)))

    ## Fast arc distance: identical to the C reference.
    ## Precondition: all values in y must be in [0, 2π).
    ## For arbitrary-range input use d() directly instead.
    dist(m::Int, n::Int, i::Int, j::Int)::Tbuf = begin
        aux = abs(mirror_index(y, m, n, M, N) - mirror_index(y, i, j, M, N))
        ifelse(aux > Tbuf(π), Tbuf(2π) - aux, aux)
    end

    ## ── First row (m = 1) ────────────────────────────────────────────────────
    let m = 1

        ## G[:, :, 1] computed from scratch
        let n = 1
            for k in (-r):r
                for l in (-t):t
                    G[k + r + 1, l + t + 1, n] = sum(
                        dist(m + k, n + l, m + i, n + j) for i in (-r):r, j in (-t):t;
                        init=zero(Tbuf),
                    )
                end
            end
            _, idx = findmin(view(G, :, :, n))
            u[m, n] = mirror_index(y, m + idx[1] - r - 1, n + idx[2] - t - 1, M, N)
        end

        ## G[:, :, n] for n ≥ 2: iterative column update
        for n in 2:N
            for k in (-r):r
                for l in (-t):(t - 1)
                    G[k + r + 1, l + t + 1, n] =
                        G[k + r + 1, l + t + 2, n - 1] + sum(
                            dist(m + k, n + l, m + i, n + t) -
                            dist(m + k, n + l, m + i, n - t - 1) for i in (-r):r;
                            init=zero(Tbuf),
                        )
                end
                ## l = t boundary
                G[k + r + 1, T, n] = sum(
                    dist(m + k, n + t, m + i, n + j) for i in (-r):r, j in (-t):t;
                    init=zero(Tbuf),
                )
            end
            _, idx = findmin(view(G, :, :, n))
            u[m, n] = mirror_index(y, m + idx[1] - r - 1, n + idx[2] - t - 1, M, N)
        end
    end

    ## ── Remaining rows (m ≥ 2) ───────────────────────────────────────────────
    for m in 2:M

        let n = 1
            ## Init Z from scratch (all (k,l) computed independently)
            for k in (-r):r
                for l in (-t):t
                    Z[k + r + 1, l + t + 1] = sum(
                        dist(m + k, n + l, m + r, n + j) -
                        dist(m + k, n + l, m - r - 1, n + j) for j in (-t):t;
                        init=zero(Tbuf),
                    )
                end
            end

            ## Update G[:, :, 1] in-place using previous m's G[:, :, 1] and Z.
            ## k ascending: G[k+r+2, :, n] (i.e. k+1) is read before G[k+r+1, :, n] is written.
            for k in (-r):(r - 1)
                for l in (-t):t
                    G[k + r + 1, l + t + 1, n] =
                        G[k + r + 2, l + t + 1, n] + Z[k + r + 1, l + t + 1]
                end
            end
            ## k = r boundary: last row computed from scratch
            for l in (-t):t
                G[R, l + t + 1, n] = sum(
                    dist(m + r, n + l, m + i, n + j) for i in (-r):r, j in (-t):t;
                    init=zero(Tbuf),
                )
            end
            _, idx = findmin(view(G, :, :, n))
            u[m, n] = mirror_index(y, m + idx[1] - r - 1, n + idx[2] - t - 1, M, N)
        end

        for n in 2:N
            ## Update Z in-place: l ascending ensures Z[k, l+t+2] (i.e. l+1)
            ## still holds the previous column's value when Z[k, l+t+1] is written.
            for k in (-r):r
                for l in (-t):(t - 1)
                    Z[k + r + 1, l + t + 1] =
                        Z[k + r + 1, l + t + 2] + dist(m + k, n + l, m + r, n + t) -
                        dist(m + k, n + l, m + r, n - t - 1) -
                        dist(m + k, n + l, m - r - 1, n + t) +
                        dist(m + k, n + l, m - r - 1, n - t - 1)
                end
                ## l = t boundary
                Z[k + r + 1, T] = sum(
                    dist(m + k, n + t, m + r, n + j) - dist(m + k, n + t, m - r - 1, n + j)
                    for j in (-t):t;
                    init=zero(Tbuf),
                )
            end

            ## Update G[:, :, n] in-place using previous m's G[:, :, n] and Z.
            for k in (-r):(r - 1)
                for l in (-t):t
                    G[k + r + 1, l + t + 1, n] =
                        G[k + r + 2, l + t + 1, n] + Z[k + r + 1, l + t + 1]
                end
            end

            ## k = r, l in -t:t-1: use previous column's G (n-1) iteratively
            for l in (-t):(t - 1)
                G[R, l + t + 1, n] =
                    G[R, l + t + 2, n - 1] + sum(
                        dist(m + r, n + l, m + i, n + t) -
                        dist(m + r, n + l, m + i, n - t - 1) for i in (-r):r;
                        init=zero(Tbuf),
                    )
            end
            ## k = r, l = t: full sum from scratch
            G[R, T, n] = sum(
                dist(m + r, n + t, m + i, n + j) for i in (-r):r, j in (-t):t;
                init=zero(Tbuf),
            )

            _, idx = findmin(view(G, :, :, n))
            u[m, n] = mirror_index(y, m + idx[1] - r - 1, n + idx[2] - t - 1, M, N)
        end
    end

    return u
end

include("demo_utils.jl")

end
