# Developer utilities for creating synthetic phase test data.
# These functions are not part of the public API and are intended for
# use in examples, benchmarks, and documentation.

"""
    add_concentric_circles!(phase_array, X, Y, center_x, center_y, radii, phases)

Add concentric rings with specified phase values to a phase array (in-place).

`radii` is a vector of `(r_inner, r_outer)` pairs; `phases` is the corresponding
phase value for each ring. Intended for generating synthetic test data.
"""
function add_concentric_circles!(phase_array, X, Y, center_x, center_y, radii, phases)
    for ((r_inner, r_outer), phase) in zip(radii, phases)
        circle_mask =
            ((X .- center_x) .^ 2 + (Y .- center_y) .^ 2 .>= r_inner^2) .&
            ((X .- center_x) .^ 2 + (Y .- center_y) .^ 2 .<= r_outer^2)
        phase_array[circle_mask] .= phase
    end
    return phase_array
end

"""
    add_phase_noise!(phase_array, noise_fraction, gaussian_std=0.1)

Corrupt a phase array with salt-and-pepper noise (`noise_fraction` of pixels replaced
by uniform random values in [0, 2π]) and additive Gaussian noise (std `gaussian_std`).
Values are wrapped back to [0, 2π] afterwards. Modifies `phase_array` in-place.
"""
function add_phase_noise!(phase_array, noise_fraction::Float64, gaussian_std::Float64=0.1)
    noise_mask = rand(size(phase_array)...) .< noise_fraction
    phase_array[noise_mask] .= rand(sum(noise_mask)) .* 2π
    phase_array .+= gaussian_std .* randn(size(phase_array)...)
    phase_array .= mod.(phase_array, 2π)
    return phase_array
end

"""
    create_test_phase(ny=512, nx=512)

Return a synthetic wrapped phase map (`ny × nx`) suitable for benchmarking and
demonstrating the circular median filter. The image contains:

- Top-left quadrant: 8×8 checkerboard alternating between π/4 and 7π/4.
- Bottom-right quadrant: uniform background at 5π/6, overlaid with 6 concentric
  rings whose phases span [0, 4π] (wrapped to [0, 2π]).
- Central horizontal band: three gradient stripes with total phase growth of
  2π, 6π, and 20π respectively, to test behaviour under steep gradients.
"""
function create_test_phase(ny::Integer=512, nx::Integer=512)
    gt_phase = zeros(Float64, ny, nx)

    cx, cy = nx ÷ 2, ny ÷ 2
    X = repeat((1:nx)', ny, 1)
    Y = repeat(1:ny, 1, nx)

    ## Top-left checkerboard (8×8 cells)
    rect_mask1 = (X .<= cx) .& (Y .<= cy)
    if any(rect_mask1)
        rect_indices = findall(rect_mask1)
        min_x = minimum(getindex.(Tuple.(rect_indices), 2))
        max_x = maximum(getindex.(Tuple.(rect_indices), 2))
        min_y = minimum(getindex.(Tuple.(rect_indices), 1))
        max_y = maximum(getindex.(Tuple.(rect_indices), 1))
        cell_width = (max_x - min_x + 1) / 8
        cell_height = (max_y - min_y + 1) / 8
        for idx in rect_indices
            row, col = Tuple(idx)
            cell_row = min(7, floor(Int, (row - min_y) / cell_height))
            cell_col = min(7, floor(Int, (col - min_x) / cell_width))
            gt_phase[row, col] = (cell_row + cell_col) % 2 == 0 ? π / 4 : 7π / 4
        end
    end

    ## Bottom-right uniform background for the concentric circles
    gt_phase[(X .>= cx) .& (Y .>= cy)] .= 5π / 6

    ## Concentric rings centred in the bottom-right quadrant
    circle_cx, circle_cy = cx + cx ÷ 2, cy + cy ÷ 2
    max_circle_radius = Int(round(0.9 * min(cx ÷ 2, cy ÷ 2)))
    num_circles = 6
    circle_radii = [
        (
            Int(round(i * max_circle_radius / num_circles)),
            Int(round((i + 1) * max_circle_radius / num_circles)),
        ) for i in 0:(num_circles - 1)
    ]
    circle_phases = [i * 4π / num_circles for i in 0:(num_circles - 1)]
    add_concentric_circles!(
        gt_phase, X, Y, circle_cx, circle_cy, circle_radii, circle_phases
    )

    ## Three horizontal gradient stripes (increasing steepness)
    for (dy_lo, dy_hi, total_phase) in ((-40, -10, 2π), (-10, 10, 6π), (10, 40, 20π))
        mask =
            (Y .>= cy + dy_lo) .& (Y .<= cy + dy_hi) .& (X .>= cx - 100) .& (X .<= cx + 100)
        stripe_phase = (X .- (cx - 100)) ./ 200 .* total_phase
        gt_phase[mask] .= stripe_phase[mask]
    end

    gt_phase .= mod.(gt_phase, 2π)
    return gt_phase
end
