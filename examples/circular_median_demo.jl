# ```@meta
# CurrentModule = CircleMedianFilter
# DocTestSetup = quote
#     using CircleMedianFilter
# end
# ```


using CircleMedianFilter
using CairoMakie
using Statistics
using Random

# # Circular Median Filter Demonstration
#
# This example demonstrates the effectiveness of circular median filtering for phase data
# containing various types of noise and discontinuities.
#
# ## Overview
#
# The circular median filter is specifically designed for circular data (like phase images)
# where standard median filtering fails due to the wraparound nature of angular data.
# This demo shows:
#
# 1. **Synthetic Phase Generation**: Creating test phase patterns with circles, stripes, and checkerboards
# 2. **Noise Addition**: Adding salt-and-pepper and Gaussian noise
# 3. **Filter Performance**: Comparing results across different filter sizes and noise levels
# 4. **Quantitative Analysis**: Measuring improvement in terms of phase error

#=
## Helper Functions

First, we define utility functions for creating test patterns and adding noise:
=#

"""
    add_concentric_circles!(phase_array, X, Y, center_x, center_y, radii, phases)

Add concentric circles with different phase values to a phase array.

## Arguments
- `phase_array`: Input phase array to be modified (modified in-place)
- `X`, `Y`: Coordinate grids for the image
- `center_x`, `center_y`: Center coordinates for the circles
- `radii`: Vector of radius pairs [(r1_inner, r1_outer), (r2_inner, r2_outer), ...]
- `phases`: Vector of phase values for each ring

## Returns
- The modified `phase_array` with circles added
"""
function add_concentric_circles!(phase_array, X, Y, center_x, center_y, radii, phases)
    for (i, ((r_inner, r_outer), phase)) in enumerate(zip(radii, phases))
        circle_mask =
            ((X .- center_x) .^ 2 + (Y .- center_y) .^ 2 .>= r_inner^2) .&
            ((X .- center_x) .^ 2 + (Y .- center_y) .^ 2 .<= r_outer^2)
        phase_array[circle_mask] .= phase
    end
    return phase_array
end

"""
    add_phase_noise!(phase_array, noise_fraction::Float64, gaussian_std::Float64=0.1)

Add salt-and-pepper noise and Gaussian noise to a phase array.

## Arguments
- `phase_array`: Input phase array to be corrupted with noise (modified in-place)
- `noise_fraction`: Fraction of pixels to corrupt with salt-and-pepper noise (0.0 to 1.0)
- `gaussian_std`: Standard deviation of Gaussian noise to add to all pixels (default: 0.1)

## Returns
- The modified `phase_array` with noise added, wrapped to [0, 2π] range
"""
function add_phase_noise!(phase_array, noise_fraction::Float64, gaussian_std::Float64=0.1)
    ## Add salt-and-pepper noise (random phase jumps)
    noise_mask = rand(size(phase_array)...) .< noise_fraction
    phase_array[noise_mask] .= rand(sum(noise_mask)) .* 2π

    ## Add Gaussian noise to all pixels
    gaussian_noise = gaussian_std * randn(size(phase_array)...)
    phase_array .+= gaussian_noise

    ## Wrap phase to [0, 2π] range
    phase_array .= mod.(phase_array, 2π)

    return phase_array
end

#=
## Creating Test Phase Data

We create a complex synthetic phase image with multiple challenging features:
=#

## Create a test phase from several concentric circles, wedge, and rectangles and add noise
gt_phase = zeros(512, 512)

## Get image dimensions and center coordinates
ny, nx = size(gt_phase)
cx, cy = nx ÷ 2, ny ÷ 2

## Create coordinate grids
x = 1:nx
y = 1:ny
X = repeat(x', ny, 1)
Y = repeat(y, 1, nx)

#=
### 1. Checkerboard Pattern (Top-Left Quadrant)

First, we add a checkerboard pattern to test the filter's behavior on sharp transitions:
=#

## Add rectangular regions first
## Rectangle 1: top-left quadrant with checkerboard pattern, 8x8 cells
rect_mask1 = (X .<= cx) .& (Y .<= cy)
if any(rect_mask1)
    ## Find the bounding box of the rectangular region
    rect_indices = findall(rect_mask1)
    min_x = minimum(getindex.(Tuple.(rect_indices), 2))
    max_x = maximum(getindex.(Tuple.(rect_indices), 2))
    min_y = minimum(getindex.(Tuple.(rect_indices), 1))
    max_y = maximum(getindex.(Tuple.(rect_indices), 1))

    ## Create 8x8 checkerboard pattern within the rectangle
    cell_width = (max_x - min_x + 1) / 8
    cell_height = (max_y - min_y + 1) / 8

    for (i, idx) in enumerate(rect_indices)
        row, col = Tuple(idx)
        ## Determine which cell this pixel belongs to
        cell_row = min(7, floor(Int, (row - min_y) / cell_height))
        cell_col = min(7, floor(Int, (col - min_x) / cell_width))

        ## Checkerboard pattern: alternate between two phase values
        if (cell_row + cell_col) % 2 == 0
            gt_phase[row, col] = π / 4      ## Light squares
        else
            gt_phase[row, col] = 7π / 4     ## Dark squares
        end
    end
end

#=
### 2. Background Region (Bottom-Right Quadrant)

Set a uniform background for the concentric circles:
=#

## Rectangle 2: bottom-right quadrant, phase = 5π/6 (background for circles)
rect_mask2 = (X .>= cx) .& (Y .>= cy)
gt_phase[rect_mask2] .= 5π / 6

#=
### 3. Concentric Circles

Add concentric circles with increasing phase values to create smooth gradients with sharp boundaries:
=#

## Add concentric circles (larger, centered in bottom-right quadrant) - AFTER rectangles
## Position circles at center of bottom-right quadrant
circle_cx, circle_cy = cx + cx ÷ 2, cy + cy ÷ 2  ## Center of bottom-right quadrant

## Calculate quadrant size and make circles 90% of quadrant radius
quadrant_radius = min(cx ÷ 2, cy ÷ 2)  ## Half of quadrant dimension
max_circle_radius = Int(round(0.9 * quadrant_radius))  ## 90% of quadrant radius

## Create 6 concentric circles with 2×2π total phase change
num_circles = 6
circle_radii = [
    (
        Int(round(i * max_circle_radius / num_circles)),
        Int(round((i + 1) * max_circle_radius / num_circles)),
    ) for i in 0:(num_circles - 1)
]
## Phase values spanning 2×2π = 4π total, evenly distributed
circle_phases = [i * 4π / num_circles for i in 0:(num_circles - 1)]
add_concentric_circles!(gt_phase, X, Y, circle_cx, circle_cy, circle_radii, circle_phases)

#=
### 4. Phase Gradient Stripes

Add horizontal stripes with different phase gradients to test handling of rapid phase changes:
=#

## Add some sharp phase boundaries to test discontinuity handling
## Three horizontal stripes with different phase gradients
## Stripe 1: Moderate gradient (2π total growth)
stripe_mask1 = (Y .>= cy - 40) .& (Y .<= cy - 10) .& (X .>= cx - 100) .& (X .<= cx + 100)
stripe_phase1 = (X .- (cx - 100)) ./ 200 .* 2π  ## Linear gradient from 0 to 2π
gt_phase[stripe_mask1] .= stripe_phase1[stripe_mask1]

## Stripe 2: Higher gradient (3 × 2π total growth)
stripe_mask2 = (Y .>= cy - 10) .& (Y .<= cy + 10) .& (X .>= cx - 100) .& (X .<= cx + 100)
stripe_phase2 = (X .- (cx - 100)) ./ 200 .* 6π  ## Linear gradient from 0 to 6π
gt_phase[stripe_mask2] .= stripe_phase2[stripe_mask2]

## Stripe 3: Very high gradient (10 × 2π total growth) - original stripe
stripe_mask3 = (Y .>= cy + 10) .& (Y .<= cy + 40) .& (X .>= cx - 100) .& (X .<= cx + 100)
stripe_phase3 = (X .- (cx - 100)) ./ 200 .* 20π  ## Linear gradient from 0 to 20π
gt_phase[stripe_mask3] .= stripe_phase3[stripe_mask3]

## Wrap all phases to [0, 2π] to create realistic phase discontinuities
gt_phase .= mod.(gt_phase, 2π)

#=
## Filter Performance Analysis

Now we test the filter performance across different noise levels and filter sizes:
=#

## Add random noise to the phase
Random.seed!(42)  ## For reproducible results

## Test different noise levels
noise_levels = [0.05, 0.25, 0.5, 0.75]
gaussian_levels = [0.1, 0.03]  ## Test two levels of Gaussian noise
filter_radii = [1, 2, 3]  ## This creates 3×3, 5×5, and 7×7 filters

println("Testing circular median filter performance...")
println("Noise levels: ", noise_levels)
println("Gaussian levels: ", gaussian_levels)
println("Filter sizes: ", [2r + 1 for r in filter_radii])

#=
### Example: Single Demonstration

Let's show a detailed example with moderate noise:
=#

## Example with moderate noise level
noise_fraction = 0.25
gaussian_std = 0.1
filter_radius = 2 # 5×5 filter

## Create noisy version
noisy_phase = copy(gt_phase)
add_phase_noise!(noisy_phase, noise_fraction, gaussian_std)

## Apply circular median filter
filtered_phase = similar(noisy_phase)
arc_distance_median_filter!(filtered_phase, noisy_phase, filter_radius, filter_radius)

#=
### Visualization

Create comparison plots to visualize the filtering results:
=#

CairoMakie.activate!(; type="png")

## Create comparison plots
fig = Figure(; size=(1200, 1000));

## Ground truth phase
ax1 = Axis(fig[1, 1]; title="Ground Truth Phase", aspect=DataAspect())
hm1 = heatmap!(ax1, gt_phase; colormap=:hsv, colorrange=(0, 2π))
Colorbar(fig[1, 2], hm1; label="Phase [rad]")

## Noisy phase
ax2 = Axis(fig[1, 3]; title="Noisy Phase", aspect=DataAspect())
hm2 = heatmap!(ax2, noisy_phase; colormap=:hsv, colorrange=(0, 2π))
Colorbar(fig[1, 4], hm2; label="Phase [rad]")

## Filtered phase
filter_size = 2 * filter_radius + 1
ax3 = Axis(
    fig[2, 1]; title="Filtered Phase ($(filter_size)×$(filter_size))", aspect=DataAspect()
)
hm3 = heatmap!(ax3, filtered_phase; colormap=:hsv, colorrange=(0, 2π))
Colorbar(fig[2, 2], hm3; label="Phase [rad]")

## Error maps
error_noisy = CircleMedianFilter.d.(gt_phase, noisy_phase)
error_filtered = CircleMedianFilter.d.(gt_phase, filtered_phase)

## Error: Filtered vs GT
ax4 = Axis(fig[2, 3]; title="Error: Filtered vs GT", aspect=DataAspect())
hm4 = heatmap!(ax4, error_filtered; colormap=:hot, colorrange=(0, π))
Colorbar(fig[2, 4], hm4; label="Error [rad]")

Label(
    fig[0, :],
    "Circular Median Filter Demo\nSalt-pepper: $(noise_fraction), Gaussian σ: $(gaussian_std), Filter: $(filter_size)×$(filter_size)";
    fontsize=16,
    halign=:center,
)

## Display the figure
fig

#=
### Quantitative Results

Calculate error statistics to quantify the improvement:
=#

## Statistics
mean_error_noisy = mean(error_noisy)
mean_error_filtered = mean(error_filtered)

println(
    "Results for salt-pepper: $(noise_fraction), Gaussian σ: $(gaussian_std), Filter $(filter_size)×$(filter_size):",
)
println("  Mean error (noisy): $(round(mean_error_noisy, digits=4)) radians")
println("  Mean error (filtered): $(round(mean_error_filtered, digits=4)) radians")
println(
    "  Improvement: $(round((mean_error_noisy - mean_error_filtered) / mean_error_noisy * 100, digits=2))%",
)

#=
## Key Observations

1. **Edge Preservation**: The circular median filter preserves sharp phase boundaries while reducing noise
2. **Circular Consistency**: Unlike standard median filtering, the results maintain circular data properties
3. **Noise Reduction**: Significant reduction in both salt-and-pepper and Gaussian noise
4. **Parameter Sensitivity**: Larger filters provide more smoothing but may over-blur fine details

This demonstration shows that circular median filtering is essential for proper processing of phase data,
providing robust noise reduction while preserving the fundamental properties of circular measurements.
=#
