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
using PrettyTables

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
## Creating Test Phase Data

We create a complex synthetic phase image with multiple challenging features:
=#

gt_phase = CircleMedianFilter.create_test_phase()

#=
## Filter Performance Analysis

Now we test the filter performance across different noise levels and filter sizes:
=#

## Add random noise to the phase
Random.seed!(42)  ## For reproducible results

## Test different noise levels
const noise_levels = [0.05, 0.25, 0.5, 0.75]
const gaussian_levels = [0.1, 0.03]  ## Test two levels of Gaussian noise
const filter_radii = [1, 2, 3]  ## This creates 3×3, 5×5, and 7×7 filters

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
CircleMedianFilter.add_phase_noise!(noisy_phase, noise_fraction, gaussian_std);

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
## Noise-Reduction Summary Table

Compute noise-reduction statistics for multiple noise settings and filter sizes,
then print a compact results table.
=#

summary_rows = NamedTuple[]

for noise_fraction in noise_levels
    for gaussian_std in gaussian_levels
        noisy = copy(gt_phase)
        CircleMedianFilter.add_phase_noise!(noisy, noise_fraction, gaussian_std)

        for r in filter_radii
            filtered = similar(noisy)
            arc_distance_median_filter!(filtered, noisy, r, r)

            err_noisy = CircleMedianFilter.d.(gt_phase, noisy)
            err_filtered = CircleMedianFilter.d.(gt_phase, filtered)

            mean_noisy = mean(err_noisy)
            mean_filtered = mean(err_filtered)
            improvement = 100 * (mean_noisy - mean_filtered) / mean_noisy

            push!(
                summary_rows,
                (
                    salt_pepper=noise_fraction,
                    gaussian_sigma=gaussian_std,
                    filter_size=2r + 1,
                    mean_error_noisy=mean_noisy,
                    mean_error_filtered=mean_filtered,
                    improvement_pct=improvement,
                ),
            )
        end
    end
end

tablestr = pretty_table(
    summary_rows;
    backend=:html,
    formatters=[fmt__round(4, [4, 5]), fmt__round(2, [1, 2, 6])],
)

Base.HTML(tablestr)



#=
## Runtime Benchmark Across Filter Sizes (Log-Scale)

Measure runtime versus filter size for two image resolutions and visualize the
results on a log-scale y-axis.
=#

function median_runtime_seconds(input_phase, r; repeats=5)
    out = similar(input_phase)

    ## Warmup to avoid first-call compilation costs in timing loop.
    arc_distance_median_filter!(out, input_phase, r, r)

    times = Float64[]
    for _ in 1:repeats
        t = @elapsed arc_distance_median_filter!(out, input_phase, r, r)
        push!(times, t)
    end
    return median(times)
end

const benchmark_radii = 1:15  ## 3x3 up to 31x31
const benchmark_sizes = [256, 512]

benchmark_inputs = Dict{Int,Matrix{Float64}}()
for n in benchmark_sizes
    phase = CircleMedianFilter.create_test_phase(n, n)
    CircleMedianFilter.add_phase_noise!(phase, 0.25, 0.1)
    benchmark_inputs[n] = phase
end

runtime_results = Dict{Int,Vector{Float64}}()
for n in benchmark_sizes
    runtime_results[n] = [
        median_runtime_seconds(benchmark_inputs[n], r) for r in benchmark_radii
    ]
end


xvals = [(2r + 1)^2 for r in benchmark_radii]
runtime_fig = Figure(; size=(900, 500))
runtime_ax = Axis(
    runtime_fig[1, 1];
    title="Runtime vs Number of Filter-Mask Elements",
    xlabel="Number of elements in filter mask",
    ylabel="Runtime [s]",
    xscale=log10,
    yscale=log10,
    xticks=xvals,
)


for n in benchmark_sizes
    scatterlines!(
        runtime_ax, xvals, runtime_results[n]; label="$(n)x$(n)", linewidth=2, markersize=8
    )
end

axislegend(runtime_ax; position=:lt)
runtime_fig

#=
## Key Observations

1. **Edge Preservation**: The circular median filter preserves sharp phase boundaries while reducing noise
2. **Circular Consistency**: Unlike standard median filtering, the results maintain circular data properties
3. **Noise Reduction**: Significant reduction in both salt-and-pepper and Gaussian noise
4. **Parameter Sensitivity**: Larger filters provide more smoothing but may over-blur fine details

This demonstration shows that circular median filtering is essential for proper processing of phase data,
providing robust noise reduction while preserving the fundamental properties of circular measurements.
=#
