# About CircleMedianFilter.jl

## Background

CircleMedianFilter.jl addresses a fundamental challenge in signal and image processing: how to apply robust filtering techniques to data that lives on a circle. Traditional median filtering, while excellent for linear data, fails when applied directly to circular data because it doesn't account for the wraparound nature of angular measurements.

## The Problem with Standard Median Filtering

Consider phase values near the wraparound point (e.g., 0.3, 6.2, 0.05 radians). A standard median filter might return a value around 0.3 radians, which is completely wrong, as the correct median should be close to 0.05 radians. This happens because standard algorithms don't understand that 6.2 radians is actually something about -0.08 radians on the unit circle.

## The Solution: Arc Distance Median

The arc distance median solves this problem by:

1. **Understanding Circular Geometry**: Using the shortest arc distance between angles instead of linear differences
2. **Preserving Circular Structure**: Ensuring that the median of circular data remains meaningful in the circular domain
3. **Maintaining Robustness**: Providing the same outlier resistance as traditional median filtering

## Mathematical Foundation

The package implements the arc distance median as defined by Storath and Weinmann. For a set of angles θ₁, θ₂, ..., θₙ, the arc distance median θ* minimizes:

```
∑ᵢ d(θᵢ, θ*)
```

where d(a, b) is the shortest angular distance between angles a and b on the unit circle.

## Research Applications

This filtering approach is particularly valuable in:

- **Remote Sensing**: Processing interferometric SAR data where phase continuity is crucial
- **Computer Vision**: Handling orientation fields in texture and flow analysis
- **Meteorology**: Filtering directional measurements like wind direction
- **Robotics**: Processing angular sensor data and orientation measurements

## Performance Characteristics

The implementation focuses on:

- **Computational Efficiency**: Optimized algorithms for large-scale data processing
- **Memory Efficiency**: In-place operations to minimize memory footprint
- **Numerical Stability**: Robust handling of edge cases and numerical precision issues

## Future Extensions

The package is designed to be extensible for:

- Different kernel shapes and sizes
- Alternative circular distance metrics
- Multi-dimensional circular data
- Adaptive filtering approaches
