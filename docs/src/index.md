# CircleMedianFilter.jl

Documentation for CircleMedianFilter.jl - A Julia implementation of fast median filtering for circular/phase data.

## Overview

CircleMedianFilter.jl provides efficient implementations of median filtering algorithms specifically designed for circular data such as phase images, orientation fields, and angular measurements. The package implements Algorithm 1 from the research paper:

**"Fast Median Filtering for Phase or Orientation Data"** by Martin Storath and Andreas Weinmann
*IEEE Transactions on Pattern Analysis and Machine Intelligence*, Vol. 40, No. 3, March 2018
DOI: [10.1109/TPAMI.2017.2692779](https://doi.org/10.1109/TPAMI.2017.2692779); see the original Matlab implementation on publication [github repo](https://github.com/mstorath/CircleMedianFilter).

## Key Features

- **Arc Distance Median**: Robust filtering that preserves the circular nature of data
- **Edge-Preserving**: Maintains sharp boundaries in phase images
- **Value-Preserving**: Preserves original data values where appropriate
- **High Performance**: Optimized implementation for large datasets
- **Memory Efficient**: In-place operations to minimize memory usage

## Applications

This package is particularly useful for:

- **Interferometric SAR**: Phase image denoising and unwrapping preprocessing
- **Optical Flow**: Processing planar flow fields with circular orientation data
- **Wind Direction Analysis**: Filtering time series of wind direction measurements
- **Computer Vision**: Processing orientation fields in texture analysis
- **Signal Processing**: Any application involving data on the unit circle

## Quick Start

```julia
using CircleMedianFilter

# Create some noisy phase data
phase_data = your_phase_array  # Should be in radians [0, 2π]

# Apply circular median filter with 3×3 kernel
filtered_data = similar(phase_data)
arc_distance_median_filter!(filtered_data, phase_data, 1, 1)  # radius_x=1, radius_y=1
```

## Index
```@index
```

## Funding

This work is part of the 14AMI project (grant agreement No 101111948). The project is supported by the Chips Joint Undertaking and its members including the top-up funding by RVO (The Netherlands Enterprise Agency).

```@raw html
<img src="assets/funding/Chips-JU.png" alt="Chips Joint Undertaking, co-funded by the European Union" height="60" style="height: 60px; width: auto; margin-right: 1em;">
```
