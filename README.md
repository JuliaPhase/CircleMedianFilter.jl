# CircleMedianFilter.jl

Part of the [Phase.jl](https://github.com/JuliaPhase/Phase.jl) ecosystem.

<!-- DOI badge: add after first Zenodo release -->

[![Build Status](https://github.com/JuliaPhase/CircleMedianFilter.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPhase/CircleMedianFilter.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaphase.github.io/CircleMedianFilter.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaphase.github.io/CircleMedianFilter.jl/dev/)

A Julia implementation of fast median filtering algorithms for circular data (phase images, orientation fields, angular measurements).

## Overview

CircleMedianFilter.jl provides efficient implementations of the arc distance median filter, specifically designed for data that lives on the unit circle. This is essential for processing phase data, where traditional median filtering fails due to the wraparound nature of angular measurements.

## Features

- **Arc Distance Median**: Robust filtering that respects circular geometry
- **Edge-Preserving**: Maintains sharp boundaries in phase images
- **High Performance**: Optimized implementation for large datasets
- **Memory Efficient**: In-place operations to minimize memory usage

## Installation

```julia
using Pkg
Pkg.add("CircleMedianFilter")
```

## Quick Start

```julia
using CircleMedianFilter

# Your phase data in radians [0, 2π]
phase_data = your_phase_array

# Apply 5×5 circular median filter
filtered_data = similar(phase_data)
arc_distance_median_filter!(filtered_data, phase_data, 2, 2)  # radius_x=2, radius_y=2
```

## Applications

- **Interferometric SAR**: Phase image denoising and preprocessing for unwrapping
- **Computer Vision**: Processing orientation fields in texture analysis
- **Optical Flow**: Filtering directional flow fields
- **Meteorology**: Processing wind direction time series
- **Signal Processing**: Any circular/angular data filtering

## Algorithm Reference

This package implements Algorithm 1 from:

**"Fast Median Filtering for Phase or Orientation Data"**
by Martin Storath and Andreas Weinmann
*IEEE Transactions on Pattern Analysis and Machine Intelligence*, Vol. 40, No. 3, March 2018
DOI: [10.1109/TPAMI.2017.2692779](https://doi.org/10.1109/TPAMI.2017.2692779)

## Documentation

For detailed documentation, examples, and API reference, visit:
- [**Stable Documentation**](https://juliaphase.github.io/CircleMedianFilter.jl/stable/)
- [**Development Documentation**](https://juliaphase.github.io/CircleMedianFilter.jl/dev/)

## Funding

This work is part of the 14AMI project (grant agreement No 101111948). The project is supported by the Chips Joint Undertaking and its members including the top-up funding by RVO (The Netherlands Enterprise Agency).

<img src="docs/src/assets/funding/Chips-JU.png" alt="Chips Joint Undertaking, co-funded by the European Union" height="60">
