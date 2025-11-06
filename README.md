# Weibull Random Field Framework for Mining Haul Road Analysis

This repository contains the supplementary code for the research article:

> **From LiDAR Intensity to Smart Maintenance: A Weibull Random-Field Framework for Diagnosis, Prognosis, and Decision Automation on Mining Haul Roads**

## Overview

This project implements a Weibull Random Field (WRF) framework for analyzing LiDAR intensity data from mining haul roads. The framework enables:

1. **Spatial segmentation** of road surfaces into homogeneous clusters
2. **Weibull marginal fitting** and diagnostic validation
3. **Spatial parameter estimation** using pairwise likelihood with Wendland correlation
4. **Health indicator computation** (ELI, TMR, EA) for degradation assessment
5. **Automated flagging** of degraded road segments

## Project Structure

```
2025_weibull_rf_roads/
├── README.md                 # This file
├── config.R                  # Configuration parameters
├── run_analysis.R            # Main workflow script
│
├── R/                        # Modular functions
│   ├── 01_data_processing.R
│   ├── 02_marginal_analysis.R
│   ├── 03_spatial_estimation.R
│   ├── 04_health_indicators.R
│   ├── 05_visualization.R
│   └── utils.R
│
├── data/
│   └── raw/                  # Input LAS file
│
├── results/                  # Analysis outputs
│   ├── sections/             # Segmented cluster CSVs
│   ├── estimation/            # Parameter estimates
│   ├── health_indicators/     # Health indicator CSVs
│   └── flagged_clusters/      # Flagged degraded clusters
│
├── figures/                  # Generated plots
│   └── results/
│
└── docs/                     # Documentation
    ├── workflow.md
    └── functions.md
```

## Data Availability

The LAS file (`minasur.las`, ~67 MB) used in this research article is not included in this repository due to its size. The data file is available **upon request** from the authors. Please contact the corresponding authors for access to the LiDAR data.

Note: The workflow can be executed from Step 2 onwards if you have pre-processed section files in `results/sections/`.

## Installation

### Required R Packages

```r
install.packages(c(
  "lidR",           # LiDAR data processing
  "GeoModels",      # Spatial modeling
  "readr",          # Data I/O
  "parallel",       # Parallel processing
  "ggplot2",       # Plotting
  "plotly",         # Interactive plots
  "MASS",           # Statistical functions
  "dplyr",          # Data manipulation
  "tidyr",          # Data reshaping
  "fields"          # Spatial statistics
))
```

### Installation from GitHub

```r
# If publishing as a package, install with:
# devtools::install_github("username/repo-name")
```

## Quick Start

### 1. Configure Paths

Edit `config.R` to set:
- Path to your LAS file
- Output directories
- Model parameters (number of clusters, optimization settings)

### 2. Run Analysis

```bash
# Run complete workflow
Rscript run_analysis.R

# Skip steps if outputs already exist
Rscript run_analysis.R --skip-segmentation --skip-estimation
```

Or in R:

```r
source("run_analysis.R")
```

### 3. Workflow Steps

The main script (`run_analysis.R`) performs:

1. **Data Processing**: Loads LAS file and segments road into clusters
2. **Marginal Analysis**: Validates Weibull marginal fits
3. **Spatial Estimation**: Estimates WRF parameters via pairwise likelihood
4. **Health Indicators**: Computes ELI, TMR, and EA indicators
5. **Visualization**: Generates diagnostic plots

## Usage Examples

### Load and Segment Data

```r
source("R/01_data_processing.R")
data <- load_las_data("data/raw/minasur.las")
data <- segment_road(data, num_clusters = 250, n_start = 25)
save_sections(data, "results/sections")
```

### Estimate Parameters

```r
source("R/03_spatial_estimation.R")
results <- estimate_all_sections(
  sections_dir = "results/sections",
  output_file = "results/estimation/param_estimation.csv",
  parallel = TRUE
)
```

### Compute Health Indicators

```r
source("R/04_health_indicators.R")
health_results <- compute_all_indicators(
  sections_dir = "results/sections",
  param_estimation_file = "results/estimation/param_estimation.csv",
  output_dir = "results/health_indicators"
)
```

## Health Indicators

The framework computes three health indicators:

1. **Energy-Loss Index (ELI)**: `ELI_k = -(1/n_k) * sum(ln(I_i/I_0))`
2. **Tail-to-Mean Ratio (TMR)**: `TMR_{0.95,k} = q_{0.95,k} / mean(I_k)`
3. **Entropy of Attenuation (EA)**: `EA_k = -sum(p_{j,k} * ln(p_{j,k}))`

Clusters are flagged as degraded if:
- Indicator value > 75th percentile across clusters
- Mean intensity < 25th percentile across clusters

## Configuration

Key parameters in `config.R`:

- `NUM_CLUSTERS`: Number of road segments (default: 250)
- `CORRMODEL`: Correlation function (default: "Wend0" for Wendland)
- `OPTIMIZER`: Optimization method (default: "nlminb")
- `N_CORES`: Number of parallel cores

See `config.R` for full list of parameters.

## Output Files

- `results/sections/LAS_sec_*.csv`: Segmented cluster data
- `results/estimation/param_estimation.csv`: WRF parameter estimates
- `results/health_indicators/health_indicators.csv`: All health indicators
- `results/flagged_clusters/flagged_*.csv`: Flagged clusters per indicator

