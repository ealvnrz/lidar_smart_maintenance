# Workflow Documentation

## Complete Workflow Description

This document describes the complete workflow for analyzing mining haul road LiDAR data using the Weibull Random Field framework.

## Step-by-Step Workflow

### Step 1: Data Processing and Segmentation

**Purpose**: Load LiDAR point cloud data and partition the road into homogeneous spatial segments.

**Process**:
1. Load LAS file using `lidR::readLAS()`
2. Extract coordinates (x, y, z), RGB values, and intensity
3. Normalize all values to [0, 1] range for clustering
4. Perform K-means clustering on normalized spatial coordinates (x, y, z)
5. Split data by cluster assignment
6. Save each cluster as a separate CSV file

**Functions**: `load_las_data()`, `segment_road()`, `save_sections()`

**Output**: `results/sections/LAS_sec_*.csv` (250 files, one per cluster)

**Parameters**:
- `NUM_CLUSTERS`: 250 (default)
- `N_START`: 25 (number of random starts for K-means)
- `CLUSTERING_SEED`: 234528 (for reproducibility)

### Step 2: Marginal Analysis (Diagnostic)

**Purpose**: Validate that intensity data follows Weibull distribution at the marginal level.

**Process**:
1. For each cluster (or sample clusters), fit Weibull distribution ignoring spatial correlation
2. Generate diagnostic plots:
   - Histogram with Weibull overlay
   - QQ plot comparing theoretical vs empirical quantiles
3. Assess goodness of fit

**Functions**: `fit_weibull_marginal()`, `diagnose_marginals()`, `check_marginal_fit()`

**Output**: Diagnostic plots (saved to `figures/results/`)

**Note**: This step is optional but recommended for validation.

### Step 3: Spatial Parameter Estimation

**Purpose**: Estimate Weibull Random Field parameters accounting for spatial dependence.

**Process**:
1. For each cluster:
   - Load section data
   - Extract coordinates and intensity
   - Fit WRF using pairwise likelihood with Wendland correlation
   - Estimate parameters: mean, scale, shape, correlation range
2. Process all clusters in parallel
3. Save parameter estimates to CSV

**Functions**: `estimate_wrf_parameters()`, `estimate_all_sections()`

**Output**: `results/estimation/param_estimation.csv`

**Parameters**:
- Correlation model: `Wend0` (Wendland with compact support)
- Optimization: `nlminb` with bounds
- Neighbors: 10 (for pairwise likelihood)

**Model**:
- Marginal distribution: Weibull
- Spatial correlation: Wendland (compactly supported)
- Estimation: Pairwise likelihood (weighted composite likelihood)

### Step 4: Health Indicator Computation

**Purpose**: Compute operational health indicators for maintenance decision-making.

**Process**:
1. Compute global reference intensity I₀ (99th percentile across all clusters)
2. For each cluster, compute:
   - **ELI** (Energy-Loss Index): Average logarithmic deviation from reference
   - **TMR** (Tail-to-Mean Ratio): Ratio of 95th percentile to mean
   - **EA** (Entropy of Attenuation): Entropy of intensity distribution
3. Flag degraded clusters:
   - Indicator > 75th percentile AND
   - Mean intensity < 25th percentile

**Functions**: `compute_eli()`, `compute_tmr()`, `compute_ea()`, `compute_all_indicators()`, `flag_degraded_clusters()`

**Output**:
- `results/health_indicators/health_indicators.csv`: All indicators
- `results/health_indicators/flagged/flagged_ELI.csv`: ELI-flagged clusters
- `results/health_indicators/flagged/flagged_TMR.csv`: TMR-flagged clusters
- `results/health_indicators/flagged/flagged_EA.csv`: EA-flagged clusters

**Formulas**:
- ELI_k = -(1/n_k) * Σ ln(I_i / I₀)
- TMR_{0.95,k} = q_{0.95,k} / mean(I_k)
- EA_k = -Σ p_{j,k} * ln(p_{j,k}) where p_j are histogram proportions

### Step 5: Visualization (Optional)

**Purpose**: Generate diagnostic and summary plots.

**Functions**: `plot_intensity_histogram()`, `plot_parameter_boxplots()`, `plot_flagged_clusters()`, etc.

**Output**: Figures saved to `figures/results/`

## Data Flow

```
LAS File
  ↓
[Load & Normalize]
  ↓
[K-means Clustering]
  ↓
Section CSVs (250 files)
  ↓
[WRF Parameter Estimation]
  ↓
Parameter Estimates CSV
  ↓
[Health Indicator Computation]
  ↓
Indicators CSV + Flagged Clusters
```

## Parallel Processing

The workflow uses parallel processing for:
- **Parameter estimation**: Processes multiple clusters simultaneously
- **Health indicator computation**: Computes indicators for multiple clusters in parallel

Number of cores: `detectCores() - 1` (configurable in `config.R`)

## Error Handling

- Failed cluster fits are logged but processing continues
- Missing sections are reported but don't stop the workflow
- Invalid parameters return NA values with warnings

## Reproducibility

- Random seeds are set in `config.R`
- All random operations use fixed seeds
- Session info is printed at the end of `run_analysis.R`

## Performance Considerations

- **Segmentation**: Single-threaded (fast, ~minutes)
- **Parameter estimation**: Parallel (slowest step, ~hours for 250 clusters)
- **Health indicators**: Parallel (moderate, ~minutes)
- **Visualization**: Single-threaded (fast, ~seconds)

## Troubleshooting

### Common Issues

1. **Out of memory**: Reduce number of clusters or process in batches
2. **Convergence failures**: Adjust starting parameters in `config.R`
3. **Missing sections**: Check that all section files exist in `results/sections/`
4. **Package errors**: Ensure all required packages are installed

### Validation Checks

- Check that all 250 section files exist
- Verify parameter estimation file has 250 rows
- Confirm health indicators file has all clusters
- Review diagnostic plots for model fit quality

## Extending the Workflow

### Adding New Indicators

1. Create function in `R/04_health_indicators.R`
2. Add to `compute_cluster_indicators()`
3. Update `compute_all_indicators()` to include new indicator

### Custom Correlation Models

1. Modify `CORRMODEL` in `config.R`
2. Update `estimate_wrf_parameters()` if needed
3. Ensure GeoModels supports the correlation function

### Different Segmentation Methods

1. Replace `segment_road()` in `R/01_data_processing.R`
2. Maintain output format (CSV files with Cluster column)
3. Update `NUM_CLUSTERS` in config

