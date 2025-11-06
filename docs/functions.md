# Function Reference

This document provides detailed documentation for all functions in the workflow.

## Data Processing (`R/01_data_processing.R`)

### `load_las_data(file_path)`

Loads and normalizes LAS point cloud data.

**Parameters**:
- `file_path` (character): Path to LAS file

**Returns**: Data frame with original and normalized values (x, y, z, intensity, RGB)

**Example**:
```r
data <- load_las_data("data/raw/minasur.las")
```

---

### `segment_road(data, num_clusters, n_start, seed)`

Performs K-means clustering to segment road into homogeneous clusters.

**Parameters**:
- `data` (data.frame): Data with normalized coordinates
- `num_clusters` (integer): Number of clusters (default: 250)
- `n_start` (integer): Random starts for K-means (default: 25)
- `seed` (integer): Random seed (default: 234528)

**Returns**: Data frame with added `Cluster` column

**Example**:
```r
data_clustered <- segment_road(data, num_clusters = 250, n_start = 25)
```

---

### `save_sections(data, output_dir)`

Saves segmented data to CSV files (one per cluster).

**Parameters**:
- `data` (data.frame): Data with Cluster column
- `output_dir` (character): Directory to save files

**Returns**: Invisible NULL

---

## Marginal Analysis (`R/02_marginal_analysis.R`)

### `fit_weibull_marginal(intensity_data)`

Fits Weibull distribution to intensity data (ignoring spatial correlation).

**Parameters**:
- `intensity_data` (numeric): Vector of intensity values

**Returns**: List with `shape`, `scale`, `fit_object`, `n_obs`

**Example**:
```r
fit <- fit_weibull_marginal(data$intensity_scaled)
```

---

### `diagnose_marginals(fit, data, title)`

Generates diagnostic plots for Weibull fit.

**Parameters**:
- `fit`: Fitted Weibull object
- `data` (numeric): Intensity data
- `title` (character): Plot title

**Returns**: List with `histogram` and `qqplot` ggplot objects

---

### `check_marginal_fit(sections_dir, cluster_id, plot)`

Validates Weibull marginal fit for a specific cluster.

**Parameters**:
- `sections_dir` (character): Directory containing sections
- `cluster_id` (integer): Cluster ID to check
- `plot` (logical): Whether to generate plots

**Returns**: List with fit results and optionally plots

---

## Spatial Estimation (`R/03_spatial_estimation.R`)

### `estimate_wrf_parameters(section_data, coords, start_params, ...)`

Estimates WRF parameters for a single section using pairwise likelihood.

**Parameters**:
- `section_data` (data.frame): Section data with x, y, intensity_scaled
- `coords` (matrix): Optional coordinates (n x 2)
- `start_params` (list): Starting parameters (mean, scale, shape)
- `fixed_params` (list): Fixed parameters (nugget, power2)
- `lower_bounds`, `upper_bounds` (list): Optimization bounds
- `corrmodel` (character): Correlation model (default: "Wend0")
- `optimizer` (character): Optimizer (default: "nlminb")
- `neighb` (integer): Number of neighbors (default: 10)

**Returns**: List with estimated parameters and fit statistics

**Example**:
```r
fit <- estimate_wrf_parameters(section_data,
                               start_params = list(mean = 0.5, scale = 10, shape = 5))
```

---

### `estimate_all_sections(sections_dir, output_file, parallel, ...)`

Estimates WRF parameters for all sections in parallel.

**Parameters**:
- `sections_dir` (character): Directory with section files
- `output_file` (character): Path to output CSV
- `parallel` (logical): Use parallel processing (default: TRUE)
- `n_cores` (integer): Number of cores (default: auto-detect)
- `max_files` (integer): Maximum files to process (default: 999999)

**Returns**: Data frame with parameter estimates

---

## Health Indicators (`R/04_health_indicators.R`)

### `compute_eli(intensity, I0)`

Computes Energy-Loss Index.

**Parameters**:
- `intensity` (numeric): Intensity values
- `I0` (numeric): Global reference intensity (99th percentile)

**Returns**: Numeric ELI value

**Formula**: ELI_k = -(1/n_k) * Σ ln(I_i / I₀)

---

### `compute_tmr(intensity, quantile)`

Computes Tail-to-Mean Ratio.

**Parameters**:
- `intensity` (numeric): Intensity values
- `quantile` (numeric): Quantile level (default: 0.95)

**Returns**: Numeric TMR value

**Formula**: TMR_{q,k} = q_{q,k} / mean(I_k)

---

### `compute_ea(intensity, bins)`

Computes Entropy of Attenuation.

**Parameters**:
- `intensity` (numeric): Intensity values
- `bins` (integer): Number of histogram bins (default: 30)

**Returns**: Numeric EA value

**Formula**: EA_k = -Σ p_{j,k} * ln(p_{j,k})

---

### `compute_all_indicators(sections_dir, param_estimation_file, output_dir, ...)`

Main function to compute all health indicators and flag degraded clusters.

**Parameters**:
- `sections_dir` (character): Directory with section files
- `param_estimation_file` (character): Path to parameter estimates
- `output_dir` (character): Directory to save results
- `I0_quantile` (numeric): Quantile for I0 (default: 0.99)
- `tmr_quantile` (numeric): Quantile for TMR (default: 0.95)
- `ea_bins` (integer): Bins for EA (default: 30)
- `flag_intensity_percentile` (numeric): Flag threshold (default: 0.75)
- `flag_mean_percentile` (numeric): Mean threshold (default: 0.25)
- `parallel` (logical): Use parallel processing (default: TRUE)

**Returns**: List with indicators and flagged clusters data frames

---

## Visualization (`R/05_visualization.R`)

### `plot_intensity_histogram(data, fit, title, ...)`

Creates histogram with Weibull overlay.

**Parameters**:
- `data` (numeric): Intensity values
- `fit`: Weibull fit object
- `title` (character): Plot title
- `bins` (integer): Number of bins (default: 15)

**Returns**: ggplot object

---

### `plot_3d_intensity(coords, intensity, title, ...)`

Creates interactive 3D intensity plot.

**Parameters**:
- `coords`: Data frame or matrix with x, y, z
- `intensity` (numeric): Intensity values
- `title` (character): Plot title
- `sample_size` (numeric): Proportion to sample (default: 0.2)

**Returns**: plotly object

---

### `plot_parameter_boxplots(param_df, param_name, ...)`

Creates boxplot for a parameter across clusters.

**Parameters**:
- `param_df` (data.frame): Parameter estimates
- `param_name` (character): Parameter name ("mean", "scale", "shape")
- `title` (character): Optional title

**Returns**: ggplot object

---

### `plot_flagged_clusters(cluster_data, flagged_ids, indicator_name, ...)`

Creates spatial map of flagged clusters.

**Parameters**:
- `cluster_data` (data.frame): Data with x, y, Cluster columns
- `flagged_ids` (numeric): Vector of flagged cluster IDs
- `indicator_name` (character): Indicator name for title

**Returns**: ggplot object

---

## Utilities (`R/utils.R`)

### `extract_cluster_id(file_path)`

Extracts cluster ID from section file path.

**Parameters**:
- `file_path` (character): Path to section CSV file

**Returns**: Numeric cluster ID

**Example**: `extract_cluster_id("results/sections/LAS_sec_42.csv")` returns `42`

---

### `load_section_data(sections_dir, cluster_id)`

Loads section data for a specific cluster.

**Parameters**:
- `sections_dir` (character): Directory with sections
- `cluster_id` (integer): Cluster ID

**Returns**: Data frame with section data

---

### `setup_parallel_cluster(n_cores)`

Sets up parallel processing cluster.

**Parameters**:
- `n_cores` (integer): Number of cores (default: auto-detect)

**Returns**: Cluster object

---

### `get_section_files(sections_dir)`

Gets list of all section files.

**Parameters**:
- `sections_dir` (character): Directory with sections

**Returns**: Character vector of file paths

---

## Configuration (`config.R`)

The configuration file defines all paths, parameters, and settings used throughout the workflow. Key variables:

- Paths: `LAS_FILE`, `RESULTS_DIR`, `FIGURES_DIR`
- Clustering: `NUM_CLUSTERS`, `N_START`, `CLUSTERING_SEED`
- Model: `MODEL`, `CORRMODEL`, `START_PARAMS`, `FIXED_PARAMS`
- Optimization: `OPTIMIZER`, `LOWER_BOUNDS`, `UPPER_BOUNDS`, `NEIGHB`
- Health indicators: `I0_QUANTILE`, `TMR_QUANTILE`, `EA_BINS`
- Parallel: `N_CORES`, `USE_PARALLEL`

---

## Main Script (`run_analysis.R`)

The main workflow script orchestrates all steps. It can be run with command-line flags:

- `--skip-segmentation`: Skip data segmentation (if sections exist)
- `--skip-estimation`: Skip parameter estimation (if estimates exist)
- `--skip-indicators`: Skip health indicator computation
- `--skip-figures`: Skip figure generation

Example:
```bash
Rscript run_analysis.R --skip-segmentation --skip-estimation
```

