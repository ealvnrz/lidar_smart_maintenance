# =============================================================================
# Configuration File for Weibull Random Field Analysis
# =============================================================================
# This file contains all configuration parameters, paths, and settings
# for the Weibull Random Field workflow.

# ---- Paths ----
# Base directory (assumes script is run from project root)
BASE_DIR <- getwd()

# Input data paths
DATA_DIR <- file.path(BASE_DIR, "data")
DATA_RAW_DIR <- file.path(DATA_DIR, "raw")
LAS_FILE <- file.path(DATA_RAW_DIR, "minasur.las")

# Output paths
RESULTS_DIR <- file.path(BASE_DIR, "results")
RESULTS_SECTIONS_DIR <- file.path(RESULTS_DIR, "sections")
RESULTS_ESTIMATION_DIR <- file.path(RESULTS_DIR, "estimation")
RESULTS_HEALTH_DIR <- file.path(RESULTS_DIR, "health_indicators")
RESULTS_FLAGGED_DIR <- file.path(RESULTS_DIR, "flagged_clusters")

FIGURES_DIR <- file.path(BASE_DIR, "figures")
FIGURES_RESULTS_DIR <- file.path(FIGURES_DIR, "results")

# ---- Clustering Parameters ----
NUM_CLUSTERS <- 250
N_START <- 25
CLUSTERING_SEED <- 234528

# ---- Model Parameters ----
MODEL <- "Weibull"
CORRMODEL <- "Wend0"

# Fixed parameters for WRF estimation
FIXED_PARAMS <- list(
  nugget = 0,
  power2 = 4
)

# Starting parameters for WRF estimation
START_PARAMS <- list(
  mean = 0.5,
  scale = 10,
  shape = 5
)

# Optimization bounds
I_BOUND <- 1000
LOWER_BOUNDS <- list(
  mean = -I_BOUND,
  shape = 0.00001,
  scale = 0
)
UPPER_BOUNDS <- list(
  mean = I_BOUND,
  shape = I_BOUND,
  scale = I_BOUND
)

# Optimizer settings
OPTIMIZER <- "nlminb"  # Options: "nlminb" or default (uses GeoModels default)
NEIGHB <- 10  # Number of neighbors for pairwise likelihood

# ---- Health Indicator Parameters ----
# Global intensity reference (99th percentile)
I0_QUANTILE <- 0.99

# TMR quantile
TMR_QUANTILE <- 0.95

# Entropy of Attenuation bins
EA_BINS <- 30

# Flagging thresholds (percentiles)
FLAG_INTENSITY_PERCENTILE <- 0.75  # Flag if intensity > 75th percentile
FLAG_MEAN_PERCENTILE <- 0.25  # Flag if mean intensity < 25th percentile

# ---- Parallel Processing ----
N_CORES <- parallel::detectCores() - 5  # Leave five cores free
USE_PARALLEL <- TRUE

# ---- Random Seeds ----
ANALYSIS_SEED <- 323622

# ---- Other Settings ----
MAX_FILES <- 999999  # Maximum number of files to process (safety limit)

# =============================================================================
# Helper function to validate all paths exist
# =============================================================================
validate_config <- function() {
  required_dirs <- c(
    DATA_RAW_DIR,
    RESULTS_DIR,
    RESULTS_SECTIONS_DIR,
    RESULTS_ESTIMATION_DIR,
    RESULTS_HEALTH_DIR,
    RESULTS_FLAGGED_DIR,
    FIGURES_DIR,
    FIGURES_RESULTS_DIR
  )
  
  missing_dirs <- required_dirs[!dir.exists(required_dirs)]
  
  if (length(missing_dirs) > 0) {
    cat("Creating missing directories:\n")
    for (dir in missing_dirs) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      cat("  Created:", dir, "\n")
    }
  }
  
  if (!file.exists(LAS_FILE)) {
    warning(paste("LAS file not found at:", LAS_FILE))
  }
  
  invisible(TRUE)
}

# Validate paths on source
validate_config()

