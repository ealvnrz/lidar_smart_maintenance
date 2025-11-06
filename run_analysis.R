# =============================================================================
# Main Workflow Script
# =============================================================================
# This script orchestrates the complete Weibull Random Field analysis pipeline

# ---- Setup ----
# Load configuration
source("config.R")

# Load required libraries
required_packages <- c(
  "lidR", "GeoModels", "readr", "parallel", "ggplot2", 
  "plotly", "MASS", "dplyr", "tidyr", "fields"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Required package not installed:", pkg))
  }
}

# Load all R functions
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (r_file in r_files) {
  source(r_file)
}

# Set random seed for reproducibility
set.seed(ANALYSIS_SEED)

# ---- Command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
skip_segmentation <- "--skip-segmentation" %in% args
skip_estimation <- "--skip-estimation" %in% args
skip_indicators <- "--skip-indicators" %in% args
skip_figures <- "--skip-figures" %in% args

# ---- Step 1: Data Processing and Segmentation ----
if (!skip_segmentation) {
  cat("\n=== Step 1: Data Processing and Segmentation ===\n")
  
  # Check if sections already exist
  if (check_sections_exist(RESULTS_SECTIONS_DIR)) {
    cat("Sections already exist. Skipping segmentation.\n")
    cat("(Use --skip-segmentation to suppress this check)\n")
  } else {
    cat("Processing LAS file and segmenting road...\n")
    
    data <- process_and_segment_data(
      las_file = LAS_FILE,
      output_dir = RESULTS_SECTIONS_DIR,
      num_clusters = NUM_CLUSTERS,
      n_start = N_START,
      seed = CLUSTERING_SEED
    )
    
    cat("Segmentation complete. Created", NUM_CLUSTERS, "sections.\n")
  }
} else {
  cat("Skipping segmentation (--skip-segmentation flag)\n")
}

# ---- Step 2: Marginal Analysis (Optional Diagnostic) ----
if (!skip_figures) {
  cat("\n=== Step 2: Marginal Analysis (Diagnostic) ===\n")
  
  # Check marginal fit for a sample cluster
  sample_cluster <- 2
  cat("Checking marginal fit for cluster", sample_cluster, "...\n")
  
  marginal_result <- check_marginal_fit(
    sections_dir = RESULTS_SECTIONS_DIR,
    cluster_id = sample_cluster,
    plot = TRUE
  )
  
  if (!is.null(marginal_result)) {
    cat("Marginal fit check complete.\n")
    
    # Save diagnostic plot
    if (!is.null(marginal_result$plots)) {
      save_plot(
        marginal_result$plots$histogram,
        file.path(FIGURES_RESULTS_DIR, "hist_weib_c_2.png")
      )
    }
  }
}

# ---- Step 3: Spatial Parameter Estimation ----
if (!skip_estimation) {
  cat("\n=== Step 3: Spatial Parameter Estimation ===\n")
  
  # Check if estimation already exists
  estimation_file <- file.path(RESULTS_ESTIMATION_DIR, "param_estimation.csv")
  
  if (file.exists(estimation_file)) {
    cat("Parameter estimation already exists. Skipping estimation.\n")
    cat("(Use --skip-estimation to suppress this check)\n")
    
    param_estimation <- readr::read_csv(estimation_file, show_col_types = FALSE)
    cat("Loaded", nrow(param_estimation), "parameter estimates.\n")
    
  } else {
    cat("Estimating WRF parameters for all sections...\n")
    
    param_estimation <- estimate_all_sections(
      sections_dir = RESULTS_SECTIONS_DIR,
      output_file = estimation_file,
      parallel = USE_PARALLEL,
      n_cores = N_CORES,
      max_files = MAX_FILES
    )
    
    cat("Parameter estimation complete.\n")
  }
} else {
  cat("Skipping estimation (--skip-estimation flag)\n")
  
  # Still try to load existing estimation
  estimation_file <- file.path(RESULTS_ESTIMATION_DIR, "param_estimation.csv")
  if (file.exists(estimation_file)) {
    param_estimation <- readr::read_csv(estimation_file, show_col_types = FALSE)
  } else {
    stop("Estimation file not found and estimation was skipped.")
  }
}

# ---- Step 4: Health Indicators ----
if (!skip_indicators) {
  cat("\n=== Step 4: Computing Health Indicators ===\n")
  
  # Check if indicators already exist
  indicators_file <- file.path(RESULTS_HEALTH_DIR, "health_indicators.csv")
  
  if (file.exists(indicators_file)) {
    cat("Health indicators already exist. Skipping computation.\n")
    cat("(Use --skip-indicators to suppress this check)\n")
    
    health_results <- list(
      indicators = readr::read_csv(indicators_file, show_col_types = FALSE)
    )
    
    # Load flagged clusters
    if (dir.exists(RESULTS_FLAGGED_DIR)) {
      health_results$flagged_eli <- readr::read_csv(
        file.path(RESULTS_FLAGGED_DIR, "flagged_ELI.csv"),
        show_col_types = FALSE
      )
      health_results$flagged_tmr <- readr::read_csv(
        file.path(RESULTS_FLAGGED_DIR, "flagged_TMR.csv"),
        show_col_types = FALSE
      )
      health_results$flagged_ea <- readr::read_csv(
        file.path(RESULTS_FLAGGED_DIR, "flagged_EA.csv"),
        show_col_types = FALSE
      )
    }
    
  } else {
    cat("Computing health indicators...\n")
    
    health_results <- compute_all_indicators(
      sections_dir = RESULTS_SECTIONS_DIR,
      param_estimation_file = estimation_file,
      output_dir = RESULTS_HEALTH_DIR,
      I0_quantile = I0_QUANTILE,
      tmr_quantile = TMR_QUANTILE,
      ea_bins = EA_BINS,
      flag_intensity_percentile = FLAG_INTENSITY_PERCENTILE,
      flag_mean_percentile = FLAG_MEAN_PERCENTILE,
      parallel = USE_PARALLEL
    )
    
    cat("Health indicators computation complete.\n")
    cat("Flagged clusters - ELI:", sum(health_results$flagged_eli$flagged, na.rm = TRUE),
        "TMR:", sum(health_results$flagged_tmr$flagged, na.rm = TRUE),
        "EA:", sum(health_results$flagged_ea$flagged, na.rm = TRUE), "\n")
  }
} else {
  cat("Skipping health indicators (--skip-indicators flag)\n")
}

# ---- Step 5: Generate Figures (Optional) ----
if (!skip_figures) {
  cat("\n=== Step 5: Generating Figures ===\n")
  
  cat("Note: Figures will be modified later. Generating basic diagnostics only.\n")
  
  # Parameter boxplots
  if (exists("param_estimation")) {
    cat("Creating parameter boxplots...\n")
    
    p_scale <- plot_parameter_boxplots(param_estimation, "scale", 
                                       fill_color = "lightcoral")
    p_shape <- plot_parameter_boxplots(param_estimation, "shape",
                                       fill_color = "lightgreen")
    
    save_plot(p_scale, file.path(FIGURES_RESULTS_DIR, "scale_boxplot.png"))
    save_plot(p_shape, file.path(FIGURES_RESULTS_DIR, "shape_boxplot.png"))
  }
  
  # Sample semivariogram (if estimation results available)
  if (exists("param_estimation")) {
    sample_cluster <- 73
    cat("Creating semivariogram for cluster", sample_cluster, "...\n")
    
    tryCatch({
      section_data <- load_section_data(RESULTS_SECTIONS_DIR, sample_cluster)
      coords <- cbind(section_data$x, section_data$y)
      
      # Note: Full semivariogram plotting requires a fitted GeoFit object
      # This would require re-fitting or loading saved fit objects
      # For now, we skip this step
      cat("Semivariogram plotting requires fit objects (to be implemented)\n")
      
    }, error = function(e) {
      cat("Could not create semivariogram:", e$message, "\n")
    })
  }
  
} else {
  cat("Skipping figures (--skip-figures flag)\n")
}

# ---- Summary ----
cat("\n=== Analysis Complete ===\n")
cat("Results saved to:", RESULTS_DIR, "\n")
cat("Figures saved to:", FIGURES_RESULTS_DIR, "\n")

# Print session info for reproducibility
cat("\n=== Session Info ===\n")
print(sessionInfo())

