# =============================================================================
# Health Indicators Functions
# =============================================================================
# Functions for computing health indicators (ELI, TMR, EA) and flagging degraded clusters

#' Compute Energy-Loss Index (ELI)
#'
#' Computes the Energy-Loss Index for a cluster, which measures the average
#' logarithmic deviation of intensities from a reference value.
#'
#' @param intensity Numeric vector of intensity values
#' @param I0 Numeric, global reference intensity (typically 99th percentile)
#' @return Numeric, ELI value
#'
#' @details
#' Formula: ELI_k = -(1/n_k) * sum(ln(I_i / I_0))
#'
#' @examples
#' eli <- compute_eli(cluster_data$intensity_scaled, I0 = 0.95)
compute_eli <- function(intensity, I0) {
  if (length(intensity) == 0) {
    stop("Intensity vector is empty")
  }
  
  if (I0 <= 0) {
    stop("I0 must be positive")
  }
  
  # Filter positive values
  intensity_positive <- intensity[intensity > 0]
  
  if (length(intensity_positive) == 0) {
    warning("No positive intensity values, returning NA")
    return(NA)
  }
  
  # Compute ELI
  n_k <- length(intensity_positive)
  eli <- -(1 / n_k) * sum(log(intensity_positive / I0))
  
  return(eli)
}

#' Compute Tail-to-Mean Ratio (TMR)
#'
#' Computes the ratio of a tail quantile to the mean intensity,
#' indicating the presence of extreme values.
#'
#' @param intensity Numeric vector of intensity values
#' @param quantile Numeric, quantile level (default: 0.95)
#' @return Numeric, TMR value
#'
#' @details
#' Formula: TMR_{q,k} = q_{q,k} / mean(I_k)
#'
#' @examples
#' tmr <- compute_tmr(cluster_data$intensity_scaled, quantile = 0.95)
compute_tmr <- function(intensity, quantile = 0.95) {
  if (length(intensity) == 0) {
    stop("Intensity vector is empty")
  }
  
  if (quantile < 0 || quantile > 1) {
    stop("Quantile must be between 0 and 1")
  }
  
  # Filter positive values
  intensity_positive <- intensity[intensity > 0]
  
  if (length(intensity_positive) == 0) {
    warning("No positive intensity values, returning NA")
    return(NA)
  }
  
  # Compute TMR
  q_quantile <- quantile(intensity_positive, probs = quantile)
  mean_intensity <- mean(intensity_positive)
  
  if (mean_intensity == 0) {
    warning("Mean intensity is zero, returning NA")
    return(NA)
  }
  
  tmr <- as.numeric(q_quantile) / mean_intensity
  
  return(tmr)
}

#' Compute Entropy of Attenuation (EA)
#'
#' Computes the entropy of the intensity distribution, measuring
#' the uncertainty or variability in the intensity values.
#'
#' @param intensity Numeric vector of intensity values
#' @param bins Integer, number of bins for histogram (default: 30)
#' @return Numeric, EA value
#'
#' @details
#' Formula: EA_k = -sum(p_{j,k} * ln(p_{j,k}))
#' where p_{j,k} are histogram proportions
#'
#' @examples
#' ea <- compute_ea(cluster_data$intensity_scaled, bins = 30)
compute_ea <- function(intensity, bins = 30) {
  if (length(intensity) == 0) {
    stop("Intensity vector is empty")
  }
  
  # Filter positive values
  intensity_positive <- intensity[intensity > 0]
  
  if (length(intensity_positive) == 0) {
    warning("No positive intensity values, returning NA")
    return(NA)
  }
  
  # Create histogram
  hist_result <- hist(intensity_positive, breaks = bins, plot = FALSE)
  
  # Compute proportions (avoiding zero counts)
  counts <- hist_result$counts
  total <- sum(counts)
  
  if (total == 0) {
    warning("No counts in histogram, returning NA")
    return(NA)
  }
  
  # Compute proportions (add small epsilon to avoid log(0))
  p_j <- (counts + 1e-10) / (total + bins * 1e-10)
  
  # Compute entropy
  p_j_positive <- p_j[p_j > 0]
  ea <- -sum(p_j_positive * log(p_j_positive))
  
  return(ea)
}

#' Compute all health indicators for a cluster
#'
#' Convenience function to compute ELI, TMR, and EA for a single cluster.
#'
#' @param intensity Numeric vector of intensity values
#' @param I0 Numeric, global reference intensity
#' @param tmr_quantile Numeric, quantile for TMR (default: 0.95)
#' @param ea_bins Integer, number of bins for EA (default: 30)
#' @return Data frame with cluster_id and all three indicators
#'
#' @examples
#' indicators <- compute_cluster_indicators(cluster_data$intensity_scaled, I0 = 0.95)
compute_cluster_indicators <- function(intensity, 
                                      I0,
                                      tmr_quantile = 0.95,
                                      ea_bins = 30) {
  eli <- compute_eli(intensity, I0)
  tmr <- compute_tmr(intensity, quantile = tmr_quantile)
  ea <- compute_ea(intensity, bins = ea_bins)
  
  result <- data.frame(
    ELI = eli,
    TMR = tmr,
    EA = ea,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Flag degraded clusters based on health indicators
#'
#' Flags clusters as degraded based on intensity and mean intensity thresholds.
#'
#' @param indicators_df Data frame with columns: cluster, mean_intensity, and indicator columns
#' @param intensity_threshold Numeric, percentile threshold for intensity (default: 0.75)
#' @param mean_intensity_threshold Numeric, percentile threshold for mean intensity (default: 0.25)
#' @param indicator_name Character string, name of indicator column to use for flagging
#' @return Data frame with added 'flagged' logical column
#'
#' @details
#' A cluster is flagged if:
#' - indicator_value > intensity_threshold percentile AND
#' - mean_intensity < mean_intensity_threshold percentile
#'
#' @examples
#' flagged_df <- flag_degraded_clusters(indicators_df, indicator_name = "ELI")
flag_degraded_clusters <- function(indicators_df,
                                  intensity_threshold = 0.75,
                                  mean_intensity_threshold = 0.25,
                                  indicator_name = "ELI") {
  
  if (!indicator_name %in% colnames(indicators_df)) {
    stop(paste("Indicator column", indicator_name, "not found in data frame"))
  }
  
  if (!"mean_intensity" %in% colnames(indicators_df)) {
    stop("mean_intensity column not found in data frame")
  }
  
  # Compute thresholds
  indicator_threshold <- quantile(indicators_df[[indicator_name]], 
                                  probs = intensity_threshold, 
                                  na.rm = TRUE)
  mean_threshold <- quantile(indicators_df$mean_intensity, 
                            probs = mean_intensity_threshold, 
                            na.rm = TRUE)
  
  # Flag clusters
  indicators_df$flagged <- (indicators_df[[indicator_name]] > indicator_threshold) &
                           (indicators_df$mean_intensity < mean_threshold)
  
  # Add threshold information
  indicators_df$indicator_threshold <- indicator_threshold
  indicators_df$mean_threshold <- mean_threshold
  
  return(indicators_df)
}

#' Compute health indicators for all clusters
#'
#' Main workflow function to compute ELI, TMR, and EA for all sections,
#' and flag degraded clusters.
#'
#' @param sections_dir Character string, directory containing section CSV files
#' @param param_estimation_file Character string, path to parameter estimation CSV
#' @param output_dir Character string, directory to save results
#' @param I0_quantile Numeric, quantile for global reference intensity (default: 0.99)
#' @param tmr_quantile Numeric, quantile for TMR (default: 0.95)
#' @param ea_bins Integer, number of bins for EA (default: 30)
#' @param flag_intensity_percentile Numeric, percentile for flagging (default: 0.75)
#' @param flag_mean_percentile Numeric, percentile for mean flagging (default: 0.25)
#' @param parallel Logical, whether to use parallel processing (default: TRUE)
#' @return List with indicators data frame and flagged clusters data frame
#'
#' @examples
#' results <- compute_all_indicators("results/sections", 
#'                                   "results/estimation/param_estimation.csv",
#'                                   "results/health_indicators")
compute_all_indicators <- function(sections_dir,
                                  param_estimation_file,
                                  output_dir,
                                  I0_quantile = 0.99,
                                  tmr_quantile = 0.95,
                                  ea_bins = 30,
                                  flag_intensity_percentile = 0.75,
                                  flag_mean_percentile = 0.25,
                                  parallel = TRUE) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Get all section files
  csv_files <- get_section_files(sections_dir)
  
  if (length(csv_files) == 0) {
    stop(paste("No section files found in", sections_dir))
  }
  
  cat("Computing health indicators for", length(csv_files), "clusters...\n")
  
  # First pass: collect all intensities to compute global I0
  cat("Computing global reference intensity (I0)...\n")
  all_intensities <- numeric(0)
  
  for (file_path in csv_files) {
    data <- readr::read_csv(file_path, show_col_types = FALSE)
    all_intensities <- c(all_intensities, data$intensity_scaled[data$intensity_scaled > 0])
  }
  
  I0 <- quantile(all_intensities, probs = I0_quantile, na.rm = TRUE)
  cat("I0 (", I0_quantile * 100, "th percentile):", I0, "\n")
  
  # Function to process a single cluster
  process_cluster <- function(file_path) {
    tryCatch({
      cluster_id <- extract_cluster_id(file_path)
      data <- readr::read_csv(file_path, show_col_types = FALSE)
      
      # Compute indicators
      indicators <- compute_cluster_indicators(
        intensity = data$intensity_scaled,
        I0 = I0,
        tmr_quantile = tmr_quantile,
        ea_bins = ea_bins
      )
      
      # Add cluster and summary statistics
      indicators$cluster <- cluster_id
      indicators$mean_intensity <- mean(data$intensity_scaled[data$intensity_scaled > 0], na.rm = TRUE)
      indicators$n_points <- nrow(data)
      
      return(indicators)
      
    }, error = function(e) {
      cat("Error processing cluster from", file_path, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Process clusters
  if (parallel && length(csv_files) > 1) {
    n_cores <- parallel::detectCores() - 1
    cl <- setup_parallel_cluster(n_cores)
    
    # Load required packages on cluster nodes
    parallel::clusterEvalQ(cl, {
      library(readr)
    })
    
    parallel::clusterExport(cl,
                           varlist = c("process_cluster",
                                      "extract_cluster_id",
                                      "compute_cluster_indicators",
                                      "compute_eli",
                                      "compute_tmr",
                                      "compute_ea",
                                      "I0",
                                      "tmr_quantile",
                                      "ea_bins"),
                           envir = environment())
    
    results_list <- parallel::parLapply(cl, csv_files, process_cluster)
    parallel::stopCluster(cl)
    
  } else {
    results_list <- lapply(csv_files, process_cluster)
  }
  
  # Combine results
  results_list <- results_list[!sapply(results_list, is.null)]
  indicators_df <- do.call(rbind, results_list)
  
  # Flag degraded clusters for each indicator
  indicators_df_eli <- flag_degraded_clusters(
    indicators_df,
    intensity_threshold = flag_intensity_percentile,
    mean_intensity_threshold = flag_mean_percentile,
    indicator_name = "ELI"
  )
  
  indicators_df_tmr <- flag_degraded_clusters(
    indicators_df,
    intensity_threshold = flag_intensity_percentile,
    mean_intensity_threshold = flag_mean_percentile,
    indicator_name = "TMR"
  )
  
  indicators_df_ea <- flag_degraded_clusters(
    indicators_df,
    intensity_threshold = flag_intensity_percentile,
    mean_intensity_threshold = flag_mean_percentile,
    indicator_name = "EA"
  )
  
  # Save results
  indicators_file <- file.path(output_dir, "health_indicators.csv")
  readr::write_csv(indicators_df, indicators_file)
  cat("Saved indicators to", indicators_file, "\n")
  
  # Save flagged clusters to separate directory
  # Get flagged_clusters directory from config
  if (!exists("RESULTS_FLAGGED_DIR")) {
    # Try to get it from config if not already loaded
    if (file.exists("config.R")) {
      source("config.R", local = TRUE)
    }
  }
  
  # Use RESULTS_FLAGGED_DIR if available, otherwise create flagged subdirectory
  if (exists("RESULTS_FLAGGED_DIR")) {
    flagged_dir <- RESULTS_FLAGGED_DIR
  } else {
    flagged_dir <- file.path(output_dir, "flagged")
  }
  
  if (!dir.exists(flagged_dir)) {
    dir.create(flagged_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  readr::write_csv(indicators_df_eli, file.path(flagged_dir, "flagged_ELI.csv"))
  readr::write_csv(indicators_df_tmr, file.path(flagged_dir, "flagged_TMR.csv"))
  readr::write_csv(indicators_df_ea, file.path(flagged_dir, "flagged_EA.csv"))
  
  cat("Saved flagged clusters to", flagged_dir, "\n")
  
  return(list(
    indicators = indicators_df,
    flagged_eli = indicators_df_eli,
    flagged_tmr = indicators_df_tmr,
    flagged_ea = indicators_df_ea
  ))
}

