# =============================================================================
# Spatial Estimation Functions
# =============================================================================
# Functions for estimating Weibull Random Field parameters using pairwise likelihood

#' Estimate WRF parameters for a single section
#'
#' Estimates Weibull Random Field parameters (mean, scale, shape, correlation range)
#' for a single road section using pairwise likelihood with Wendland correlation.
#'
#' @param section_data Data frame with section data (must contain x, y, intensity_scaled)
#' @param coords Matrix or data frame with coordinates (n x 2). If NULL, uses section_data$x, section_data$y
#' @param start_params List with starting parameters (mean, scale, shape)
#' @param fixed_params List with fixed parameters (nugget, power2)
#' @param lower_bounds List with lower bounds for optimization
#' @param upper_bounds List with upper bounds for optimization
#' @param corrmodel Character string, correlation model (default: "Wend0")
#' @param model Character string, distribution model (default: "Weibull")
#' @param optimizer Character string, optimizer to use (default: "nlminb")
#' @param neighb Integer, number of neighbors for pairwise likelihood (default: 10)
#' @return List with estimated parameters and fit statistics
#'
#' @examples
#' fit <- estimate_wrf_parameters(section_data, 
#'                                 start_params = list(mean = 0.5, scale = 10, shape = 5))
estimate_wrf_parameters <- function(section_data,
                                    coords = NULL,
                                    start_params = list(mean = 0.5, scale = 10, shape = 5),
                                    fixed_params = list(nugget = 0, power2 = 4),
                                    lower_bounds = list(mean = -1000, shape = 0.00001, scale = 0),
                                    upper_bounds = list(mean = 1000, shape = 1000, scale = 1000),
                                    corrmodel = "Wend0",
                                    model = "Weibull",
                                    optimizer = "nlminb",
                                    neighb = 10) {
  
  # Extract coordinates if not provided
  if (is.null(coords)) {
    if (!all(c("x", "y") %in% colnames(section_data))) {
      stop("section_data must contain x and y columns, or coords must be provided")
    }
    coords <- cbind(section_data$x, section_data$y)
  }
  
  # Check intensity data
  if (!"intensity_scaled" %in% colnames(section_data)) {
    stop("section_data must contain intensity_scaled column")
  }
  
  intensity_data <- section_data$intensity_scaled
  
  # Prepare GeoFit arguments
  fit_args <- list(
    data = intensity_data,
    coordx = coords,
    corrmodel = corrmodel,
    model = model,
    start = start_params,
    fixed = fixed_params,
    likelihood = "Marginal",
    type = "Pairwise",
    sensitivity = TRUE,
    neighb = neighb
  )
  
  # Add optimizer and bounds if specified
  if (!is.null(optimizer)) {
    fit_args$optimizer <- optimizer
    fit_args$lower <- lower_bounds
    fit_args$upper <- upper_bounds
  }
  
  # Fit the model
  tryCatch({
    fit <- do.call(GeoModels::GeoFit, fit_args)
    
    # Extract results
    result <- list(
      cluster = section_data$Cluster[1],
      mean = fit$param$mean,
      scale = fit$param$scale,
      shape = fit$param$shape,
      logCompLik = fit$logCompLik,
      fit_object = fit,
      success = TRUE
    )
    
    return(result)
    
  }, error = function(e) {
    warning(paste("Error fitting cluster", section_data$Cluster[1], ":", e$message))
    return(list(
      cluster = section_data$Cluster[1],
      mean = NA,
      scale = NA,
      shape = NA,
      logCompLik = NA,
      success = FALSE,
      error = e$message
    ))
  })
}

#' Get estimation configuration from config.R
#'
#' Convenience function to retrieve estimation settings.
#'
#' @return List with estimation configuration
get_estimation_config <- function() {
  # Source config if not already loaded
  if (!exists("START_PARAMS")) {
    source("config.R")
  }
  
  config <- list(
    start_params = START_PARAMS,
    fixed_params = FIXED_PARAMS,
    lower_bounds = LOWER_BOUNDS,
    upper_bounds = UPPER_BOUNDS,
    corrmodel = CORRMODEL,
    model = MODEL,
    optimizer = OPTIMIZER,
    neighb = NEIGHB
  )
  
  return(config)
}

#' Estimate WRF parameters for all sections
#'
#' Processes all section files in parallel and estimates WRF parameters.
#' Results are saved to a CSV file.
#'
#' @param sections_dir Character string, directory containing section CSV files
#' @param output_file Character string, path to output CSV file
#' @param parallel Logical, whether to use parallel processing (default: TRUE)
#' @param n_cores Integer, number of cores for parallel processing (default: auto-detect)
#' @param max_files Integer, maximum number of files to process (safety limit)
#' @return Data frame with parameter estimates for all clusters
#'
#' @examples
#' results <- estimate_all_sections("results/sections", "results/estimation/param_estimation.csv")
estimate_all_sections <- function(sections_dir,
                                  output_file,
                                  parallel = TRUE,
                                  n_cores = NULL,
                                  max_files = 999999) {
  
  # Get section files
  csv_files <- get_section_files(sections_dir)
  
  if (length(csv_files) == 0) {
    stop(paste("No section files found in", sections_dir))
  }
  
  # Limit number of files
  if (length(csv_files) > max_files) {
    warning(paste("Limiting processing to", max_files, "files"))
    csv_files <- csv_files[1:max_files]
  }
  
  # Get estimation configuration
  config <- get_estimation_config()
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Initialize output file
  if (file.exists(output_file)) {
    file.remove(output_file)
  }
  
  # Function to process a single file
  process_file <- function(file_path) {
    tryCatch({
      data <- readr::read_csv(file_path, show_col_types = FALSE)
      
      # Estimate parameters
      result <- estimate_wrf_parameters(
        section_data = data,
        start_params = config$start_params,
        fixed_params = config$fixed_params,
        lower_bounds = config$lower_bounds,
        upper_bounds = config$upper_bounds,
        corrmodel = config$corrmodel,
        model = config$model,
        optimizer = config$optimizer,
        neighb = config$neighb
      )
      
      # Create result row
      new_row <- data.frame(
        cluster = result$cluster,
        mean = result$mean,
        scale = result$scale,
        shape = result$shape,
        logCompLik = result$logCompLik,
        stringsAsFactors = FALSE
      )
      
      # Append to file
      write.table(new_row, 
                 file = output_file, 
                 sep = ",", 
                 row.names = FALSE, 
                 col.names = !file.exists(output_file), 
                 append = file.exists(output_file))
      
      cat("Processed cluster:", result$cluster, "\n")
      
      return(new_row)
      
    }, error = function(e) {
      cat("Error processing file:", file_path, "-", e$message, "\n")
      return(NULL)
    })
  }
  
  # Process files
  if (parallel && length(csv_files) > 1) {
    # Set up parallel processing
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    
    cl <- setup_parallel_cluster(n_cores)
    
    # Load required packages on cluster nodes
    parallel::clusterEvalQ(cl, {
      library(GeoModels)
      library(readr)
    })
    
    # Export necessary variables
    parallel::clusterExport(cl, 
                           varlist = c("process_file", 
                                      "estimate_wrf_parameters",
                                      "get_estimation_config",
                                      "config",
                                      "output_file"),
                           envir = environment())
    
    # Process in parallel
    results <- parallel::parLapply(cl, csv_files, process_file)
    
    # Stop cluster
    parallel::stopCluster(cl)
    
  } else {
    # Process sequentially
    results <- lapply(csv_files, process_file)
  }
  
  # Read and return results
  if (file.exists(output_file)) {
    results_df <- readr::read_csv(output_file, show_col_types = FALSE)
    cat("Estimated parameters for", nrow(results_df), "clusters\n")
    return(results_df)
  } else {
    warning("No results file created")
    return(NULL)
  }
}

