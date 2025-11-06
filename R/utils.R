# =============================================================================
# Utility Functions
# =============================================================================
# Helper functions used across the workflow

#' Extract cluster ID from file path
#'
#' @param file_path Character string, path to section CSV file
#' @return Numeric cluster ID
#'
#' @examples
#' extract_cluster_id("results/sections/LAS_sec_42.csv")  # Returns 42
extract_cluster_id <- function(file_path) {
  file_name <- basename(file_path)
  cluster_id <- gsub("LAS_sec_|\\.csv", "", file_name)
  cluster_id <- as.numeric(cluster_id)
  
  if (is.na(cluster_id)) {
    stop(paste("Could not extract cluster ID from file:", file_path))
  }
  
  return(cluster_id)
}

#' Load section data for a specific cluster
#'
#' @param sections_dir Character string, directory containing section CSV files
#' @param cluster_id Numeric, cluster ID to load
#' @return Data frame with section data
#'
#' @examples
#' section_data <- load_section_data("results/sections", cluster_id = 42)
load_section_data <- function(sections_dir, cluster_id) {
  file_path <- file.path(sections_dir, paste0("LAS_sec_", cluster_id, ".csv"))
  
  if (!file.exists(file_path)) {
    stop(paste("Section file not found:", file_path))
  }
  
  data <- readr::read_csv(file_path, show_col_types = FALSE)
  return(data)
}

#' Setup parallel processing cluster
#'
#' @param n_cores Integer, number of cores to use (default from config)
#' @return Cluster object for parallel processing
#'
#' @examples
#' cl <- setup_parallel_cluster(n_cores = 4)
#' # ... do parallel work ...
#' parallel::stopCluster(cl)
setup_parallel_cluster <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
  }
  
  if (n_cores < 1) {
    n_cores <- 1
  }
  
  cl <- parallel::makeCluster(n_cores)
  return(cl)
}

#' Validate that required paths exist
#'
#' @param paths Character vector of paths to check
#' @return Logical, TRUE if all paths exist, FALSE otherwise
validate_paths <- function(paths) {
  all_exist <- all(file.exists(paths) | dir.exists(paths))
  
  if (!all_exist) {
    missing <- paths[!(file.exists(paths) | dir.exists(paths))]
    warning(paste("Missing paths:", paste(missing, collapse = ", ")))
  }
  
  return(all_exist)
}

#' Check if section files exist
#'
#' @param sections_dir Character string, directory containing sections
#' @param expected_count Integer, expected number of sections (default: 250)
#' @return Logical, TRUE if all expected files exist
check_sections_exist <- function(sections_dir, expected_count = 250) {
  pattern <- paste0("LAS_sec_", 1:expected_count, "\\.csv")
  files_exist <- sapply(pattern, function(p) {
    any(grepl(p, list.files(sections_dir)))
  })
  
  all_exist <- all(files_exist)
  
  if (!all_exist) {
    missing <- which(!files_exist)
    warning(paste("Missing sections:", paste(missing, collapse = ", ")))
  }
  
  return(all_exist)
}

#' Get list of all section files
#'
#' @param sections_dir Character string, directory containing sections
#' @return Character vector of full file paths
get_section_files <- function(sections_dir) {
  files <- list.files(sections_dir, pattern = "LAS_sec_.*\\.csv$", 
                      full.names = TRUE)
  return(sort(files))
}

