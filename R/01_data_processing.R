# =============================================================================
# Data Processing Functions
# =============================================================================
# Functions for loading LAS data, normalization, and road segmentation

#' Load and normalize LAS point cloud data
#'
#' Loads a LAS file and extracts coordinates, RGB, and intensity data.
#' Normalizes all values to [0, 1] range for clustering.
#'
#' @param file_path Character string, path to LAS file
#' @return Data frame with original and normalized values
#'
#' @examples
#' data <- load_las_data("data/raw/minasur.las")
load_las_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("LAS file not found:", file_path))
  }
  
  # Load LAS file
  db <- lidR::readLAS(file_path, select = "xyzRGBi")
  
  # Extract coordinates and attributes
  x <- db[['X']]
  y <- db[['Y']]
  z <- db[['Z']]
  intensity <- db[['Intensity']]
  r <- db[['R']]
  g <- db[['G']]
  b <- db[['B']]
  
  # Normalize to [0, 1]
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  y_scaled <- (y - min(y)) / (max(y) - min(y))
  z_scaled <- (z - min(z)) / (max(z) - min(z))
  intensity_scaled <- (intensity - min(intensity)) / (max(intensity) - min(intensity))
  r_scaled <- (r - min(r)) / (max(r) - min(r))
  g_scaled <- (g - min(g)) / (max(g) - min(g))
  b_scaled <- (b - min(b)) / (max(b) - min(b))
  
  # Create data frame
  data <- data.frame(
    x = x,
    y = y,
    z = z,
    x_scaled = x_scaled,
    y_scaled = y_scaled,
    z_scaled = z_scaled,
    intensity = intensity,
    intensity_scaled = intensity_scaled,
    r = r_scaled,
    g = g_scaled,
    b = b_scaled
  )
  
  return(data)
}

#' Segment road into clusters using K-means
#'
#' Performs K-means clustering on normalized spatial coordinates to partition
#' the road into homogeneous segments.
#'
#' @param data Data frame with normalized coordinates (x_scaled, y_scaled, z_scaled)
#' @param num_clusters Integer, number of clusters (default: 250)
#' @param n_start Integer, number of random starts for K-means (default: 25)
#' @param seed Integer, random seed for reproducibility
#' @return Data frame with added Cluster column
#'
#' @examples
#' data_clustered <- segment_road(data, num_clusters = 250, n_start = 25, seed = 234528)
segment_road <- function(data, num_clusters = 250, n_start = 25, seed = 234528) {
  if (!all(c("x_scaled", "y_scaled", "z_scaled") %in% colnames(data))) {
    stop("Data must contain x_scaled, y_scaled, and z_scaled columns")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Prepare data for clustering
  data2 <- data[, c("x_scaled", "y_scaled", "z_scaled")]
  
  # Perform K-means clustering
  kmeans_result <- kmeans(data2, centers = num_clusters, nstart = n_start)
  
  # Add cluster assignments
  data$Cluster <- kmeans_result$cluster
  data$Cluster_R <- as.numeric(data$Cluster)
  
  return(data)
}

#' Save segmented sections to CSV files
#'
#' Splits clustered data by cluster ID and saves each cluster as a separate CSV file.
#'
#' @param data Data frame with Cluster column
#' @param output_dir Character string, directory to save section files
#' @return Invisible NULL
#'
#' @examples
#' save_sections(data, "results/sections")
save_sections <- function(data, output_dir) {
  if (!"Cluster" %in% colnames(data)) {
    stop("Data must contain a Cluster column")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Split data by cluster
  data_split <- split(data, data$Cluster)
  
  # Save each cluster as a CSV file
  for (cluster_value in names(data_split)) {
    file_name <- file.path(output_dir, paste0("LAS_sec_", cluster_value, ".csv"))
    utils::write.csv(data_split[[cluster_value]], file = file_name, row.names = FALSE)
  }
  
  cat("Saved", length(data_split), "section files to", output_dir, "\n")
  invisible(NULL)
}

#' Complete data processing pipeline
#'
#' Loads LAS data, segments into clusters, and saves sections.
#' This is a convenience function that combines load_las_data, segment_road, and save_sections.
#'
#' @param las_file Character string, path to LAS file
#' @param output_dir Character string, directory to save sections
#' @param num_clusters Integer, number of clusters
#' @param n_start Integer, number of random starts for K-means
#' @param seed Integer, random seed
#' @return Data frame with clustered data
#'
#' @examples
#' data <- process_and_segment_data("data/raw/minasur.las", "results/sections")
process_and_segment_data <- function(las_file, output_dir, 
                                     num_clusters = 250, 
                                     n_start = 25, 
                                     seed = 234528) {
  cat("Loading LAS data...\n")
  data <- load_las_data(las_file)
  
  cat("Segmenting road into", num_clusters, "clusters...\n")
  data <- segment_road(data, num_clusters = num_clusters, 
                       n_start = n_start, seed = seed)
  
  cat("Saving sections...\n")
  save_sections(data, output_dir)
  
  return(data)
}

