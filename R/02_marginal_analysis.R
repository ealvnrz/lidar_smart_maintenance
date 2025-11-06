# =============================================================================
# Marginal Analysis Functions
# =============================================================================
# Functions for fitting Weibull marginals and diagnostic checks

#' Fit Weibull distribution to intensity data (ignoring spatial correlation)
#'
#' Fits a Weibull distribution to intensity data using maximum likelihood,
#' treating observations as independent (marginal fit).
#'
#' @param intensity_data Numeric vector of intensity values (must be > 0)
#' @return List with fitted parameters and fit object
#'
#' @examples
#' fit <- fit_weibull_marginal(data$intensity_scaled)
fit_weibull_marginal <- function(intensity_data) {
  # Filter out zero and negative values
  intensity_positive <- intensity_data[intensity_data > 0]
  
  if (length(intensity_positive) == 0) {
    stop("No positive intensity values found for Weibull fitting")
  }
  
  # Fit Weibull distribution using MASS::fitdistr
  fit_weibull <- MASS::fitdistr(intensity_positive, "weibull")
  
  # Extract parameters
  shape <- fit_weibull$estimate["shape"]
  scale <- fit_weibull$estimate["scale"]
  
  result <- list(
    shape = as.numeric(shape),
    scale = as.numeric(scale),
    fit_object = fit_weibull,
    n_obs = length(intensity_positive)
  )
  
  return(result)
}

#' Generate marginal diagnostic plots
#'
#' Creates QQ plot and density comparison for Weibull marginal fit.
#'
#' @param fit Fitted Weibull object from fit_weibull_marginal
#' @param data Numeric vector of intensity data
#' @param title Character string, plot title
#' @return List of plot objects (histogram and QQ plot)
#'
#' @examples
#' plots <- diagnose_marginals(fit, data$intensity_scaled, "Cluster 2")
diagnose_marginals <- function(fit, data, title = "Weibull Marginal Fit") {
  # Filter positive values
  data_positive <- data[data > 0]
  
  # Create histogram with Weibull overlay
  hist_data <- data.frame(intensity = data_positive)
  
  p1 <- ggplot2::ggplot(hist_data, ggplot2::aes(x = intensity)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..), 
                            bins = 15, 
                            fill = "lightblue", 
                            color = "darkblue") +
    ggplot2::stat_function(fun = dweibull, 
                          args = list(shape = fit$shape, scale = fit$scale), 
                          color = "red", 
                          size = 1.5) +
    ggplot2::labs(title = title, 
                 x = "Intensity", 
                 y = "Density")
  
  # QQ plot using GeoModels
  # Note: This requires creating a GeoFit object for GeoQQ
  # For now, we'll create a simple QQ plot manually
  theoretical_quantiles <- qweibull(ppoints(length(data_positive)), 
                                   shape = fit$shape, 
                                   scale = fit$scale)
  sample_quantiles <- sort(data_positive)
  
  qq_data <- data.frame(theoretical = theoretical_quantiles,
                       sample = sample_quantiles)
  
  p2 <- ggplot2::ggplot(qq_data, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
    ggplot2::labs(title = paste("QQ Plot:", title),
                 x = "Theoretical Quantiles",
                 y = "Sample Quantiles")
  
  return(list(histogram = p1, qqplot = p2))
}

#' Check marginal fit for a specific cluster
#'
#' Loads section data and validates Weibull marginal fit.
#' Note: Requires utils.R to be loaded for load_section_data().
#'
#' @param sections_dir Character string, directory containing section files
#' @param cluster_id Integer, cluster ID to check
#' @param plot Logical, whether to generate diagnostic plots
#' @return List with fit results and optionally plots
#'
#' @examples
#' result <- check_marginal_fit("results/sections", cluster_id = 2, plot = TRUE)
check_marginal_fit <- function(sections_dir, cluster_id, plot = FALSE) {
  # Load section data (requires utils.R)
  if (!exists("load_section_data")) {
    stop("load_section_data() not found. Please source R/utils.R first.")
  }
  data <- load_section_data(sections_dir, cluster_id)
  
  # Filter positive intensity values
  intensity_positive <- data$intensity_scaled[data$intensity_scaled > 0]
  
  if (length(intensity_positive) == 0) {
    warning(paste("No positive intensity values in cluster", cluster_id))
    return(NULL)
  }
  
  # Fit Weibull
  fit <- fit_weibull_marginal(data$intensity_scaled)
  
  result <- list(
    cluster_id = cluster_id,
    fit = fit,
    n_obs = length(intensity_positive),
    mean_intensity = mean(intensity_positive),
    median_intensity = median(intensity_positive)
  )
  
  # Generate plots if requested
  if (plot) {
    title <- paste("Cluster", cluster_id, "Weibull Marginal")
    plots <- diagnose_marginals(fit, data$intensity_scaled, title)
    result$plots <- plots
  }
  
  return(result)
}

#' Fit Weibull marginals using GeoModels (for spatial-aware diagnostics)
#'
#' Uses GeoModels package to fit Weibull marginals with spatial structure.
#' This is useful for more advanced diagnostics that account for spatial dependence.
#'
#' @param intensity_data Numeric vector of intensity values
#' @param coords Matrix or data frame with coordinates (n x 2)
#' @param start_params List with starting parameters (mean, shape)
#' @return GeoFit object from GeoModels
#'
#' @examples
#' coords <- cbind(data$x, data$y)
#' fit <- fit_weibull_marginal_geospatial(data$intensity_scaled, coords)
fit_weibull_marginal_geospatial <- function(intensity_data, coords, 
                                            start_params = list(mean = 0, shape = 7)) {
  if (nrow(coords) != length(intensity_data)) {
    stop("Coordinates and intensity data must have same length")
  }
  
  # Fit using GeoModels with Independence likelihood
  fit <- GeoModels::GeoFit(
    data = intensity_data,
    coordx = coords,
    model = "Weibull",
    start = start_params,
    likelihood = "Marginal",
    type = "Independence"
  )
  
  return(fit)
}

#' Compute residuals from GeoFit object
#'
#' @param fit GeoFit object from GeoModels
#' @return GeoResiduals object
#'
#' @examples
#' residuals <- compute_geospatial_residuals(fit)
compute_geospatial_residuals <- function(fit) {
  residuals <- GeoModels::GeoResiduals(fit)
  return(residuals)
}

