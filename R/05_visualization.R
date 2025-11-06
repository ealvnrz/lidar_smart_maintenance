# =============================================================================
# Visualization Functions
# =============================================================================
# Functions for creating plots and visualizations

#' Plot intensity histogram with Weibull overlay
#'
#' Creates a histogram of intensity values with fitted Weibull distribution overlay.
#'
#' @param data Numeric vector of intensity values
#' @param fit List with Weibull fit parameters (shape, scale) or fit object
#' @param title Character string, plot title
#' @param bins Integer, number of histogram bins (default: 15)
#' @param fill_color Character string, fill color for histogram (default: "lightblue")
#' @param line_color Character string, color for Weibull curve (default: "red")
#' @return ggplot object
#'
#' @examples
#' fit <- fit_weibull_marginal(data$intensity_scaled)
#' p <- plot_intensity_histogram(data$intensity_scaled, fit, "Cluster 2")
plot_intensity_histogram <- function(data, fit, title = "Intensity Histogram",
                                     bins = 15,
                                     fill_color = "lightblue",
                                     line_color = "red") {
  # Filter positive values
  data_positive <- data[data > 0]
  
  if (length(data_positive) == 0) {
    stop("No positive intensity values to plot")
  }
  
  # Extract shape and scale from fit
  if (is.list(fit)) {
    if ("shape" %in% names(fit) && "scale" %in% names(fit)) {
      shape <- fit$shape
      scale <- fit$scale
    } else if ("fit_object" %in% names(fit)) {
      # Extract from MASS fitdistr object
      shape <- fit$fit_object$estimate["shape"]
      scale <- fit$fit_object$estimate["scale"]
    } else {
      stop("Could not extract shape and scale from fit object")
    }
  } else {
    stop("fit must be a list with shape and scale parameters")
  }
  
  # Create plot
  p <- ggplot2::ggplot(data.frame(intensity = data_positive), 
                      ggplot2::aes(x = intensity)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..), 
                            bins = bins, 
                            fill = fill_color, 
                            color = "darkblue") +
    ggplot2::stat_function(fun = dweibull, 
                          args = list(shape = shape, scale = scale), 
                          color = line_color, 
                          size = 1.5) +
    ggplot2::labs(title = title, 
                 x = "Intensity", 
                 y = "Density") +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Plot 3D intensity visualization
#'
#' Creates an interactive 3D plot of intensity values using plotly.
#'
#' @param coords Data frame or matrix with x, y, z coordinates
#' @param intensity Numeric vector of intensity values
#' @param title Character string, plot title
#' @param sample_size Numeric, proportion of points to sample (default: 0.2)
#' @param marker_size Numeric, size of markers (default: 2)
#' @param colorscale Character string, plotly colorscale (default: "Viridis")
#' @return plotly object
#'
#' @examples
#' p <- plot_3d_intensity(data.frame(x=x, y=y, z=z), intensity, "LiDAR Intensity")
plot_3d_intensity <- function(coords, intensity, title = "3D Intensity Plot",
                              sample_size = 0.2,
                              marker_size = 2,
                              colorscale = "Viridis") {
  
  # Prepare data
  if (is.data.frame(coords)) {
    x <- coords$x
    y <- coords$y
    z <- coords$z
  } else if (is.matrix(coords) && ncol(coords) >= 3) {
    x <- coords[, 1]
    y <- coords[, 2]
    z <- coords[, 3]
  } else {
    stop("coords must be a data frame with x, y, z columns or a matrix with 3 columns")
  }
  
  # Sample data if needed
  n_points <- length(x)
  if (sample_size < 1 && n_points > 1000) {
    sample_indices <- sample(1:n_points, size = sample_size * n_points)
    x <- x[sample_indices]
    y <- y[sample_indices]
    z <- z[sample_indices]
    intensity <- intensity[sample_indices]
  }
  
  # Create plot
  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      x = ~x,
      y = ~y,
      z = ~z,
      color = ~intensity,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = marker_size,
        colorscale = colorscale,
        showscale = TRUE
      ),
      name = "Intensity"
    ) %>%
    plotly::layout(
      title = title,
      scene = list(
        xaxis = list(title = "X"),
        yaxis = list(title = "Y"),
        zaxis = list(title = "Z")
      )
    )
  
  return(p)
}

#' Plot parameter boxplots
#'
#' Creates boxplots for WRF parameters (mean, scale, shape) across clusters.
#'
#' @param param_df Data frame with parameter estimates (must contain mean, scale, shape columns)
#' @param param_name Character string, parameter name to plot (default: "scale")
#' @param title Character string, plot title (default: auto-generated)
#' @param fill_color Character string, fill color (default: "lightcoral")
#' @return ggplot object
#'
#' @examples
#' p <- plot_parameter_boxplots(param_df, "scale")
plot_parameter_boxplots <- function(param_df, param_name = "scale",
                                   title = NULL,
                                   fill_color = "lightcoral") {
  
  if (!param_name %in% colnames(param_df)) {
    stop(paste("Parameter", param_name, "not found in data frame"))
  }
  
  if (is.null(title)) {
    title <- paste("'", toupper(param_name), "' parameter boxplot", sep = "")
  }
  
  p <- ggplot2::ggplot(param_df, ggplot2::aes(y = .data[[param_name]])) +
    ggplot2::geom_boxplot(fill = fill_color, color = "black") +
    ggplot2::labs(title = title, x = NULL, y = toupper(param_name)) +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Plot semivariogram
#'
#' Creates a semivariogram plot comparing empirical and model-fitted semivariograms.
#'
#' @param fit GeoFit object from GeoModels
#' @param coords Matrix or data frame with coordinates (n x 2)
#' @param maxdist Numeric, maximum distance for semivariogram (default: max distance / 3)
#' @param title Character string, plot title
#' @return Plot object (from GeoModels)
#'
#' @examples
#' vario <- plot_semivariogram(fit, coords, maxdist = 100)
plot_semivariogram <- function(fit, coords, maxdist = NULL, title = "Semivariogram") {
  
  if (is.null(maxdist)) {
    maxdist <- max(dist(coords)) / 3
  }
  
  # Compute residuals
  residuals <- GeoModels::GeoResiduals(fit)
  
  # Compute empirical semivariogram
  vario <- GeoModels::GeoVariogram(
    data = residuals$data,
    coordx = coords,
    maxdist = maxdist
  )
  
  # Plot covariogram
  p <- GeoModels::GeoCovariogram(
    residuals,
    show.vario = TRUE,
    vario = vario,
    pch = 20
  )
  
  return(p)
}

#' Plot flagged clusters map
#'
#' Creates a spatial map showing which clusters are flagged as degraded.
#'
#' @param cluster_data Data frame with cluster data (must contain x, y, Cluster columns)
#' @param flagged_ids Numeric vector, cluster IDs that are flagged
#' @param indicator_name Character string, name of indicator used for flagging
#' @param title Character string, plot title
#' @param sample_size Numeric, proportion of points to sample (default: 0.1)
#' @return ggplot object
#'
#' @examples
#' p <- plot_flagged_clusters(all_clusters, flagged_ids = c(15, 29, 31), "ELI")
plot_flagged_clusters <- function(cluster_data,
                                 flagged_ids,
                                 indicator_name = "Indicator",
                                 title = NULL,
                                 sample_size = 0.1) {
  
  if (!all(c("x", "y", "Cluster") %in% colnames(cluster_data))) {
    stop("cluster_data must contain x, y, and Cluster columns")
  }
  
  # Add flagged status
  cluster_data$flagged <- cluster_data$Cluster %in% flagged_ids
  
  # Sample data for plotting
  if (sample_size < 1) {
    n_points <- nrow(cluster_data)
    sample_indices <- sample(1:n_points, size = sample_size * n_points)
    cluster_data <- cluster_data[sample_indices, ]
  }
  
  if (is.null(title)) {
    title <- paste("Flagged clusters by", indicator_name)
  }
  
  # Create plot
  p <- ggplot2::ggplot(cluster_data, ggplot2::aes(x = x, y = y, color = flagged)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "gray70", "TRUE" = "red"),
      labels = c("FALSE" = "Normal", "TRUE" = "Flagged"),
      name = "Status"
    ) +
    ggplot2::labs(
      title = title,
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom"
    )
  
  return(p)
}

#' Create multi-panel parameter boxplots
#'
#' Creates boxplots for all three parameters (mean, scale, shape) in a single plot.
#'
#' @param param_df Data frame with parameter estimates
#' @return ggplot object with three panels
#'
#' @examples
#' p <- plot_all_parameter_boxplots(param_df)
plot_all_parameter_boxplots <- function(param_df) {
  
  # Reshape data for plotting
  param_long <- tidyr::pivot_longer(
    param_df,
    cols = c(mean, scale, shape),
    names_to = "parameter",
    values_to = "value"
  )
  
  # Create plot
  p <- ggplot2::ggplot(param_long, ggplot2::aes(x = parameter, y = value)) +
    ggplot2::geom_boxplot(fill = "lightblue", color = "black") +
    ggplot2::labs(
      title = "Parameter dispersion across segments",
      x = "Parameter",
      y = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~parameter, scales = "free")
  
  return(p)
}

#' Save plot to file
#'
#' Convenience function to save plots with consistent settings.
#'
#' @param plot Plot object (ggplot or plotly)
#' @param file_path Character string, path to save file
#' @param width Numeric, plot width in inches (default: 10)
#' @param height Numeric, plot height in inches (default: 6)
#' @param dpi Numeric, resolution (default: 300)
#' @return Invisible NULL
#'
#' @examples
#' save_plot(p, "figures/results/histogram.png")
save_plot <- function(plot, file_path, width = 10, height = 6, dpi = 300) {
  
  # Create directory if needed
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save based on plot type
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(file_path, plot, width = width, height = height, dpi = dpi)
  } else if (inherits(plot, "plotly")) {
    htmlwidgets::saveWidget(plot, file_path)
  } else {
    stop("Unsupported plot type")
  }
  
  cat("Saved plot to", file_path, "\n")
  invisible(NULL)
}

