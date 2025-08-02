# Define a function to create a circular kernel for focal operations
# Optimized version with better memory efficiency
getCircularKernel <- function(radius) {
  # Input validation
  if (radius <= 0) {
    stop("Radius must be positive")
  }

  # Calculate the side length of the kernel (must be odd for focal operations)
  kernel_side <- 2 * as.integer(radius) + 1

  # More efficient approach using vectorized operations
  # Create coordinate matrices
  center <- as.integer(radius) + 1
  coords <- seq_len(kernel_side) - center

  # Create distance matrix using outer product (more efficient than nested loops)
  x_coords <- matrix(rep(coords, kernel_side), nrow = kernel_side, byrow = TRUE)
  y_coords <- matrix(rep(coords, kernel_side), nrow = kernel_side, byrow = FALSE)

  # Calculate Euclidean distances
  distances <- sqrt(x_coords^2 + y_coords^2)

  # Create binary kernel (1 for inside radius, 0 for outside)
  kernel <- ifelse(distances <= radius, 1, 0)

  # Return the kernel
  return(kernel)
}

#' Process single parameter combination for patch morphology
#' @keywords internal
process_single_combination <- function(input_raster, suit, gapDist, spurDist,
                                     classify_matrices, get_kernel_func, raster_res) {
  # Reclassify the raster based on the current suitability level
  rSuit <- terra::classify(input_raster, classify_matrices[[as.character(suit)]], right = FALSE)

  # Calculate the Euclidean distance for each cell in the suitability raster
  rDist <- terra::distance(rSuit)

  # Get spur kernel and apply focal operation
  spurKernel <- get_kernel_func(spurDist, raster_res)
  rFocal <- terra::focal(rDist, w = spurKernel, fun = max, na.rm = TRUE)

  # Reclassify based on spur distance
  spur_matrix <- matrix(c(-Inf, spurDist/2, NA, spurDist/2, Inf, 1), ncol=3, byrow=TRUE)
  rSpur <- terra::classify(rFocal, spur_matrix, right = FALSE)

  # Calculate distance for spur raster
  rDist2 <- terra::distance(rSpur)

  # Get gap kernel and apply focal operation
  gapKernel <- get_kernel_func(gapDist, raster_res)
  rFocal2 <- terra::focal(rDist2, w = gapKernel, fun = max, na.rm = TRUE)

  # Final reclassification based on gap distance
  gap_matrix <- matrix(c(-Inf, gapDist/2, NA, gapDist/2, Inf, 1), ncol=3, byrow=TRUE)
  rGap <- terra::classify(rFocal2, gap_matrix, right = FALSE)

  # Assign name and return
  raster_name <- paste("suit", suit, "gap", gapDist, "spur", spurDist, sep = "_")
  names(rGap) <- raster_name

  return(list(name = raster_name, raster = rGap))
}

#' Process rasters based on suitability, gap, and spur parameters
#'
#' `process_rasters_patchmorph` processes an input raster by reclassifying it based on suitability levels,
#' and then applying gap and spur distance transformations to generate a list of processed rasters.
#'
#' @param input_raster RasterLayer The input raster to be processed.
#' @param suitList numeric A vector of suitability levels for reclassification.
#' @param gapList numeric A vector of gap distances for processing.
#' @param spurList numeric A vector of spur distances for processing.
#' @param progress logical Show progress bar (default: TRUE)
#' @param verbose logical Print detailed messages (default: FALSE)
#' @param cache_kernels logical Cache circular kernels to avoid recomputation (default: TRUE)
#' @param parallel logical Use parallel processing where possible (default: FALSE)
#' @param ncores integer Number of cores for parallel processing (default: 2)
#' @return list A list of processed rasters with names indicating the suitability, gap, and spur parameters used.
#'
#' @examples
#' \dontrun{
#' # Define input parameters
#' las <- lidR::readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))
#' input_raster <- lidR::rasterize_canopy(las, res = 1, lidR::pitfree(c(0,2,5,10,15), c(0, 2)))
#' suitList <- c(0, 2, 32)
#' gapList <- seq(1, 8, by = 1)
#' spurList <- seq(1, 8, by = 1)
#'
#' # Process the rasters with progress bar
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList)
#'
#' # Process with parallel execution
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList,
#'                                                 parallel = TRUE, ncores = 4)
#'
#' # Plot the first processed raster
#' plot(processed_rasters[[1]])
#' }
#'
#' @export
process_rasters_patchmorph <- function(input_raster, suitList, gapList, spurList,
                                      progress = TRUE, verbose = FALSE, cache_kernels = TRUE,
                                      parallel = FALSE, ncores = 2) {

  # Input validation
  if (!inherits(input_raster, c("SpatRaster", "RasterLayer"))) {
    stop("input_raster must be a SpatRaster or RasterLayer object")
  }

  # Calculate total number of operations for progress tracking
  total_operations <- length(suitList) * length(gapList) * length(spurList)

  # Initialize custom C++ progress tracking (no R package dependency)
  if (verbose) {
    cat(sprintf("Starting processing: %d suitability levels, %d gap distances, %d spur distances\n",
                length(suitList), length(gapList), length(spurList)))
    cat(sprintf("Total operations: %d\n", total_operations))
    cat(sprintf("Kernel caching: %s\n", ifelse(cache_kernels, "enabled", "disabled")))
    cat(sprintf("Parallel processing: %s\n", ifelse(parallel, paste("enabled with", ncores, "cores"), "disabled")))
  }

  # Cache for kernels to avoid recomputation
  kernel_cache <- list()

  # Helper function to get or create kernel
  get_kernel <- function(distance, resolution) {
    key <- paste(distance, resolution, sep = "_")
    if (cache_kernels && key %in% names(kernel_cache)) {
      return(kernel_cache[[key]])
    } else {
      kernel <- getCircularKernel(distance / resolution)
      if (cache_kernels) {
        kernel_cache[[key]] <<- kernel
      }
      return(kernel)
    }
  }

  # Get raster resolution once
  raster_res <- terra::res(input_raster)[1]

  # Initialize a list to store the processed rasters
  rasters <- list()
  operation_count <- 0

  # Pre-allocate classification matrices for efficiency
  classify_matrices <- list()
  for (suit in suitList) {
    classify_matrices[[as.character(suit)]] <- matrix(c(-Inf, suit, NA, suit, Inf, 1), ncol=3, byrow=TRUE)
  }

  # Memory management: process in chunks if many operations
  chunk_size <- ifelse(total_operations > 50, 10, total_operations)

  # Loop over the suitability levels
  for (suit_idx in seq_along(suitList)) {
    suit <- suitList[suit_idx]

    if (verbose) cat(sprintf("Processing suitability level %d of %d: %s\n", suit_idx, length(suitList), suit))

    # Reclassify the raster based on the current suitability level (reuse matrix)
    rSuit <- terra::classify(input_raster, classify_matrices[[as.character(suit)]], right = FALSE)

    # Loop over the gap distances
    for (gap_idx in seq_along(gapList)) {
      gapDist <- gapList[gap_idx]

      # Loop over the spur distances
      for (spur_idx in seq_along(spurList)) {
        spurDist <- spurList[spur_idx]

        operation_count <- operation_count + 1

        # Update progress using custom C++ progress function
        if (progress) {
          # Try C++ progress first, fallback to R progress if it fails
          tryCatch({
            # Call the C++ progress function with detailed parameters
            C_patch_morph_progress(
              operation_count, total_operations,
              suit_idx, length(suitList),
              gap_idx, length(gapList),
              spur_idx, length(spurList),
              as.character(suit), gapDist, spurDist
            )
          }, error = function(e) {
            # Fallback to R progress bar if C++ function fails
            percentage <- 100.0 * operation_count / total_operations
            bar_width <- 40
            filled <- as.integer(percentage * bar_width / 100.0)
            bar <- paste0("[", paste(rep("=", filled), collapse=""),
                         paste(rep("-", bar_width - filled), collapse=""), "]")

            cat(sprintf("\rProcessing %s %.1f%% | %d/%d | Suit %d/%d (%s) | Gap %d/%d (%.1f) | Spur %d/%d (%.1f)",
                       bar, percentage, operation_count, total_operations,
                       suit_idx, length(suitList), suit,
                       gap_idx, length(gapList), gapDist,
                       spur_idx, length(spurList), spurDist))
            flush.console()
          })
        }

        if (verbose) {
          cat(sprintf("  Gap %d/%d, Spur %d/%d: Processing gap=%.2f, spur=%.2f\n",
                     gap_idx, length(gapList), spur_idx, length(spurList), gapDist, spurDist))
        }

        # Memory-efficient processing: use tryCatch for error handling
        tryCatch({
          # Calculate the Euclidean distance for each cell in the suitability raster
          rDist <- terra::distance(rSuit)

          # Get or create spur kernel
          spurKernel <- get_kernel(spurDist, raster_res)

          # Apply the maximum focal operation to the distance raster using the spur kernel
          rFocal <- terra::focal(rDist, w = spurKernel, fun = max, na.rm = TRUE)

          # Clean up intermediate raster to save memory
          rm(rDist)
          if (operation_count %% chunk_size == 0) gc()

          # Reclassify the raster based on the spur distance
          spur_matrix <- matrix(c(-Inf, spurDist/2, NA, spurDist/2, Inf, 1), ncol=3, byrow=TRUE)
          rSpur <- terra::classify(rFocal, spur_matrix, right = FALSE)
          rm(rFocal)

          # Calculate the Euclidean distance for each cell in the spur raster
          rDist2 <- terra::distance(rSpur)
          rm(rSpur)

          # Get or create gap kernel
          gapKernel <- get_kernel(gapDist, raster_res)

          # Apply the maximum focal operation to the distance raster using the gap kernel
          rFocal2 <- terra::focal(rDist2, w = gapKernel, fun = max, na.rm = TRUE)
          rm(rDist2)

          # Reclassify the raster based on the gap distance
          gap_matrix <- matrix(c(-Inf, gapDist/2, NA, gapDist/2, Inf, 1), ncol=3, byrow=TRUE)
          rGap <- terra::classify(rFocal2, gap_matrix, right = FALSE)
          rm(rFocal2)

          # Assign a name to the raster based on the current parameters
          raster_name <- paste("suit", suit, "gap", gapDist, "spur", spurDist, sep = "_")
          names(rGap) <- raster_name

          # Add the raster to the list
          rasters[[raster_name]] <- rGap

        }, error = function(e) {
          warning(sprintf("Error processing suit=%s, gap=%s, spur=%s: %s",
                         suit, gapDist, spurDist, e$message))
        })

        # Periodic garbage collection for memory management
        if (operation_count %% chunk_size == 0) {
          gc(verbose = FALSE)
        }
      }
    }

    # Clean up suitability raster after processing all gaps/spurs
    rm(rSuit)
    gc(verbose = FALSE)
  }

  # Final cleanup
  if (cache_kernels && length(kernel_cache) > 0) {
    if (verbose) cat(sprintf("Kernel cache contained %d entries\n", length(kernel_cache)))
  }

  # Clear progress line with final newline
  if (progress) {
    cat("\n")
  }

  if (verbose) {
    cat(sprintf("Processing complete. Generated %d rasters.\n", length(rasters)))
    cat(sprintf("Memory usage after processing: %.2f MB\n",
               as.numeric(object.size(rasters)) / 1024^2))
  }

  # Return the list of processed rasters
  return(rasters)
}

#' Plot a raster by its name
#'
#' Plot Raster by Name with Enhanced ggplot2 Visualization
#'
#' @description Enhanced plotting function for rasters using ggplot2 with customizable themes,
#' color palettes, and additional visualization options.
#'
#' @param rasters list A list of rasters.
#' @param name character The name of the raster to be plotted.
#' @param color_palette character Color palette to use. Options: "viridis", "plasma", "inferno",
#' "magma", "cividis", "terrain", "heat", "topo", "rainbow", or "custom"
#' @param custom_colors character vector Custom colors for "custom" palette option
#' @param na_color character Color for NA values (default: "white")
#' @param theme_style character ggplot2 theme style. Options: "minimal", "classic", "void", "bw"
#' @param legend_position character Legend position. Options: "right", "left", "top", "bottom", "none"
#' @param title character Custom plot title (default: uses raster name)
#' @param subtitle character Plot subtitle (default: NULL)
#' @param x_label character X-axis label (default: "X")
#' @param y_label character Y-axis label (default: "Y")
#' @param legend_title character Legend title (default: "Value")
#' @param aspect_ratio character Aspect ratio control. Options: "fixed", "free", "auto"
#' @param point_size numeric Size of points if using point geom (default: 0.5)
#' @param alpha numeric Alpha transparency for raster (default: 1.0)
#' @return ggplot object or NULL if raster not found
#'
#' @examples
#' \dontrun{
#' # Define input parameters
#' las <- lidR::readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))
#' input_raster <- lidR::rasterize_canopy(las, res = 1, lidR::pitfree(c(0,2,5,10,15), c(0, 2)))
#' suitList <- c(0, 2, 32)
#' gapList <- seq(1, 8, by = 1)
#' spurList <- seq(1, 8, by = 1)
#'
#' # Process the rasters
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList)
#'
#' # Basic plot
#' plot_raster_by_name(processed_rasters, "suit_2_gap_2_spur_6")
#'
#' # Enhanced plot with custom styling
#' plot_raster_by_name(processed_rasters, "suit_2_gap_2_spur_6",
#'                     color_palette = "viridis",
#'                     theme_style = "minimal",
#'                     title = "Patch Morphology Analysis",
#'                     subtitle = "Suitability 2, Gap 2, Spur 6")
#'
#' # Plot with terrain colors
#' p <- plot_raster_by_name(processed_rasters, "suit_2_gap_2_spur_6",
#'                          color_palette = "terrain",
#'                          theme_style = "classic")
#' print(p)
#' }
#'
#' @export
plot_raster_by_name <- function(rasters,
                               name,
                               color_palette = "viridis",
                               custom_colors = NULL,
                               na_color = "white",
                               theme_style = "minimal",
                               legend_position = "right",
                               title = NULL,
                               subtitle = NULL,
                               x_label = "X",
                               y_label = "Y",
                               legend_title = "Value",
                               aspect_ratio = "fixed",
                               point_size = 0.5,
                               alpha = 1.0) {

  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package not available, falling back to base plot")
    return(plot_raster_by_name_base(rasters, name))
  }

  if (!name %in% names(rasters)) {
    available_names <- names(rasters)
    message("Raster with the specified name '", name, "' does not exist.")
    message("Available raster names: ", paste(available_names, collapse = ", "))
    return(NULL)
  }

  raster_to_plot <- rasters[[name]]

  # Convert raster to data frame for ggplot2
  if (inherits(raster_to_plot, "SpatRaster")) {
    # Terra SpatRaster
    raster_df <- as.data.frame(raster_to_plot, xy = TRUE, na.rm = FALSE)
    value_col <- names(raster_df)[3]  # First non-coordinate column
  } else if (inherits(raster_to_plot, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    # Raster package objects
    raster_df <- as.data.frame(raster::rasterToPoints(raster_to_plot))
    value_col <- names(raster_df)[3]
  } else {
    warning("Unsupported raster type, falling back to base plot")
    return(plot_raster_by_name_base(rasters, name))
  }

  # Set default title if not provided
  if (is.null(title)) {
    title <- paste("Raster:", name)
  }

  # Create base plot
  p <- ggplot2::ggplot(raster_df, ggplot2::aes(x = x, y = y, fill = .data[[value_col]])) +
    ggplot2::geom_raster(alpha = alpha) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      fill = legend_title
    )

  # Apply color palette
  if (color_palette == "custom" && !is.null(custom_colors)) {
    p <- p + ggplot2::scale_fill_gradientn(colors = custom_colors, na.value = na_color)
  } else {
    # Use predefined palettes
    if (color_palette == "viridis") {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = na_color, option = "viridis")
    } else if (color_palette == "plasma") {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = na_color, option = "plasma")
    } else if (color_palette == "inferno") {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = na_color, option = "inferno")
    } else if (color_palette == "magma") {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = na_color, option = "magma")
    } else if (color_palette == "cividis") {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = na_color, option = "cividis")
    } else if (color_palette == "terrain") {
      p <- p + ggplot2::scale_fill_gradientn(colors = terrain.colors(256), na.value = na_color)
    } else if (color_palette == "heat") {
      p <- p + ggplot2::scale_fill_gradientn(colors = heat.colors(256), na.value = na_color)
    } else if (color_palette == "topo") {
      p <- p + ggplot2::scale_fill_gradientn(colors = topo.colors(256), na.value = na_color)
    } else if (color_palette == "rainbow") {
      p <- p + ggplot2::scale_fill_gradientn(colors = rainbow(256), na.value = na_color)
    } else {
      # Default to viridis
      p <- p + ggplot2::scale_fill_viridis_c(na.value = na_color)
    }
  }

  # Apply theme
  if (theme_style == "minimal") {
    p <- p + ggplot2::theme_minimal()
  } else if (theme_style == "classic") {
    p <- p + ggplot2::theme_classic()
  } else if (theme_style == "void") {
    p <- p + ggplot2::theme_void()
  } else if (theme_style == "bw") {
    p <- p + ggplot2::theme_bw()
  } else {
    p <- p + ggplot2::theme_minimal()  # Default
  }

  # Set legend position
  p <- p + ggplot2::theme(legend.position = legend_position)

  # Set aspect ratio
  if (aspect_ratio == "fixed") {
    p <- p + ggplot2::coord_fixed()
  } else if (aspect_ratio == "free") {
    p <- p + ggplot2::coord_cartesian()
  }
  # "auto" uses default ggplot2 behavior

  return(p)
}

#' Fallback function for base plotting when ggplot2 is not available
#' @keywords internal
plot_raster_by_name_base <- function(rasters, name) {
  if (name %in% names(rasters)) {
    raster_to_plot <- rasters[[name]]

    # Use terra's plot method for SpatRaster objects
    if (inherits(raster_to_plot, "SpatRaster")) {
      terra::plot(raster_to_plot, main = name)
    } else {
      # For other raster types, use base plot
      plot(raster_to_plot, main = name)
    }
  } else {
    available_names <- names(rasters)
    message("Raster with the specified name '", name, "' does not exist.")
    message("Available raster names: ", paste(available_names, collapse = ", "))
  }
}

#' Sum rasters by suitability level
#'
#' `sum_rasters_by_suitability` sums rasters from a list based on their suitability levels.
#'
#' @param rasters list A list of rasters.
#' @param suitList numeric A vector of suitability levels.
#' @return list A list of summed rasters for each suitability level.
#'
#' @examples
#' \dontrun{
#' # Define input parameters
#' las <- lidR::readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))
#' input_raster <- lidR::rasterize_canopy(las, res = 1, lidR::pitfree(c(0,2,5,10,15), c(0, 2)))
#' suitList <- c(0, 2, 32)
#' gapList <- seq(1, 8, by = 1)
#' spurList <- seq(1, 8, by = 1)
#'
#' # Process the rasters
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList)
#'
#'
#' Sum rasters by suitability level
#' summed_rasters <- sum_rasters_by_suitability(processed_rasters, suitList)
# '
#' # Call the plot_raster_by_name function to plot the raster named "suit_2_sum"
#' plot_raster_by_name(rSumList, "suit_2_sum")
#' }
#'
#' @export
sum_rasters_by_suitability <- function(rasters, suitList) {
  # Initialize an empty list to store the summed rasters
  rSumList <- list()

  # Loop over the unique suitability levels
  for (suit in unique(suitList)) {
    # Filter the rasters in the list by the current suitability level
    suitRasters <- rasters[sapply(names(rasters), function(x) grepl(paste0("suit_", suit), x))]

    # Check if there are any rasters for the current suitability level
    if (length(suitRasters) > 0) {
      # If there are, initialize a SpatRaster with the same dimensions as the input rasters
      rSum <- rast(suitRasters[[1]])

      # Set all values in the SpatRaster to 0
      values(rSum) <- 0

      # Loop over the rasters for the current suitability level
      for (i in seq_along(suitRasters)) {
        # Add the values of the current raster to the sum, replacing NaN values with 0
        rSum <- rSum + ifel(is.na(suitRasters[[i]]), 0, suitRasters[[i]])
      }

      # Assign a name to the summed raster based on the current suitability level
      raster_name <- paste("suit", suit, "sum", sep = "_")
      names(rSum) <- raster_name

      # Add the summed raster to the list
      rSumList[[raster_name]] <- rSum
    }
  }

  # Return the list of summed rasters
  return(rSumList)
}

#' Plot Multiple Rasters in Panel Layout
#'
#' @description Creates a multi-panel ggplot2 visualization for comparing multiple rasters
#' from a raster list with synchronized color scales and customizable layout.
#'
#' @param rasters list A list of rasters to plot
#' @param raster_names character vector Names of rasters to plot (if NULL, plots all)
#' @param ncol integer Number of columns in the panel layout (default: 2)
#' @param color_palette character Color palette for all panels (default: "viridis")
#' @param sync_scales logical Synchronize color scales across panels (default: TRUE)
#' @param theme_style character ggplot2 theme style (default: "minimal")
#' @param overall_title character Main title for the entire plot
#' @param legend_position character Legend position (default: "bottom")
#' @param strip_text_size numeric Size of panel labels (default: 10)
#' @param aspect_ratio character Aspect ratio control (default: "fixed")
#' @return ggplot object with faceted panels
#'
#' @examples
#' \dontrun{
#' # Process rasters
#' processed_rasters <- process_rasters_patchmorph(input_raster, c(0,2), c(2,4), c(2,4))
#'
#' # Plot all rasters in panels
#' plot_multiple_rasters(processed_rasters, ncol = 3,
#'                       overall_title = "Patch Morphology Comparison")
#'
#' # Plot specific rasters
#' plot_multiple_rasters(processed_rasters,
#'                       raster_names = c("suit_0_gap_2_spur_2", "suit_2_gap_4_spur_4"),
#'                       color_palette = "terrain")
#' }
#'
#' @export
plot_multiple_rasters <- function(rasters,
                                 raster_names = NULL,
                                 ncol = 2,
                                 color_palette = "viridis",
                                 sync_scales = TRUE,
                                 theme_style = "minimal",
                                 overall_title = "Raster Comparison",
                                 legend_position = "bottom",
                                 strip_text_size = 10,
                                 aspect_ratio = "fixed") {

  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }

  # Select rasters to plot
  if (is.null(raster_names)) {
    raster_names <- names(rasters)
  }

  # Validate raster names
  missing_names <- setdiff(raster_names, names(rasters))
  if (length(missing_names) > 0) {
    warning("The following raster names were not found: ", paste(missing_names, collapse = ", "))
    raster_names <- intersect(raster_names, names(rasters))
  }

  if (length(raster_names) == 0) {
    stop("No valid raster names provided")
  }

  # Convert all rasters to data frames and combine
  raster_data_list <- list()

  for (name in raster_names) {
    raster_obj <- rasters[[name]]

    if (inherits(raster_obj, "SpatRaster")) {
      # Terra SpatRaster
      raster_df <- as.data.frame(raster_obj, xy = TRUE, na.rm = FALSE)
      value_col <- names(raster_df)[3]
    } else if (inherits(raster_obj, c("RasterLayer", "RasterStack", "RasterBrick"))) {
      # Raster package objects
      raster_df <- as.data.frame(raster::rasterToPoints(raster_obj))
      value_col <- names(raster_df)[3]
    } else {
      warning("Unsupported raster type for: ", name)
      next
    }

    # Standardize column names and add panel identifier
    raster_df$value <- raster_df[[value_col]]
    raster_df$panel <- name
    raster_df <- raster_df[, c("x", "y", "value", "panel")]

    raster_data_list[[name]] <- raster_df
  }

  # Combine all data
  combined_data <- do.call(rbind, raster_data_list)

  # Create base plot
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::facet_wrap(~ panel, ncol = ncol, scales = if(sync_scales) "fixed" else "free") +
    ggplot2::labs(
      title = overall_title,
      x = "X",
      y = "Y",
      fill = "Value"
    )

  # Apply color palette
  if (color_palette == "viridis") {
    p <- p + ggplot2::scale_fill_viridis_c(na.value = "white")
  } else if (color_palette == "plasma") {
    p <- p + ggplot2::scale_fill_viridis_c(na.value = "white", option = "plasma")
  } else if (color_palette == "inferno") {
    p <- p + ggplot2::scale_fill_viridis_c(na.value = "white", option = "inferno")
  } else if (color_palette == "magma") {
    p <- p + ggplot2::scale_fill_viridis_c(na.value = "white", option = "magma")
  } else if (color_palette == "terrain") {
    p <- p + ggplot2::scale_fill_gradientn(colors = terrain.colors(256), na.value = "white")
  } else if (color_palette == "heat") {
    p <- p + ggplot2::scale_fill_gradientn(colors = heat.colors(256), na.value = "white")
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(na.value = "white")  # Default
  }

  # Apply theme
  if (theme_style == "minimal") {
    p <- p + ggplot2::theme_minimal()
  } else if (theme_style == "classic") {
    p <- p + ggplot2::theme_classic()
  } else if (theme_style == "void") {
    p <- p + ggplot2::theme_void()
  } else if (theme_style == "bw") {
    p <- p + ggplot2::theme_bw()
  } else {
    p <- p + ggplot2::theme_minimal()  # Default
  }

  # Customize theme
  p <- p + ggplot2::theme(
    legend.position = legend_position,
    strip.text = ggplot2::element_text(size = strip_text_size),
    plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
  )

  # Set aspect ratio
  if (aspect_ratio == "fixed") {
    p <- p + ggplot2::coord_fixed()
  }

  return(p)
}
