#' Comprehensive Tree Detection with Two-Phase RANSAC Filtering
#'
#' @description This function provides a complete tree detection workflow using the optimized
#' get_raster_eigen_treelocs_optimized function with C_geometric_features_simple internally.
#' It performs initial tree detection, fits cylinders using RANSAC, and then applies a second
#' pass for trees with poor RANSAC fits, returning an sf object of tree locations.
#' @name comprehensive_tree_detection
#'
#' @param las LAS object containing the normalized point cloud
#' @param neigh_sizes Neighborhood sizes for analysis c(eigen_radius, small_density_radius, large_density_radius)
#' @param res Grid resolution for raster operations
#' @param grid_slice_min Minimum height for processing slice
#' @param grid_slice_max Maximum height for processing slice
#' @param dens_threshold Density threshold for initial tree detection
#' @param vert_threshold Verticality threshold for initial tree detection
#' @param minimum_area Minimum area for tree candidate regions
#' @param pt_spacing Point spacing for downsampling
#' @param SDvert Standard deviation threshold for verticality filtering
#' @param ransac_error_threshold Maximum acceptable RANSAC error for cylinder fitting
#' @param cylinder_search_radius Radius for extracting points around tree locations for cylinder fitting
#' @param min_cylinder_points Minimum points required for cylinder fitting
#' @param n_ransac_samples Number of RANSAC samples for cylinder fitting
#' @param ransac_confidence RANSAC confidence level
#' @param ransac_inliers RANSAC inlier percentage
#' @param max_angle Maximum angle for cylinder fitting (degrees)
#' @param n_best Number of best cylinders to consider
#' @param second_pass_vert_threshold Higher verticality threshold for second pass
#' @param second_pass_aniso_threshold Anisotropy threshold for second pass
#' @param ncpu Number of CPU cores to use
#' @param verbose Print progress messages
#'
#' @return sf object with tree locations and attributes including TreeID, X, Y, Z, Radius, Error, and Quality
#'
#' @examples
#' \dontrun{
#' # Load example data
#' LASfile <- system.file("extdata", "DensePatchA.laz", package="spanner")
#' las <- readLAS(LASfile, select = "xyzc")
#'
#' # Detect trees with comprehensive two-phase approach
#' trees <- comprehensive_tree_detection(
#'   las = las,
#'   neigh_sizes = c(0.5, 0.25, 0.75),
#'   res = 0.05,
#'   ransac_error_threshold = 0.1,
#'   verbose = TRUE
#' )
#'
#' # Plot results
#' plot(las, color = "Z", backend = "lidRviewer")
#' plot(trees$geometry, add = TRUE, col = "red", pch = 16)
#' }
#'
#' @export
comprehensive_tree_detection <- function(las,
                                        neigh_sizes = c(0.5, 0.25, 0.75),
                                        res = 0.05,
                                        grid_slice_min = 1.2,
                                        grid_slice_max = 1.8,
                                        dens_threshold = 0.15,
                                        vert_threshold = 0.6,
                                        minimum_area = 0.1,
                                        pt_spacing = 0.025,
                                        SDvert = 0.15,
                                        ransac_error_threshold = 0.15,
                                        cylinder_search_radius = 0.4,
                                        min_cylinder_points = 15,
                                        n_ransac_samples = 10,
                                        ransac_confidence = 0.95,
                                        ransac_inliers = 0.8,
                                        max_angle = 20,
                                        n_best = 20,
                                        second_pass_vert_threshold = 0.75,
                                        second_pass_aniso_threshold = 0.7,
                                        ncpu = parallel::detectCores(logical = FALSE),
                                        verbose = TRUE) {

  if (verbose) message("=== Comprehensive Tree Detection with Two-Phase RANSAC ===")

  # Phase 1: Initial tree detection using optimized raster-based approach
  if (verbose) message("Phase 1: Initial tree detection...")

  initial_trees <- get_raster_eigen_treelocs_optimized(
    las = las,
    neigh_sizes = neigh_sizes,
    res = res,
    grid_slice_min = grid_slice_min,
    grid_slice_max = grid_slice_max,
    dens_threshold = dens_threshold,
    vert_threshold = vert_threshold,
    minimum_area = minimum_area,
    pt_spacing = pt_spacing,
    SDvert = SDvert,
    ncpu = ncpu
  )

  if (nrow(initial_trees) == 0) {
    warning("No trees detected in initial phase")
    return(create_empty_sf_result())
  }

  if (verbose) message(sprintf("Phase 1 complete: %d tree candidates detected", nrow(initial_trees)))

  # Phase 2: Cylinder fitting and error assessment
  if (verbose) message("Phase 2: RANSAC cylinder fitting...")

  # Create processing slice for cylinder fitting
  slice_extra <- lidR::filter_poi(las, Z >= grid_slice_min, Z <= grid_slice_max)

  # Compute geometric features for the slice if not already done
  if (!any(c("verticality", "anisotropy") %in% names(slice_extra@data))) {
    if (verbose) message("Computing geometric features for cylinder fitting...")
    geom_features <- spanner:::C_geometric_features_simple(
      slice_extra,
      radius = neigh_sizes[1],
      max_neighbors = 50,
      ncpu = ncpu
    )
    slice_extra <- lidR::add_lasattribute(slice_extra, geom_features$verticality, "verticality", "verticality")
    slice_extra <- lidR::add_lasattribute(slice_extra, geom_features$anisotropy, "anisotropy", "anisotropy")
    slice_extra@data$verticality[is.na(slice_extra@data$verticality)] <- 0.5
    slice_extra@data$anisotropy[is.na(slice_extra@data$anisotropy)] <- 0.5
  }

  # Fit cylinders to each tree candidate
  tree_results <- list()
  good_fits <- c()
  poor_fits <- c()

  for (i in 1:nrow(initial_trees)) {
    tree_x <- initial_trees$X[i]
    tree_y <- initial_trees$Y[i]

    # Extract points around tree location
    tree_points <- slice_extra@data[
      sqrt((slice_extra@data$X - tree_x)^2 + (slice_extra@data$Y - tree_y)^2) <= cylinder_search_radius,
    ]

    if (nrow(tree_points) >= min_cylinder_points) {
      # Fit cylinder using RANSAC
      cylinder_fit <- tryCatch({
        spanner::cylinderFit(
          lidR::LAS(tree_points, header = slice_extra@header, crs = lidR::st_crs(slice_extra)),
          method = "ransac",
          n = n_ransac_samples,
          inliers = ransac_inliers,
          conf = ransac_confidence,
          max_angle = max_angle,
          n_best = n_best
        )
      }, error = function(e) NULL)

      if (!is.null(cylinder_fit) && !is.na(cylinder_fit$err)) {
        # Store results
        tree_result <- data.frame(
          TreeID = i,
          X = cylinder_fit$px,
          Y = cylinder_fit$py,
          Z = cylinder_fit$pz,
          Radius = cylinder_fit$radius,
          Error = cylinder_fit$err,
          PointCount = nrow(tree_points),
          Quality = if (cylinder_fit$err <= ransac_error_threshold) "Good" else "Poor",
          stringsAsFactors = FALSE
        )
        tree_results[[i]] <- tree_result

        # Classify fit quality
        if (cylinder_fit$err <= ransac_error_threshold) {
          good_fits <- c(good_fits, i)
        } else {
          poor_fits <- c(poor_fits, i)
        }
      }
    }
  }

  if (verbose) {
    message(sprintf("Cylinder fitting complete: %d good fits, %d poor fits",
                   length(good_fits), length(poor_fits)))
  }

  # Phase 3: Second pass for poor fits using enhanced geometric filtering
  if (length(poor_fits) > 0) {
    if (verbose) message("Phase 3: Enhanced processing for poor fits...")

    for (tree_idx in poor_fits) {
      tree_x <- initial_trees$X[tree_idx]
      tree_y <- initial_trees$Y[tree_idx]

      # Extract larger region for enhanced processing
      enhanced_points <- slice_extra@data[
        sqrt((slice_extra@data$X - tree_x)^2 + (slice_extra@data$Y - tree_y)^2) <= cylinder_search_radius * 1.5,
      ]

      if (nrow(enhanced_points) >= min_cylinder_points) {
        # Apply stricter geometric filtering
        filtered_points <- enhanced_points[
          enhanced_points$verticality >= second_pass_vert_threshold &
          enhanced_points$anisotropy >= second_pass_aniso_threshold,
        ]

        if (nrow(filtered_points) >= min_cylinder_points) {
          # Refit cylinder with filtered points
          enhanced_fit <- tryCatch({
            spanner::cylinderFit(
              lidR::LAS(filtered_points, header = slice_extra@header, crs = lidR::st_crs(slice_extra)),
              method = "ransac",
              n = n_ransac_samples,
              inliers = ransac_inliers * 0.9,  # Slightly stricter
              conf = ransac_confidence,
              max_angle = max_angle * 0.8,     # Stricter angle constraint
              n_best = n_best
            )
          }, error = function(e) NULL)

          if (!is.null(enhanced_fit) && !is.na(enhanced_fit$err)) {
            # Update result with enhanced fit
            tree_results[[tree_idx]] <- data.frame(
              TreeID = tree_idx,
              X = enhanced_fit$px,
              Y = enhanced_fit$py,
              Z = enhanced_fit$pz,
              Radius = enhanced_fit$radius,
              Error = enhanced_fit$err,
              PointCount = nrow(filtered_points),
              Quality = if (enhanced_fit$err <= ransac_error_threshold) "Enhanced" else "Poor",
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    if (verbose) message("Phase 3 complete: Enhanced processing applied to poor fits")
  }

  # Combine results and create sf object
  if (length(tree_results) == 0) {
    warning("No valid cylinder fits obtained")
    return(create_empty_sf_result())
  }

  final_results <- dplyr::bind_rows(tree_results)
  final_results <- final_results[!is.na(final_results$X), ]

  if (nrow(final_results) == 0) {
    warning("No valid trees after cylinder fitting")
    return(create_empty_sf_result())
  }

  # Create sf object
  trees_sf <- sf::st_as_sf(final_results,
                          coords = c("X", "Y", "Z"),
                          crs = lidR::st_crs(las))

  # Add coordinate columns back for easier access
  coords <- sf::st_coordinates(trees_sf)
  trees_sf$X <- coords[, 1]
  trees_sf$Y <- coords[, 2]
  trees_sf$Z <- coords[, 3]

  if (verbose) {
    quality_summary <- table(trees_sf$Quality)
    message(sprintf("Final results: %d trees detected", nrow(trees_sf)))
    for (quality in names(quality_summary)) {
      message(sprintf("  %s fits: %d", quality, quality_summary[quality]))
    }
  }

  return(trees_sf)
}

#' Create empty sf result for consistent return type
#' @keywords internal
create_empty_sf_result <- function() {
  empty_data <- data.frame(
    TreeID = integer(0),
    Radius = numeric(0),
    Error = numeric(0),
    PointCount = integer(0),
    Quality = character(0),
    X = numeric(0),
    Y = numeric(0),
    Z = numeric(0),
    stringsAsFactors = FALSE
  )

  sf::st_as_sf(empty_data, coords = c("X", "Y", "Z"))
}

#' Validate Detected Trees
#'
#' @description Validates the quality of detected trees based on geometric and statistical criteria
#' @param trees_sf sf object from comprehensive_tree_detection
#' @param min_radius Minimum acceptable radius
#' @param max_radius Maximum acceptable radius
#' @param max_error Maximum acceptable RANSAC error
#' @param min_points Minimum points per tree
#' @return Filtered sf object with only valid trees
#' @export
validate_detected_trees <- function(trees_sf,
                                   min_radius = 0.05,
                                   max_radius = 1.0,
                                   max_error = 0.2,
                                   min_points = 10) {

  if (nrow(trees_sf) == 0) return(trees_sf)

  valid_trees <- trees_sf[
    trees_sf$Radius >= min_radius &
    trees_sf$Radius <= max_radius &
    trees_sf$Error <= max_error &
    trees_sf$PointCount >= min_points,
  ]

  return(valid_trees)
}

#' Plot Tree Detection Results
#'
#' @description Creates a visualization of detected trees overlaid on the point cloud
#' @param las Original LAS object
#' @param trees_sf sf object from comprehensive_tree_detection
#' @param color_by Color trees by this attribute
#' @param point_size Size of tree markers
#' @return ggplot object
#' @export
plot_tree_detection <- function(las, trees_sf, color_by = "Quality", point_size = 3) {
  require(ggplot2)

  if (nrow(trees_sf) == 0) {
    warning("No trees to plot")
    return(NULL)
  }

  # Create canopy height model for background
  chm <- lidR::rasterize_canopy(las, res = 0.2, lidR::p2r())
  chm_df <- as.data.frame(chm, xy = TRUE, na.rm = TRUE)

  # Get tree coordinates
  tree_coords <- sf::st_coordinates(trees_sf)
  plot_data <- data.frame(
    X = tree_coords[, 1],
    Y = tree_coords[, 2],
    Color = trees_sf[[color_by]],
    Radius = trees_sf$Radius,
    Error = trees_sf$Error
  )

  p <- ggplot() +
    geom_raster(data = chm_df, aes(x = x, y = y, fill = Z), alpha = 0.7) +
    scale_fill_viridis_c(name = "Height\n(m)", na.value = "transparent") +
    geom_point(data = plot_data,
               aes(x = X, y = Y, color = Color, size = Radius),
               alpha = 0.8) +
    scale_size_continuous(name = "Radius\n(m)", range = c(1, point_size)) +
    labs(title = "Tree Detection Results",
         subtitle = sprintf("%d trees detected", nrow(trees_sf)),
         x = "X (m)", y = "Y (m)") +
    theme_minimal() +
    coord_fixed()

  if (color_by == "Quality") {
    p <- p + scale_color_manual(values = c("Good" = "green", "Enhanced" = "orange", "Poor" = "red"),
                               name = "Fit Quality")
  }

  return(p)
}
