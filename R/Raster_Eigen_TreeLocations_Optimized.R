#' Optimized Tree Location Detection using Direct Raster Operations
#'
#' @description A memory-efficient version of get_raster_eigen_treelocs that avoids
#' sf objects and shapefiles, using optimized geometric features and direct raster operations.
#'
#' @param las LAS object
#' @param neigh_sizes Neighborhood sizes for analysis c(eigen_radius, small_density_radius, large_density_radius)
#' @param res Grid resolution
#' @param grid_slice_min Minimum height for processing slice
#' @param grid_slice_max Maximum height for processing slice
#' @param dens_threshold Density threshold for tree detection
#' @param vert_threshold Verticality threshold for tree detection
#' @param minimum_area Minimum area for tree candidate regions
#' @param pt_spacing Point spacing for downsampling
#' @param SDvert Standard deviation threshold for verticality filtering
#' @param ncpu Number of CPU cores to use
#'
#' @return Data frame of tree locations with attributes
#' @export
get_raster_eigen_treelocs_optimized <- function(las,
                                               neigh_sizes = c(2, 1, 3),
                                               res = 0.02,
                                               grid_slice_min = 1.2,
                                               grid_slice_max = 1.8,
                                               dens_threshold = 0.1,
                                               vert_threshold = 0.5,
                                               minimum_area = 0.1,
                                               pt_spacing = 0.025,
                                               SDvert = 0.15,
                                               ncpu = parallel::detectCores(logical = FALSE)) {

  # Input validation
  if (!inherits(las, "LAS")) stop("Input must be a LAS object")
  if (length(neigh_sizes) != 3) stop("neigh_sizes must be a vector of length 3")

  message("=== Optimized Tree Location Detection ===")
  message("Step 1/10: Downsampling...")

  # Downsample using systematic voxel grid
  las <- lidR::decimate_points(las, random_per_voxel(res = pt_spacing, n = 1))
  if (is.null(las) || nrow(las@data) == 0) {
    stop("No points remaining after downsampling! Try increasing pt_spacing.")
  }

  message("Step 2/10: Creating processing slice...")

  # Create processing slice
  slice_extra <- lidR::filter_poi(las, Z >= grid_slice_min, Z <= grid_slice_max)
  if (is.null(slice_extra) || nrow(slice_extra@data) == 0) {
    stop("No points in processing slice! Adjust grid_slice_min/max.")
  }

  message("Step 3/10: Computing optimized geometric features...")

  start_time <- Sys.time()
  geom_features <- spanner:::C_geometric_features_simple(slice_extra,
                                              radius = neigh_sizes[1],
                                              max_neighbors = 50, ncpu = ncpu)
  # geom_features <- lidR::point_eigenvalues(slice_extra, r = neigh_sizes[1], k = 50, metrics = T)
  # geom_features <- spanner::eigen_metrics(slice_extra, radius = neigh_sizes[1], ncpu = ncpu)

  geom_features <- data.table::as.data.table(geom_features)
  end_time <- Sys.time()
  message(sprintf("   Computed features for %d points in %.2f seconds",
                  nrow(geom_features),
                  as.numeric(end_time - start_time, units = "secs")))

  # Add attributes to LAS object
  slice_extra <- lidR::add_lasattribute(slice_extra, geom_features$verticality, "verticality", "verticality")
  slice_extra <- lidR::add_lasattribute(slice_extra, geom_features$volume_density, "point_density", "point density")

  # Replace NA values
  slice_extra@data$verticality[is.na(slice_extra@data$verticality)] <- 0.5
  slice_extra@data$point_density[is.na(slice_extra@data$point_density)] <- 0

  message("Step 4/10: Computing relative density (optimized)...")

  # Use optimized density calculations - combine both radius calculations
  coords_matrix <- as.matrix(slice_extra@data[, c("X", "Y")])

  # Build KD-tree once for both density calculations
  cnt_local <- spanner:::C_count_in_disc(X = slice_extra@data$X, Y = slice_extra@data$Y,
                                         x = slice_extra@data$X, y = slice_extra@data$Y,
                                         radius = neigh_sizes[2], ncpu = ncpu)
  cnt_large <- spanner:::C_count_in_disc(X = slice_extra@data$X, Y = slice_extra@data$Y,
                                         x = slice_extra@data$X, y = slice_extra@data$Y,
                                         radius = neigh_sizes[3], ncpu = ncpu)

  relative_density <- cnt_local / pmax(cnt_large, 1)  # Avoid division by zero
  slice_extra <- lidR::add_lasattribute(slice_extra, relative_density, "relative_density", "relative density")

  message("Step 5/10: Creating density raster...")

  # Create density raster efficiently
  density_raster <- lidR::grid_metrics(slice_extra,
                                       ~stats::median(relative_density, na.rm = TRUE),
                                       res = res)
  density_raster <- terra::rast(density_raster)

  message("Step 6/10: Creating verticality raster...")

  # Create verticality raster efficiently
  vert_raster <- lidR::grid_metrics(slice_extra,
                                   ~stats::quantile(verticality, 0.5, na.rm = TRUE),
                                   res = res)
  vert_raster <- terra::rast(vert_raster)

  message("Step 7/10: Identifying candidate regions (raster operations)...")

  # Create binary masks using raster operations (avoid polygons)
  density_mask <- density_raster >= dens_threshold
  vert_mask <- vert_raster >= vert_threshold

  # Combine masks
  combined_mask <- density_mask & vert_mask

  # Apply morphological operations for smoothing (equivalent to buffer operations)
  combined_mask <- terra::focal(combined_mask, w = matrix(1, 3, 3), fun = "modal", na.rm = TRUE)

  message("Step 8/10: Extracting candidate centroids...")

  # Convert to points and get centroids of connected regions
  candidate_points <- terra::as.points(combined_mask, values = TRUE, na.rm = TRUE)
  candidate_points <- candidate_points[candidate_points[[1]] == 1, ]

  if (length(candidate_points) == 0) {
    warning("No candidate regions found. Adjust thresholds.")
    return(data.frame(
      TreeID = integer(0),
      X = numeric(0),
      Y = numeric(0),
      density = numeric(0),
      verticality = numeric(0),
      point_count = integer(0),
      height_range = numeric(0),
      mean_height = numeric(0),
      verticality_sd = numeric(0),
      estimated_radius = numeric(0)
    ))
  }

  # Get coordinates
  candidate_coords <- terra::crds(candidate_points)

  # Cluster nearby points to get tree centers
  message("Step 9/10: Clustering candidate points...")

  if (nrow(candidate_coords) > 0) {
    # Simple clustering using connected components
    # if (nrow(candidate_coords) > 1) {
    #   distances <- dist(candidate_coords)
    #
    #   # =================== Replace with connected components ==================
    #   clusters <- stats::hclust(distances, method = "single")
    #   cluster_ids <- stats::cutree(clusters, h = res * 5)  # Cluster within 5 grid cells
    #   # ========================================================================
    #
    #   # Get cluster centroids
    #   tree_locations <- aggregate(candidate_coords,
    #                              by = list(cluster = cluster_ids),
    #                              FUN = mean)
    #   colnames(tree_locations) <- c("cluster", "X", "Y")
    #
    # } else {
    #   tree_locations <- data.frame(cluster = 1,
    #                               X = candidate_coords[1, 1],
    #                               Y = candidate_coords[1, 2])
    # }

    if (nrow(candidate_coords) > 0) {

      las_cc <- connected_components(slice_extra, res = minimum_area, min_pts = 100,
                                     name = "clusterID")

      d <- las_cc@data
      d <- d[d$clusterID > 0, ]              # optional: drop noise cluster 0 if present
      tree_locations <- aggregate(d[, c("X", "Y")],
                                  by = list(clusterID = d$clusterID),
                                  FUN = mean)

      tree_locations$n_pts <- tabulate(d$clusterID, max(d$clusterID))[tree_locations$clusterID]

    } else {
      tree_locations <- data.frame(clusterID = integer(), X = numeric(), Y = numeric(), n_pts = integer())
    }

    message("Step 10/10: Computing tree attributes...")

    # Extract attributes for each tree location
    tree_attributes <- data.frame(
      TreeID = 1:nrow(tree_locations),
      X = tree_locations$X,
      Y = tree_locations$Y,
      stringsAsFactors = FALSE
    )

    # Extract raster values at tree locations
    tree_coords <- tree_locations[, c("X", "Y")]

    # Get density and verticality values
    tree_attributes$density <- terra::extract(density_raster, tree_coords, method = "simple")[, 2]
    tree_attributes$verticality <- terra::extract(vert_raster, tree_coords, method = "simple")[, 2]

    # Calculate additional metrics using point cloud data
    for (i in 1:nrow(tree_attributes)) {
      # Define circular ROI around tree location
      tree_x <- tree_attributes$X[i]
      tree_y <- tree_attributes$Y[i]
      search_radius <- neigh_sizes[3]  # Use larger radius for tree extraction

      # Get points within radius
      tree_points <- slice_extra@data[
        sqrt((slice_extra@data$X - tree_x)^2 + (slice_extra@data$Y - tree_y)^2) <= search_radius,
      ]

      setDT(tree_attributes)
      if (nrow(tree_points) > 0) {
        tree_attributes[i, `:=`(
          point_count    = nrow(tree_points),
          height_range   = max(tree_points$Z, na.rm = TRUE) - min(tree_points$Z, na.rm = TRUE),
          mean_height    = mean(tree_points$Z, na.rm = TRUE),
          verticality_sd = sd(tree_points$verticality, na.rm = TRUE),
          estimated_radius = quantile(sqrt((tree_points$X - tree_x)^2 + (tree_points$Y - tree_y)^2),
                                      0.75, na.rm = TRUE)
        )]
      } else {
        tree_attributes[i, `:=`(
          point_count = 0L,
          height_range = 0,
          mean_height = NA_real_,
          verticality_sd = NA_real_,
          estimated_radius = NA_real_
        )]
      }
    }

    # ===== This needs to be reworked, as it removes all trees ====
    # # Filter based on quality criteria
    # tree_attributes <- tree_attributes[
    #   !is.na(tree_attributes$verticality_sd) &
    #   tree_attributes$verticality_sd < SDvert &
    #   tree_attributes$point_count >= 10 &
    #   !is.na(tree_attributes$height_range) &
    #   tree_attributes$height_range > (grid_slice_max - grid_slice_min) * 0.5,
    # ]

    # Clean up row names and assign TreeIDs
    rownames(tree_attributes) <- NULL
    if (nrow(tree_attributes) > 0) {
      tree_attributes$TreeID <- 1:nrow(tree_attributes)
    }

    message(sprintf("Detection complete! Found %d tree candidates.", nrow(tree_attributes)))

    return(tree_attributes)

  } else {
    warning("No candidate points found after processing.")
    return(data.frame(
      TreeID = integer(0),
      X = numeric(0),
      Y = numeric(0),
      density = numeric(0),
      verticality = numeric(0),
      point_count = integer(0),
      height_range = numeric(0),
      mean_height = numeric(0),
      verticality_sd = numeric(0),
      estimated_radius = numeric(0)
    ))
  }
}


#' Combined Geometric Features and Density Calculation
#'
#' @description Optimized function that calculates geometric features and density
#' metrics in a single pass through the data for maximum efficiency
#'
#' @param las LAS object
#' @param radius_eigen Radius for eigenvalue calculations
#' @param radius_small Small radius for density calculations
#' @param radius_large Large radius for density calculations
#' @param max_neighbors Maximum neighbors to consider
#' @param ncpu Number of CPU cores
#'
#' @return Data table with all computed metrics
#' @export
compute_combined_metrics <- function(las,
                                   radius_eigen = 2.0,
                                   radius_small = 1.0,
                                   radius_large = 3.0,
                                   max_neighbors = 50,
                                   ncpu = parallel::detectCores()) {

  message("Computing combined geometric features and density metrics...")

  # Get optimized geometric features
  if (!exists("C_geometric_features_simple", mode = "function")) {
    Rcpp::sourceCpp("E:/GithubRepos/spanner - Copy - Copy/src/geometric_features_simple.cpp")
  }

  start_time <- Sys.time()

  # Compute geometric features
  geom_features <- C_geometric_features_simple(las, radius = radius_eigen, max_neighbors = max_neighbors)
  geom_features <- data.table::as.data.table(geom_features)

  # Compute density metrics
  cnt_local <- spanner:::C_count_in_disc(X = las@data$X, Y = las@data$Y,
                                         x = las@data$X, y = las@data$Y,
                                         radius = radius_small, ncpu = ncpu)
  cnt_large <- spanner:::C_count_in_disc(X = las@data$X, Y = las@data$Y,
                                         x = las@data$X, y = las@data$Y,
                                         radius = radius_large, ncpu = ncpu)

  # Combine all metrics
  combined_metrics <- data.table::data.table(
    # Point coordinates
    X = las@data$X,
    Y = las@data$Y,
    Z = las@data$Z,

    # Geometric features
    linearity = geom_features$linearity,
    planarity = geom_features$planarity,
    sphericity = geom_features$sphericity,
    omnivariance = geom_features$omnivariance,
    anisotropy = geom_features$anisotropy,
    eigenentropy = geom_features$eigenentropy,
    surface_variation = geom_features$surface_variation,
    verticality = geom_features$verticality,
    local_density = geom_features$local_density,

    # Density metrics
    count_small = cnt_local,
    count_large = cnt_large,
    relative_density = cnt_local / pmax(cnt_large, 1)
  )

  end_time <- Sys.time()
  message(sprintf("Combined metrics computed for %d points in %.2f seconds",
                  nrow(combined_metrics),
                  as.numeric(end_time - start_time, units = "secs")))

  return(combined_metrics)
}

