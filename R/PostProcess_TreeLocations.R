# Function to estimate the crown base height of a tree
estimate_crown_base_height <- function(z, threshold = 0.05) {
  quantiles <- quantile(z, probs = seq(0, 1, 0.05))
  quantile_diffs <- diff(quantiles)
  asymptote_index <- which(quantile_diffs > threshold)[1]
   if (!is.na(asymptote_index)) {
    crown_base_height <- quantiles[asymptote_index + 1]
  } else {
    crown_base_height <- max(z)
  }
  return(crown_base_height)
}

# Function to calculate the convex hull and its volume and surface area
calculate_hull <- function(xyz) {
  if (is.list(xyz) && !is.data.frame(xyz)) {
    p <- as.matrix(do.call("rbind", xyz))
  } else {
    p <- as.matrix(xyz)
  }
  # Wrap the convhulln call in tryCatch to handle errors
  ch <- tryCatch({
    geometry::convhulln(p, "FA")
  }, error = function(e) {
    # Return NA values if an error occurs
    return(list(vol = NA_real_, area = NA_real_))
  })
  return(list(crownvolume = ch$vol, crownsurface = ch$area))
}

# Function to fit a convex hull to points above the crown base height and calculate its volume
fit_convex_hull_and_volume <- function(x, y, z) {
  crown_base_height <- estimate_crown_base_height(z)
  points_above_crown_base <- data.frame(X = x, Y = y, Z = z)[z > crown_base_height, ]
  # Check if there are enough points to calculate a convex hull
  if (nrow(points_above_crown_base) < 4) {
    return(NA_real_)  # Return NA if there are not enough points
  }
  # Calculate the convex hull
  hull <- calculate_hull(points_above_crown_base[, c("X", "Y", "Z")])
  # Check if the hull calculation returned NA values
  if (is.na(hull$crownvolume)) {
    return(NA_real_)  # Return NA if the convex hull calculation failed
  }
  # Return the calculated volume
  return(hull$crownvolume)
}

#' Obtain tree information by processing point cloud data
#'
#' `process_tree_data` processes the output of `get_raster_eigen_treelocs` and `segment_graph` to add information
#' about the height, crown area, and diameter for each unique TreeID. It also has an optional parameter to return
#' an `sf` object representing the convex hulls for each tree.
#'
#' For terrestrial and mobile lidar datasets, tree locations and estimates of DBH are provided
#' by rasterizing individual point cloud values of relative neighborhood density (at 0.3 and
#' 1 m radius) and verticality within a slice of the normalized point cloud around breast height
#' (1.34 m). The algorithm then uses defined threshold values to classify the resulting rasters
#' and create unique polygons from the resulting classified raster. These point-density and
#' verticality polygons were selected by their intersection with one another, resulting in a
#' final set of polygons which were used to clip out regions of the point cloud that were most
#' likely to represent tree boles. A RANSAC cylinder fitting algorithm was then used to estimate
#' the fit of a cylinder to individual bole points. Cylinder centers and radius were used as inputs
#' to an individual tree segmentation.
#'
#' @param treeData An `sf` object containing the following tree information: `TreeID`,
#' `X`, `Y`, `Z`, `Radius`, and `Error`, output from the `get_raster_eigen_treelocs` function.
#' @param segmentedLAS A LAS object that is the output from `segment_graph`.
#' @param return_sf logical: If TRUE, returns an `sf` object representing the convex hulls for each tree.
#' @param seg_table `data.frame` or `NULL`. Optional output of [segment_bole()].
#'   When supplied, per-tree bole volume is merged into the result via
#'   [compute_bole_volume()] (unless `vol_table` is also provided). Required
#'   columns: `TreeID`, `Segment`, `Z_low`, `Z_high`, `Radius`, `Valid`.
#' @param vol_table `data.frame` or `NULL`. Optional output of
#'   [compute_bole_volume()]. When supplied alongside `seg_table`, volume
#'   columns are merged directly without recomputing. Required columns:
#'   `TreeID`, `Volume_m3`.
#' @param qual_table `data.frame` / `sf` or `NULL`. Optional output of
#'   [assess_tree_quality()]. Quality columns (`quality_class`, `fit_score`,
#'   etc.) are merged into the result. Required columns: `TreeID`,
#'   `quality_class`, `fit_score`.
#' @return sf object An updated `sf` object with the original columns plus:
#'   \describe{
#'     \item{height}{numeric: Height of the highest point for each TreeID.}
#'     \item{crown_area}{numeric: Area of the convex hull for each TreeID.}
#'     \item{crown_base_height}{numeric: Height to the base of the live crown for each TreeID.}
#'     \item{crown_volume}{numeric: Volume of the convex hull for the crown of each TreeID.}
#'     \item{diameter}{numeric: Diameter of the tree, calculated as twice the Radius.}
#'   }
#'   If `return_sf` is TRUE, returns an `sf` object where the geometry is the convex hulls for each tree.
#'   If `return_sf` is FALSE, returns an `sf` object with point geometries using treeData.
#'
#' @examples
#'
#' \donttest{
#' # Set the number of threads to use in lidR
#' set_lidr_threads(8)
#'
#' LASfile = system.file("extdata", "TLS_Clip.laz", package="spanner")
#' las = readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' # Don't forget to make sure the las object has a projection
#' sf::st_crs(las) <- 26912
#'
#' # Pre-process the example lidar dataset by classifying the ground  and noise points
#' # using lidR::csf(), normalizing it, and removing outlier points
#' # using lidR::ivf()
#' # las = classify_ground(las, csf(sloop_smooth = FALSE,
#' #                                 class_threshold = 0.5,
#' #                                cloth_resolution = 0.5, rigidness = 1L,
#' #                                 iterations = 500L, time_step = 0.65))
#' # las = normalize_height(las, tin())
#' # las = classify_noise(las, ivf(0.25, 3))
#' # las = filter_poi(las, Classification != LASNOISE)
#'
#' # Plot the non-ground points, colored by height
#' # plot(filter_poi(las, Classification != 2), color = "Z")
#'
#' # Find individual tree locations and attribute data
#' # find tree locations and attribute data
#' myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.025, pt_spacing = 0.0254,
#'                                        dens_threshold = 0.25,
#'                                        neigh_sizes = c(0.25, 0.15, 0.66),
#'                                        eigen_threshold = 0.75,
#'                                        grid_slice_min = 1,
#'                                        grid_slice_max = 2,
#'                                        minimum_polygon_area = 0.005,
#'                                        cylinder_fit_type = "ransac",
#'                                        max_dia = 1,
#'                                        SDvert = 0.33,
#'                                        n_pts = 20,
#'                                        n_best = 25,
#'                                        inliers = 0.9,
#'                                        conf = 0.99,
#'                                        max_angle = 20)
#'
#' # Plot results if trees were found
#' if (!is.null(myTreeLocs) && nrow(myTreeLocs) > 0) {
#'   plot(lidR::rasterize_canopy(las, res = 0.2, p2r()))
#'   symbols(sf::st_coordinates(myTreeLocs)[,1], sf::st_coordinates(myTreeLocs)[,2],
#'           circles = myTreeLocs$Radius^2*3.14, inches = FALSE, add = TRUE, bg = 'black')
#' } else {
#'   message("No tree locations were found. Try adjusting the parameters.")
#' }
#'
#' # Segment the point cloud
#' # For areas with interlocking crowns and trees of different sizes,
#' # enable metabolic scaling to prevent height overestimation
#' myTreeGraph = segment_graph(las = las, tree.locations = myTreeLocs, k = 50,
#'                              distance.threshold = 0.5,
#'                              use.metabolic.scale = FALSE,
#'                              ptcloud_slice_min = 1,
#'                              ptcloud_slice_max = 2,
#'                              subsample.graph = 0.1,
#'                              return.dense = TRUE)
#'
#' # Plot it in 3D colored by treeID
#' plot(myTreeGraph, color = "treeID", pal=spanner_pal())
#'
#' # Process the data
#' processed_data <- process_tree_data(myTreeLocs, myTreeGraph, return_sf = TRUE)
#'
#' # Print the processed data
#' print(processed_data$data)
#' # Print the sf object if return_sf is TRUE
#' if (!is.null(processed_data$sf)) {
#'   print(processed_data$sf)
#' }
#' }
#'
#' @export
process_tree_data <- function(treeData,
                              segmentedLAS,
                              return_sf  = FALSE,
                              seg_table  = NULL,
                              vol_table  = NULL,
                              qual_table = NULL) {

  # ---- Input validation -----------------------------------------------------
  if (!inherits(segmentedLAS, "LAS"))
    stop("'segmentedLAS' must be a LAS object — run segment_graph() first.")
  if (!"treeID" %in% names(segmentedLAS@data))
    stop("'segmentedLAS' must contain a 'treeID' column — run segment_graph() first.")
  if (!inherits(treeData, "sf") && !is.data.frame(treeData))
    stop("'treeData' must be an sf object from get_raster_eigen_treelocs().")
  if (!"TreeID" %in% names(treeData))
    stop("'treeData' must have a 'TreeID' column.")

  if (!is.null(seg_table)) {
    if (!is.data.frame(seg_table))
      stop("'seg_table' must be a data.frame from segment_bole().")
    req_seg <- c("TreeID", "Segment", "Z_low", "Z_high", "Radius", "Valid")
    miss    <- setdiff(req_seg, names(seg_table))
    if (length(miss) > 0)
      stop("'seg_table' is missing columns: ", paste(miss, collapse = ", "),
           "\nRun segment_bole() before process_tree_data().")
  }

  if (!is.null(vol_table)) {
    if (!is.data.frame(vol_table))
      stop("'vol_table' must be a data.frame from compute_bole_volume().")
    req_vol <- c("TreeID", "Volume_m3")
    miss    <- setdiff(req_vol, names(vol_table))
    if (length(miss) > 0)
      stop("'vol_table' is missing columns: ", paste(miss, collapse = ", "),
           "\nRun compute_bole_volume() before process_tree_data().")
  }

  if (!is.null(qual_table)) {
    qual_df_check <- if (inherits(qual_table, "sf"))
      sf::st_drop_geometry(qual_table) else qual_table
    if (!is.data.frame(qual_df_check))
      stop("'qual_table' must be an sf or data.frame from assess_tree_quality().")
    req_qual <- c("TreeID", "quality_class", "fit_score")
    miss     <- setdiff(req_qual, names(qual_df_check))
    if (length(miss) > 0)
      stop("'qual_table' is missing columns: ", paste(miss, collapse = ", "),
           "\nRun assess_tree_quality() before process_tree_data().")
  }

  # ---- TreeID diagnostics --------------------------------------------------
  segmented_tree_ids <- unique(segmentedLAS@data$treeID[
    !is.na(segmentedLAS@data$treeID)])
  tree_data_ids <- unique(treeData$TreeID)

  missing_in_segmented <- setdiff(tree_data_ids, segmented_tree_ids)
  missing_in_treedata  <- setdiff(segmented_tree_ids, tree_data_ids)

  if (length(missing_in_segmented) > 0)
    message("TreeIDs in treeData but not in segmentedLAS: ",
            paste(missing_in_segmented, collapse = ", "))
  if (length(missing_in_treedata) > 0)
    message("TreeIDs in segmentedLAS but not in treeData: ",
            paste(missing_in_treedata, collapse = ", "))
  if (!all(treeData$TreeID %in% segmented_tree_ids) ||
      length(tree_data_ids) != length(segmented_tree_ids))
    warning("TreeIDs do not match between treeData and segmentedLAS.")

  unique_tree_ids <- unique(treeData$TreeID)
  n_trees <- length(unique_tree_ids)

  # ---- Pre-split LAS data.table once (avoids repeated filter_poi scans) ----
  las_dt    <- segmentedLAS@data
  tid_split <- split(seq_len(nrow(las_dt)),
                     las_dt$treeID,
                     drop = FALSE)

  # ---- Per-tree geometry metrics -------------------------------------------
  heights            <- numeric(n_trees)
  crown_areas        <- numeric(n_trees)
  crown_base_heights <- numeric(n_trees)
  crown_volumes      <- numeric(n_trees)
  convex_hulls       <- vector("list", n_trees)

  for (i in seq_along(unique_tree_ids)) {
    tree_id <- unique_tree_ids[i]
    row_idx <- tid_split[[as.character(tree_id)]]

    if (is.null(row_idx) || length(row_idx) == 0L) {
      heights[i] <- crown_areas[i] <- crown_base_heights[i] <-
        crown_volumes[i] <- NA_real_
      next
    }

    tX <- las_dt$X[row_idx]
    tY <- las_dt$Y[row_idx]
    tZ <- las_dt$Z[row_idx]

    heights[i] <- max(tZ, na.rm = TRUE)

    if (length(tX) > 2L) {
      convex_hull       <- sf::st_convex_hull(
        sf::st_union(sf::st_sfc(sf::st_multipoint(cbind(tX, tY)))))
      convex_hulls[[i]] <- convex_hull
      crown_areas[i]    <- as.numeric(sf::st_area(convex_hull))
    } else {
      crown_areas[i] <- NA_real_
    }

    crown_base_heights[i] <- estimate_crown_base_height(tZ)
    crown_volumes[i]      <- fit_convex_hull_and_volume(tX, tY, tZ)
  }

  treeData$height            <- heights
  treeData$crown_area        <- crown_areas
  treeData$diameter          <- treeData$Radius * 2
  treeData$crown_base_height <- crown_base_heights
  treeData$crown_volume      <- crown_volumes

  # ---- Merge optional pipeline-2 tables ------------------------------------
  if (!is.null(vol_table)) {
    keep_vol <- intersect(
      names(vol_table),
      c("TreeID", "Volume_m3", "N_segments", "Mean_radius",
        "Height_modeled", "Basal_area_m2", "Mean_RANSAC_err"))
    vol_merge <- vol_table[, keep_vol, drop = FALSE]
    treeData  <- merge(treeData, vol_merge,
                       by = "TreeID", all.x = TRUE, sort = FALSE)
  }

  if (!is.null(qual_table)) {
    qual_plain <- if (inherits(qual_table, "sf"))
      sf::st_drop_geometry(qual_table) else qual_table
    # Exclude geometry and any columns already present in treeData
    keep_qual <- setdiff(names(qual_plain),
                         c("geometry", setdiff(names(treeData), "TreeID")))
    keep_qual <- union("TreeID", keep_qual)
    treeData  <- merge(treeData, qual_plain[, keep_qual, drop = FALSE],
                       by = "TreeID", all.x = TRUE, sort = FALSE)
  }

  # ---- Build return object --------------------------------------------------
  if (return_sf) {
    named_hulls <- stats::setNames(convex_hulls, as.character(unique_tree_ids))
    hulls_sf <- do.call(rbind, lapply(seq_along(unique_tree_ids), function(i) {
      tid <- unique_tree_ids[i]
      if (!is.null(named_hulls[[as.character(tid)]]))
        sf::st_sf(TreeID = tid,
                  geometry = named_hulls[[as.character(tid)]])
    }))
    result          <- sf::st_sf(treeData)
    result$geometry <- hulls_sf$geometry
  } else {
    result <- sf::st_as_sf(treeData, coords = c("X", "Y"),
                           crs = sf::st_crs(segmentedLAS))
  }

  return(result)
}
