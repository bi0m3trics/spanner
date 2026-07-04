# ============================================================================
# spanner v2.0 – TLS_Pipeline.R
#
# Explicit staged TLS/MLS pipeline.
#
# Stage 1  find_trees()      Detect tree locations → sf
# Stage 2  tree_points()    Assign points to trees → LAS + treeID
# Stage 3  stem_points()    Classify stem points → LAS + Stem cols
# Stage 4  stem_segments()  Fit per-segment geometry → data.frame
# Stage 5  tree_inventory()  Build tree-level inventory → sf
#
# Backwards-compatible wrapper: get_raster_eigen_treelocs() + segment_graph()
# remain fully functional; the new stages call them internally or alongside
# the new Hough-based algorithms.
# ============================================================================


# ----------------------------------------------------------------------------
#' Detect tree locations in a point cloud (stage 1)
#'
#' Applies the chosen detection method to a normalized LAS object and returns
#' an `sf` POINT object — one point per detected tree — that is compatible
#' with all downstream pipeline stages.
#'
#' @details
#' Two method factories are available:
#' \describe{
#'   \item{[use_raster_eigen()]}{Wraps [get_raster_eigen_treelocs()].
#'     Rasterizes point-cloud density and verticality, intersects the
#'     resulting polygons, and fits RANSAC cylinders to derive DBH estimates.
#'     Highest precision for clean TLS data; computationally heavier.}
#'   \item{[use_hough()]}{Runs the circular Hough Transform over stacked
#'     height slices and selects vertically persistent circle centres as tree
#'     positions.  Faster and more tolerant of occlusion; radius estimates
#'     are less precise than the raster-eigen approach.}
#' }
#'
#' Both methods return an `sf` object with columns `TreeID` (integer),
#' `Radius` (numeric, m), and `Error` (numeric or `NA`).
#'
#' @param las    LAS object.  Must be height-normalized (ground = Z ≈ 0).
#' @param method A method closure created by [use_raster_eigen()] (default)
#'   or [use_hough()].
#'
#' @return An `sf` POINT object with columns `TreeID`, `Radius`, `Error` in
#'   the same coordinate reference system as `las`.
#'
#' @examples
#' \donttest{
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#'
#' # --- Method A: raster-eigen (high-quality TLS) ---
#' map_re <- find_trees(las, method = use_raster_eigen(
#'   res = 0.025, pt_spacing = 0.0254, dens_threshold = 0.25,
#'   neigh_sizes = c(0.25, 0.15, 0.66), eigen_threshold = 0.75,
#'   grid_slice_min = 1, grid_slice_max = 2,
#'   minimum_polygon_area = 0.005, cylinder_fit_type = "ransac",
#'   max_dia = 1, SDvert = 0.33, n_pts = 20, n_best = 25,
#'   inliers = 0.9, conf = 0.99, max_angle = 20
#' ))
#'
#' # --- Method B: Hough slice-stack (fast, handles occlusion) ---
#' map_hg <- find_trees(las, method = use_hough(
#'   min_h = 1, max_h = 3, h_step = 0.5,
#'   max_d = 0.5, min_density = 0.1, min_votes = 3
#' ))
#' }
#'
#' @seealso [use_raster_eigen()], [use_hough()], [tree_points()],
#'   [get_raster_eigen_treelocs()]
#' @export
find_trees <- function(las, method = use_raster_eigen()) {
  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object", call. = FALSE)
  if (!inherits(method, "spanner_map_mtd"))
    stop("'method' must be created by use_raster_eigen() or use_hough()",
         call. = FALSE)
  method(las)
}


# ----------------------------------------------------------------------------
#' Assign point-cloud points to individual trees (stage 2)
#'
#' Given tree seed locations from [find_trees()], assigns every point a
#' `treeID` integer column using the chosen spatial assignment method.
#'
#' @details
#' Three assignment methods are available:
#' \describe{
#'   \item{[assign_voronoi()]}{Nearest-tree Voronoi partition — every point
#'     is assigned, regardless of distance.  Fast (O(N × T) brute force for
#'     typical T < 500).  Best when the entire plot should be fully
#'     partitioned.}
#'   \item{[assign_crop()]}{Fixed-radius (or square) crop around each tree
#'     centre.  Conservative — only points close to the stem receive a
#'     `treeID`; crown/canopy points may stay at 0.  Good pre-filter before
#'     stem fitting.}
#'   \item{[assign_graph()]}{Wraps [segment_graph()] — builds a k-NN point
#'     graph and propagates IDs via shortest paths.  Most spatially accurate
#'     for crown delineation; also the most computationally expensive.}
#' }
#'
#' @param las    LAS object.  Height-normalized.
#' @param map    `sf` POINT object from [find_trees()] (or directly from
#'   [get_raster_eigen_treelocs()]).
#' @param method A method closure created by [assign_voronoi()] (default),
#'   [assign_crop()], or [assign_graph()].
#'
#' @return A LAS object identical to `las` with a new integer column
#'   `treeID`.  Points not assigned to any tree receive `treeID = 0`.
#'
#' @examples
#' \donttest{
#' map <- find_trees(las, use_raster_eigen())
#'
#' # Voronoi – every point gets a tree ID
#' las_v <- tree_points(las, map, method = assign_voronoi())
#'
#' # Crop – only points within 1.5 m of each stem
#' las_c <- tree_points(las, map, method = assign_crop(l = 1.5))
#'
#' # Graph – full crown-level graph segmentation
#' las_g <- tree_points(las, map,
#'            method = assign_graph(k = 50, distance.threshold = 0.5,
#'                                  subsample.graph = 0.1,
#'                                  return.dense = TRUE))
#' }
#'
#' @seealso [find_trees()], [assign_voronoi()], [assign_crop()],
#'   [assign_graph()], [segment_graph()]
#' @export
tree_points <- function(las, map, method = assign_graph()) {
  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object", call. = FALSE)
  if (!inherits(map, "sf"))
    stop("'map' must be an sf object (output of find_trees())", call. = FALSE)
  if (!inherits(method, "spanner_tpt_mtd"))
    stop("'method' must be created by assign_voronoi(), assign_crop(), ",
         "or assign_graph()", call. = FALSE)
  method(las, map)
}


# ----------------------------------------------------------------------------
#' Classify stem points within a segmented point cloud (stage 3)
#'
#' Labels each point as stem or non-stem using the chosen method.
#' This stage is optional; [stem_segments()] can operate on any LAS that has
#' a `treeID` column, but stem-point pre-filtering greatly improves fit quality.
#'
#' @details
#' Two methods are available:
#' \describe{
#'   \item{[stem_eigen()]}{Wraps [classify_lw_points()].  Uses local PCA
#'     eigenvalue metrics (Linearity, Verticality) with optional LeWoS
#'     label propagation.  Returns `Stem`, `StemProb`, `Branch`,
#'     `BranchProb`, `Other` columns.  Best on dense TLS; pass
#'     `neigh_k = 30L` and lower `verticality_threshold` for MLS.}
#'   \item{[stem_hough()]}{Uses `C_hough_stem_plot()` to track circular
#'     cross-sections upward per tree and label points inside the tracked
#'     cylinder.  Returns `Stem`, `Segment`, `Radius`, `Votes` columns.
#'     More robust to mixed-return and MLS noise; does not separate branches.}
#' }
#'
#' @param las    LAS object with a `treeID` column (output of [tree_points()]).
#' @param map    `sf` POINT object from [find_trees()].
#' @param method A method closure created by [stem_eigen()] (default) or
#'   [stem_hough()].
#'
#' @return A LAS object identical to `las` with classification columns added
#'   (see the individual method documentation for column names).
#'
#' @examples
#' \donttest{
#' # Eigen + LeWoS label propagation (TLS-optimised; add neigh_k=30 for MLS)
#' las_e <- stem_points(las_v, map, method = stem_eigen(neigh_k = 20L))
#'
#' # Hough slice tracker (preferred for MLS)
#' las_h <- stem_points(las_v, map,
#'            method = stem_hough(h_base_lo = 1, h_base_hi = 2.5,
#'                                max_d = 0.5, min_votes = 5L))
#' }
#'
#' @seealso [tree_points()], [stem_eigen()], [stem_hough()],
#'   [classify_lw_points()]
#' @export
stem_points <- function(las, map, method = stem_eigen()) {
  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object", call. = FALSE)
  if (!"treeID" %in% names(las@data))
    stop("'las' must have a 'treeID' column — run tree_points() first.",
         call. = FALSE)
  if (!inherits(map, "sf"))
    stop("'map' must be an sf object (output of find_trees())", call. = FALSE)
  if (!inherits(method, "spanner_stm_mtd"))
    stop("'method' must be created by stem_eigen() or stem_hough()", call. = FALSE)
  method(las, map)
}


# ----------------------------------------------------------------------------
#' Fit per-segment stem geometry (stage 4)
#'
#' Divides each tree stem into height slices and fits a geometric primitive
#' (circle or cylinder) to each slice, returning a tidy `data.frame` with
#' one row per (TreeID × segment).
#'
#' @details
#' Two methods are available (see individual documentation for parameters):
#' \describe{
#'   \item{[seg_ransac_cylinder()]}{3-D RANSAC cylinder per slice — captures
#'     stem lean and axis orientation.  Default.}
#'   \item{[seg_irls_cylinder()]}{3-D IRLS cylinder — fast for clean data;
#'     less robust to outliers.}
#' }
#'
#' @param las    LAS object with `treeID` column (and optionally `Stem` /
#'   `Stem` column from [stem_points()]).
#' @param map    `sf` POINT object from [find_trees()].
#' @param method A method closure created by [seg_ransac_cylinder()] (default)
#'   or [seg_irls_cylinder()].
#'
#' @return A `data.frame` with one row per (TreeID, Segment) containing at
#'   minimum: `TreeID` (int), `Segment` (int), `Z_low`, `Z_high`, `Z_mid`
#'   (num), `Radius` (num, m), `Diameter` (num, m), `RANSAC_err` (num,
#'   RMS fit residual when available), `N_pts` (int), `Valid` (logical),
#'   `Fit_method` (chr).
#'
#' @examples
#' \donttest{
#' segs_r <- stem_segments(las_e, map,
#'             method = seg_ransac_cylinder(dz = 0.5, max_angle = 20))
#'
#' segs_i <- stem_segments(las_e, map, method = seg_irls_cylinder())
#' }
#'
#' @seealso [stem_points()], [seg_ransac_cylinder()], [seg_irls_cylinder()]
#' @export
stem_segments <- function(las, map, method = seg_ransac_cylinder()) {
  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object", call. = FALSE)
  if (!"treeID" %in% names(las@data))
    stop("'las' must have a 'treeID' column — run tree_points() first.",
         call. = FALSE)
  if (!inherits(map, "sf"))
    stop("'map' must be an sf object (output of find_trees())", call. = FALSE)
  if (!inherits(method, "spanner_sgs_mtd"))
    stop("'method' must be created by seg_ransac_cylinder() or seg_irls_cylinder()",
         call. = FALSE)
  method(las, map)
}


# ----------------------------------------------------------------------------
#' Compute a TLS/MLS tree inventory (stage 5)
#'
#' Final stage of the spanner v2.0 pipeline.  Wraps [process_tree_data()] to
#' compute per-tree height, crown, and diameter metrics from a segmented point
#' cloud and optional segment-fitting table.
#'
#' @details
#' At minimum, `map` and `segmented_las` (with `treeID` column) are required.
#' Supplying `seg_table` (output of [stem_segments()] or [segment_stem()])
#' enables accurate DBH and stem volume computation.
#'
#' The function normalises the `seg_table` structure so that cylinder-method
#' outputs are accepted by [process_tree_data()] without modification.
#'
#' @param map            `sf` POINT object from [find_trees()].
#' @param segmented_las  LAS object from [tree_points()] (or [stem_points()]).
#'   Must have a `treeID` column.
#' @param seg_table      `data.frame` or `NULL`.  Output of [stem_segments()]
#'   or [segment_stem()].  Required columns: `TreeID`, `Segment`, `Z_low`,
#'   `Z_high`, `Radius`, `Valid`.  Default `NULL`.
#' @param vol_table      `data.frame` or `NULL`.  Output of
#'   [compute_stem_volume()].  When supplied, stem volume is merged directly.
#'   Default `NULL`.
#' @param qual_table     `sf` or `NULL`.  Output of [assess_tree_quality()].
#'   Default `NULL`.
#' @param return_sf      logical.  `TRUE` = return an `sf` object with convex-
#'   hull crown polygons as geometry.  `FALSE` (default) = return an `sf`
#'   object with point geometries.
#' @param dbh_z_lo       numeric.  Lower height (m) for DBH band extraction.
#'   Currently informational — used by downstream functions.  Default `1.2`.
#' @param dbh_z_hi       numeric.  Upper height (m) for DBH band.
#'   Default `1.5`.
#'
#' @return An `sf` object with the tree-seed columns from `map` plus:
#'   `height`, `crown_area`, `crown_base_height`, `crown_volume`, `diameter`.
#'   Additional columns from `seg_table`, `vol_table`, and `qual_table` when
#'   supplied.
#'
#' @examples
#' \donttest{
#' # Complete 5-stage pipeline
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#'
#' map   <- find_trees(las, use_raster_eigen())
#' las_t <- tree_points(las,   map, assign_voronoi())
#' las_s <- stem_points(las_t, map, stem_eigen())
#' segs  <- stem_segments(las_s, map, seg_ransac_cylinder())
#' inv   <- tree_inventory(map, las_t, seg_table = segs)
#' print(inv)
#' }
#'
#' @seealso [find_trees()], [tree_points()], [stem_points()],
#'   [stem_segments()], [process_tree_data()], [segment_stem()],
#'   [assess_tree_quality()]
#' @export
tree_inventory <- function(map,
                           segmented_las,
                           seg_table   = NULL,
                           vol_table   = NULL,
                           qual_table  = NULL,
                           return_sf   = FALSE,
                           dbh_z_lo    = 1.2,
                           dbh_z_hi    = 1.5) {

  if (!inherits(map, "sf"))
    stop("'map' must be an sf object (output of find_trees())", call. = FALSE)
  if (!inherits(segmented_las, "LAS"))
    stop("'segmented_las' must be a LAS object", call. = FALSE)
  if (!"treeID" %in% names(segmented_las@data))
    stop("'segmented_las' must have a 'treeID' column — run tree_points() first.",
         call. = FALSE)

  # Normalise seg_table: cylinder methods omit RANSAC_err / N_inliers / Inlier_frac
  # but process_tree_data() needs the segment_stem() column set.
  if (!is.null(seg_table)) {
    if (!is.data.frame(seg_table))
      stop("'seg_table' must be a data.frame (output of stem_segments() or segment_stem())",
           call. = FALSE)
    if (!"RANSAC_err" %in% names(seg_table))  seg_table$RANSAC_err   <- NA_real_
    if (!"N_inliers"  %in% names(seg_table))  seg_table$N_inliers    <- NA_integer_
    if (!"Inlier_frac"%in% names(seg_table))  seg_table$Inlier_frac  <- NA_real_
    if (!"X_ctr"      %in% names(seg_table)) {
      seg_table$X_ctr <- NA_real_
      seg_table$Y_ctr <- NA_real_
    }
  }

  process_tree_data(
    treeData     = map,
    segmentedLAS = segmented_las,
    return_sf    = return_sf,
    seg_table    = seg_table,
    vol_table    = vol_table,
    qual_table   = qual_table
  )
}
