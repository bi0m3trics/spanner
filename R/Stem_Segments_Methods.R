# ============================================================================
# spanner v2.0 – Stem_Segments_Methods.R
#
# Method factories for the stem_segments() stage.
# Each factory returns a closure of class "spanner_sgs_mtd" that accepts
# (las, map) and returns a data.frame of per-tree, per-segment geometry.
#
# Two methods are provided:
#   seg_ransac_cylinder() – 3-D RANSAC cylinder per height slice
#   seg_irls_cylinder()   – 3-D IRLS cylinder (fast, less robust to outliers)
# ============================================================================


# Internal helper: resolve the stem column and XY search radius.
.sgs_resolve <- function(las, map, use_stem_column, search_radius) {
  stem_col <- if ("Stem" %in% names(las@data)) "Stem" else NULL
  has_stem <- use_stem_column && !is.null(stem_col)

  sr <- if (is.null(search_radius)) {
    max(map$Radius, na.rm = TRUE) * 3.0
  } else {
    search_radius
  }
  sr <- max(sr, 0.2)
  list(stem_col = stem_col, has_stem = has_stem, sr = sr)
}


# Internal helper: batch cylinder fitting using plot-level C++ functions.
#
# Replaces the old per-slice R loop.  Segment IDs are pre-computed from Z
# heights (vectorised), then a single C++ call processes all trees x slices.
# 'cpp_fn' must accept (las_matrix, tId, segs, rads, ...) and return the
# nested list produced by ransacPlotCylinders / irlsPlotCylinders.
#
# The cylinder parameter vector layout (0-indexed in C++):
#   [phi, theta, x_ctr, y_ctr, radius, ..., seg_idx, tree_id]
# After do.call(rbind, ...) column 5 (1-indexed) = radius,
# second-to-last = seg_idx (0-based), last = tree_id.
.sgs_plot_cylinders <- function(las, map, dz, overlap, max_radius, min_pts,
                                 use_stem_column, search_radius,
                                 fit_method_name, cpp_fn) {
  ctx    <- .sgs_resolve(las, map, use_stem_column, search_radius)
  tl_ids <- as.integer(map$TreeID)

  # --- select qualifying points ------------------------------------------
  tid_vec <- las@data$treeID
  keep    <- !is.na(tid_vec) & tid_vec %in% tl_ids & tid_vec > 0L
  if (ctx$has_stem)
    keep <- keep & !is.na(las@data[[ctx$stem_col]]) & las@data[[ctx$stem_col]]
  if (sum(keep) < min_pts) return(data.frame())

  Z_sub  <- las@data$Z[keep]
  X_sub  <- las@data$X[keep]
  Y_sub  <- las@data$Y[keep]
  tid_sub <- as.numeric(tid_vec[keep])   # C++ wants numeric (unsigned int)

  # --- per-point segment index (0-based) from Z height -------------------
  z_start <- 0.1            # skip near-ground noise
  step    <- dz * (1.0 - overlap)
  seg_idx <- pmax(0.0, floor((Z_sub - z_start) / step))

  # --- per-point initial radius from tree map ----------------------------
  rad_lookup <- setNames(map$Radius, as.character(tl_ids))
  rad_vec    <- rad_lookup[as.character(as.integer(tid_sub))]
  rad_vec[is.na(rad_vec)] <- 0.15   # fallback for unmatched

  # --- single C++ call for all trees x slices ---------------------------
  las_mat <- cbind(X_sub, Y_sub, Z_sub)
  raw     <- cpp_fn(las_mat, tid_sub, seg_idx, as.numeric(rad_vec))

  if (is.null(raw) || length(raw) == 0) return(data.frame())

  # Flatten nested list: [[tree]][[seg]] -> matrix rows
  flat <- do.call(rbind, lapply(raw, function(t) {
    if (length(t) == 0) return(NULL)
    do.call(rbind, t)
  }))
  if (is.null(flat) || nrow(flat) == 0) return(data.frame())
  flat <- as.data.frame(flat)

  nc        <- ncol(flat)
  tree_col  <- flat[[nc]]
  seg_col   <- flat[[nc - 1]]    # 0-based segment index
  radius_col <- flat[[5]]         # element 5 (1-indexed) = radius
  err_col <- if (nc >= 6L) flat[[6]] else rep(NA_real_, nrow(flat))

  slice_key <- paste(as.integer(tid_sub), as.integer(seg_idx), sep = ":")
  n_pts_by_slice <- table(slice_key)
  out_key <- paste(as.integer(tree_col), as.integer(seg_col), sep = ":")
  n_pts <- as.integer(n_pts_by_slice[out_key])
  n_pts[is.na(n_pts)] <- NA_integer_

  fit_err <- ifelse(!is.na(err_col) & !is.na(n_pts) & n_pts > 0L,
                    sqrt(pmax(err_col, 0) / n_pts),
                    NA_real_)

  z_lo <- z_start + seg_col * step
  z_hi <- z_lo + dz

  out <- data.frame(
    TreeID     = as.integer(tree_col),
    Segment    = as.integer(seg_col) + 1L,
    Z_low      = z_lo,
    Z_high     = z_hi,
    Z_mid      = (z_lo + z_hi) / 2,
    Radius     = radius_col,
    Diameter   = 2 * radius_col,
    RANSAC_err = fit_err,
    Fit_method = fit_method_name,
    N_pts      = n_pts,
    Valid      = !is.na(radius_col) & radius_col > 0 & radius_col <= max_radius,
    stringsAsFactors = FALSE
  )
  out[out$Valid, , drop = FALSE]
}


# ----------------------------------------------------------------------------
#' Stem-segment method: RANSAC cylinder fit per height slice
#'
#' Fits a 3-D cylinder to each height slice using RANSAC.  Unlike the
#' circle-fit approach, cylinders capture stem lean and provide an axis
#' orientation.  Calls `cppCylinderFit()` with `method = "ransac"` per slice.
#'
#' @param dz              numeric.  Slice thickness (m).  Default `0.5`.
#' @param overlap         numeric \[0, 1).  Slice overlap fraction.  Default
#'   `0.1`.
#' @param n_ransac        integer.  Points sampled per RANSAC iteration.
#'   Default `10L`.
#' @param conf            numeric.  RANSAC confidence level.  Default `0.95`.
#' @param inlier_frac     numeric.  Expected inlier fraction.  Default `0.90`.
#' @param max_angle       numeric.  Max tolerated deviation from vertical
#'   (degrees).  Default `30`.
#' @param n_best          integer.  Number of best fits to aggregate.
#'   Default `20L`.
#' @param max_radius      numeric.  Max accepted cylinder radius (m).
#'   Default `1.0`.
#' @param min_pts         integer.  Min points in a slice to attempt a fit.
#'   Default `10L`.
#' @param use_stem_column logical.  Use the `Stem` column if present
#'   to restrict which points enter the fitting.  Default `TRUE`.
#' @param search_radius   numeric or `NULL`.  Max XY search radius (m) from
#'   the tree seed.  `NULL` = `max(Radius) × 3`.  Default `NULL`.
#'
#' @return A closure of class `"spanner_sgs_mtd"` / `"function"` with
#'   attribute `method_name = "ransac_cylinder"`.
#'
#' @seealso [stem_segments()], [seg_irls_cylinder()]
#'
#' @examples
#' \donttest{
#' segs <- stem_segments(las_s, map,
#'           method = seg_ransac_cylinder(dz = 0.5, max_angle = 20))
#' }
#'
#' @export
seg_ransac_cylinder <- function(dz              = 0.5,
                                 overlap         = 0.1,
                                 n_ransac        = 10L,
                                 conf            = 0.95,
                                 inlier_frac     = 0.90,
                                 max_angle       = 30,
                                 n_best          = 20L,
                                 max_radius      = 1.0,
                                 min_pts         = 10L,
                                 use_stem_column = TRUE,
                                 search_radius   = NULL) {

  func <- function(las, map) {
    if (!inherits(las, "LAS"))
      stop("seg_ransac_cylinder: 'las' must be a LAS object", call. = FALSE)
    if (!"treeID" %in% names(las@data))
      stop("seg_ransac_cylinder: 'las' must have a 'treeID' column — run tree_points() first.",
           call. = FALSE)
    .sgs_plot_cylinders(
      las, map, dz, overlap, max_radius, min_pts,
      use_stem_column, search_radius,
      "ransac_cylinder",
      function(mat, tids, segs, rads)
        C_ransac_plot_cylinders(mat, tids, segs, rads,
                                as.integer(n_ransac), conf, inlier_frac, 0.05)
    )
  }

  structure(func,
            class       = c("spanner_sgs_mtd", "function"),
            method_name = "ransac_cylinder")
}


# ----------------------------------------------------------------------------
#' Stem-segment method: IRLS cylinder fit per height slice
#'
#' Fits 3-D cylinders using Iteratively Reweighted Least Squares (IRLS).
#' Faster than RANSAC for clean, low-noise stem data; less robust to outliers
#' and partial occlusion.  Calls `cppCylinderFit()` with `method = "irls"`.
#'
#' @param dz              numeric.  Slice thickness (m).  Default `0.5`.
#' @param overlap         numeric \[0, 1).  Slice overlap fraction.  Default
#'   `0.1`.
#' @param max_radius      numeric.  Max accepted cylinder radius (m).
#'   Default `1.0`.
#' @param min_pts         integer.  Min points in a slice.  Default `10L`.
#' @param use_stem_column logical.  Default `TRUE`.
#' @param search_radius   numeric or `NULL`.  Default `NULL`.
#'
#' @return A closure of class `"spanner_sgs_mtd"` / `"function"` with
#'   attribute `method_name = "irls_cylinder"`.
#'
#' @seealso [stem_segments()], [seg_ransac_cylinder()]
#'
#' @examples
#' \donttest{
#' segs <- stem_segments(las_s, map, method = seg_irls_cylinder())
#' }
#'
#' @export
seg_irls_cylinder <- function(dz              = 0.5,
                               overlap         = 0.1,
                               max_radius      = 1.0,
                               min_pts         = 10L,
                               use_stem_column = TRUE,
                               search_radius   = NULL) {

  func <- function(las, map) {
    if (!inherits(las, "LAS"))
      stop("seg_irls_cylinder: 'las' must be a LAS object", call. = FALSE)
    if (!"treeID" %in% names(las@data))
      stop("seg_irls_cylinder: 'las' must have a 'treeID' column — run tree_points() first.",
           call. = FALSE)
    .sgs_plot_cylinders(
      las, map, dz, overlap, max_radius, min_pts,
      use_stem_column, search_radius,
      "irls_cylinder",
      function(mat, tids, segs, rads)
        C_irls_plot_cylinders(mat, tids, segs, rads, 100L, 0.05)
    )
  }

  structure(func,
            class       = c("spanner_sgs_mtd", "function"),
            method_name = "irls_cylinder")
}


# ----------------------------------------------------------------------------
#' Stem-segment method: brute-force RANSAC cylinder fit per height slice
#'
#' Fits 3-D cylinders using a brute-force orientation search combined with
#' RANSAC — the most robust option for strongly leaning or heavily occluded
#' stems.  Calls `cppCylinderFit()` with `method = "bf"`.
#'
#' @param dz              numeric.  Slice thickness (m).  Default `0.5`.
#' @param overlap         numeric \[0, 1).  Slice overlap fraction.  Default
#'   `0.1`.
#' @param n_ransac        integer.  RANSAC samples per iteration.
#'   Default `10L`.
#' @param conf            numeric.  RANSAC confidence level.  Default `0.95`.
#' @param inlier_frac     numeric.  Expected inlier fraction.  Default `0.90`.
#' @param max_angle       numeric.  Max tolerated deviation from vertical
#'   (degrees).  Default `30`.
#' @param n_best          integer.  Number of best fits to aggregate.
#'   Default `20L`.
#' @param max_radius      numeric.  Max accepted cylinder radius (m).
#'   Default `1.0`.
#' @param min_pts         integer.  Min slice points.  Default `10L`.
#' @param use_stem_column logical.  Default `TRUE`.
#' @param search_radius   numeric or `NULL`.  Default `NULL`.
#'
#' @return A closure of class `"spanner_sgs_mtd"` / `"function"` with
#'   attribute `method_name = "bf_cylinder"`.
#'
#' @seealso [stem_segments()], [seg_ransac_cylinder()], [seg_irls_cylinder()]
#'
#' @examples
#' \donttest{
#' segs <- stem_segments(las_s, map,
#'           method = seg_bf_cylinder(max_angle = 45, n_best = 25L))
#' }
#'
#' @export
seg_bf_cylinder <- function(dz              = 0.5,
                             overlap         = 0.1,
                             n_ransac        = 10L,
                             conf            = 0.95,
                             inlier_frac     = 0.90,
                             max_angle       = 30,
                             n_best          = 20L,
                             max_radius      = 1.0,
                             min_pts         = 10L,
                             use_stem_column = TRUE,
                             search_radius   = NULL) {

  fit_fn <- function(sl) {
    cppCylinderFit(sl, method   = "bf",
                       n        = as.integer(n_ransac),
                       p        = conf,
                       inliers  = inlier_frac,
                       max_angle= max_angle,
                       n_best   = as.integer(n_best))
  }

  func <- function(las, map) {
    if (!inherits(las, "LAS"))
      stop("seg_bf_cylinder: 'las' must be a LAS object", call. = FALSE)
    if (!"treeID" %in% names(las@data))
      stop("seg_bf_cylinder: 'las' must have a 'treeID' column — run tree_points() first.",
           call. = FALSE)
    .sgs_cylinder_loop(las, map, dz, overlap, max_radius, min_pts,
                        use_stem_column, search_radius,
                        "bf_cylinder", fit_fn)
  }

  structure(func,
            class       = c("spanner_sgs_mtd", "function"),
            method_name = "bf_cylinder")
}
