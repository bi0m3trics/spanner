# ============================================================================
# spanner – Bole_Segmentation.R
#
# Fit circles/cylinders up individual tree boles and compute stem volumes.
#
# Public functions:
#   segment_bole()        – per-tree, per-slice circle fitting → segment table
#   compute_bole_volume() – rollup segment table to per-tree volumes
# ============================================================================


# ----------------------------------------------------------------------------
#' Segment tree boles by fitting circles across height slices
#'
#' `segment_bole` divides each segmented tree's bole into overlapping height
#' slices, fits a 2-D RANSAC circle (or Hough-seeded RANSAC circle) to the
#' bole points in each slice, and returns a tidy `data.frame` with one row per
#' (TreeID × slice). The result feeds directly into [compute_bole_volume()]
#' and [assess_tree_quality()].
#'
#' @details
#' Three fitting algorithms are available:
#' \describe{
#'   \item{`"ransac"` (default)}{Calls the compiled C++ RANSAC circle fitter
#'     (`C_ransac_bole_slices`) for all trees simultaneously.
#'     No external dependencies. Fast and robust to outlier points.}
#'   \item{`"hough"`}{Uses `houghPrimitives::detect_tree_boles()` to label
#'     cylinder points per tree, then derives per-slice circle parameters from
#'     the labelled points. Requires the **houghPrimitives** package. Best
#'     for noisy data where pure RANSAC may drift.}
#'   \item{`"hough_seed_ransac"`}{Calls `"hough"` to obtain an initial centre
#'     and radius estimate per slice, then narrows the RANSAC search to points
#'     within `2 × hough_radius` of the Hough centre before refitting.
#'     Combines the global optimality of the Hough accumulator with the
#'     precision of RANSAC. Requires **houghPrimitives**.}
#' }
#'
#' If the `las` object has a `Bole` column (from [classify_lw_points()]) or a
#' legacy `Stem` column (from [classify_stem_points()])
#' and `use_stem_column = TRUE`, only bole/stem-labelled points are used for
#' fitting; otherwise all points within `search_radius` of the tree seed are
#' used.
#'
#' @param las LAS object. Must have a `treeID` column ([segment_graph()]
#'   output). Optionally has a `Bole` column ([classify_lw_points()] output) or
#'   a legacy `Stem` column ([classify_stem_points()] output).
#' @param tree_locs sf object from [get_raster_eigen_treelocs()].
#' @param algorithm character. Fitting algorithm: `"ransac"` (default),
#'   `"hough"`, or `"hough_seed_ransac"`.
#' @param z_min numeric. Lowest height (m above ground) to model. Default `0.1`.
#' @param z_max numeric. Highest height to model. `NULL` = use per-tree
#'   maximum Z. Default `NULL`.
#' @param dz numeric. Slice thickness (m). Default `0.5`.
#' @param overlap numeric \[0, 1). Fractional overlap between adjacent slices.
#'   Default `0.1`.
#' @param use_stem_column logical. Use the `Bole` column (preferred) or the
#'   legacy `Stem` column (if present) to
#'   restrict which points are fitted. Default `TRUE`.
#' @param search_radius numeric. Maximum XY distance (m) from tree axis for
#'   point inclusion when `Stem` column is absent or `use_stem_column = FALSE`.
#'   `NULL` = `max(tree_locs$Radius) * 3`. Default `NULL`.
#' @param inlier_tol numeric. Distance from the fitted circle surface (m) that
#'   counts as an inlier. Default `0.03`.
#' @param n_ransac integer. Points sampled per RANSAC iteration. Default `20L`.
#' @param conf numeric. RANSAC confidence level. Default `0.99`.
#' @param inlier_frac numeric. Expected inlier fraction. Default `0.85`.
#' @param n_best integer. Number of best RANSAC solutions to aggregate.
#'   Default `25L`.
#' @param max_radius numeric. Reject fitted circles with radius > this value
#'   (m). Default `1.0`.
#' @param min_pts integer. Minimum points in a slice to attempt a fit.
#'   Default `10L`.
#' @param hough_min_radius numeric. Minimum Hough cylinder radius (m).
#'   Default `0.03`.
#' @param hough_max_radius numeric. Maximum Hough cylinder radius (m).
#'   Default `0.60`.
#' @param hough_min_votes integer. Minimum Hough accumulator votes.
#'   Default `15L`.
#' @param voxel_thin numeric or `NULL`. If non-NULL, the LAS is voxel-thinned
#'   to this resolution (m) using [lidR::decimate_points()] with
#'   [lidR::barycenter_per_voxel()] before circle fitting. The barycenter
#'   method preserves sub-voxel geometry for bole cross-sections. Useful for
#'   very dense TLS to reduce C++ computation time. Default `NULL`.
#'
#' @return A `data.frame` with one row per (TreeID, Segment) and columns:
#'   \describe{
#'     \item{`TreeID`}{integer.}
#'     \item{`Segment`}{integer. 1 = lowest slice.}
#'     \item{`Z_low`, `Z_high`, `Z_mid`}{numeric. Slice height bounds/centre.}
#'     \item{`X_ctr`, `Y_ctr`}{numeric. Fitted circle centre.}
#'     \item{`Radius`}{numeric. Fitted circle radius (m).}
#'     \item{`Diameter`}{numeric. `2 × Radius` (m).}
#'     \item{`RANSAC_err`}{numeric. Mean inlier residual (m).}
#'     \item{`N_pts`}{integer. Points in slice.}
#'     \item{`N_inliers`}{integer. RANSAC inlier count.}
#'     \item{`Inlier_frac`}{numeric. `N_inliers / N_pts`.}
#'     \item{`Fit_method`}{character. Algorithm used.}
#'     \item{`Valid`}{logical. `FALSE` if too few points or radius out of range.}
#'   }
#'
#' @seealso [classify_stem_points()], [compute_bole_volume()],
#'   [assess_tree_quality()]
#'
#' @references
#' Fischler, M. A., & Bolles, R. C. (1981). Random sample consensus: a
#' paradigm for model fitting with applications to image analysis and
#' automated cartography. *Communications of the ACM*, 24(6), 381–395.
#' \doi{10.1145/358669.358692}
#'
#' @examples
#' \donttest{
#' set_lidr_threads(4)
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#'
#' tree_locs <- get_raster_eigen_treelocs(las, res = 0.025, pt_spacing = 0.0254,
#'                                        dens_threshold = 0.25,
#'                                        neigh_sizes = c(0.25, 0.15, 0.66),
#'                                        eigen_threshold = 0.75,
#'                                        grid_slice_min = 1, grid_slice_max = 2,
#'                                        minimum_polygon_area = 0.005,
#'                                        cylinder_fit_type = "ransac",
#'                                        max_dia = 1, SDvert = 0.33)
#'
#' las_seg  <- segment_graph(las, tree_locs, k = 50,
#'                           distance.threshold = 0.5,
#'                           use.metabolic.scale = FALSE,
#'                           ptcloud_slice_min = 0.5, ptcloud_slice_max = 2.5,
#'                           subsample.graph = 0.1, return.dense = TRUE)
#'
#' las_bole <- classify_stem_points(las_seg, tree_locs, method = "eigen")
#'
#' seg_tbl  <- segment_bole(las_bole, tree_locs)
#' head(seg_tbl)
#' }
#'
#' @export
segment_bole <- function(las,
                         tree_locs,
                         algorithm       = c("ransac", "hough",
                                             "hough_seed_ransac"),
                         z_min           = 0.1,
                         z_max           = NULL,
                         dz              = 0.5,
                         overlap         = 0.1,
                         use_stem_column = TRUE,
                         search_radius   = NULL,
                         inlier_tol      = 0.03,
                         n_ransac        = 20L,
                         conf            = 0.99,
                         inlier_frac     = 0.85,
                         n_best          = 25L,
                         max_radius      = 1.0,
                         min_pts         = 10L,
                         hough_min_radius = 0.03,
                         hough_max_radius = 0.60,
                         hough_min_votes  = 15L,
                         voxel_thin       = NULL,
                         ncpu             = 1L) {

  algorithm <- match.arg(algorithm)

  # ---- Validation -----------------------------------------------------------
  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object (output of segment_graph).")
  if (!"treeID" %in% names(las@data))
    stop("'las' must contain a 'treeID' column — run segment_graph() first.")
  if (!inherits(tree_locs, "sf"))
    stop("'tree_locs' must be an sf object from get_raster_eigen_treelocs().")

  if (algorithm %in% c("hough", "hough_seed_ransac")) {
    if (!requireNamespace("houghPrimitives", quietly = TRUE))
      stop("algorithm '", algorithm, "' requires the houghPrimitives package.\n",
           "Install with: devtools::install_github('bi0m3trics/houghPrimitives')\n",
           "Or use algorithm = 'ransac'.")
  }

  # ---- Geometry setup -------------------------------------------------------
  tl_xy  <- sf::st_coordinates(tree_locs)
  tl_ids <- as.integer(tree_locs$TreeID)

  if (is.null(search_radius)) {
    search_radius <- max(tree_locs$Radius, na.rm = TRUE) * 3.0
    search_radius <- max(search_radius, 0.20)
  }

  # Prefer "Bole" column (new API), fall back to legacy "Stem" column
  stem_col <- if ("Bole" %in% names(las@data)) "Bole"
              else if ("Stem" %in% names(las@data)) "Stem"
              else NULL
  has_stem <- use_stem_column && !is.null(stem_col)

  # ---- Optional voxel thinning (preserves bole geometry via barycenter) -----
  if (!is.null(voxel_thin) && voxel_thin > 0) {
    message("  [segment_bole] Voxel-thinning to barycenter (res = ",
            voxel_thin, " m) ...")
    las <- lidR::decimate_points(las,
                                 lidR::barycenter_per_voxel(res = voxel_thin))
  }

  # Compute per-tree z_max if not supplied
  if (is.null(z_max)) {
    z_max_use <- tapply(las@data$Z, las@data$treeID, max, na.rm = TRUE)
  } else {
    z_max_use <- setNames(rep(z_max, length(tl_ids)), tl_ids)
  }

  # ---- Choose algorithm path ------------------------------------------------
  if (algorithm == "ransac") {
    seg_tbl <- .segment_ransac(
      las, tl_xy, tl_ids, z_min, z_max_use,
      dz, overlap, search_radius, inlier_tol,
      n_ransac, conf, inlier_frac, n_best, max_radius, min_pts,
      has_stem, stem_col, as.integer(max(1L, ncpu))
    )
    seg_tbl$Fit_method <- "ransac"

  } else if (algorithm == "hough") {
    seg_tbl <- .segment_hough(
      las, tl_xy, tl_ids, z_min, z_max_use, dz, overlap,
      search_radius, inlier_tol,
      hough_min_radius, hough_max_radius, hough_min_votes,
      max_radius, min_pts, has_stem, seed_ransac = FALSE,
      n_ransac, conf, inlier_frac, n_best, stem_col
    )
    seg_tbl$Fit_method <- "hough"

  } else {  # hough_seed_ransac
    seg_tbl <- .segment_hough(
      las, tl_xy, tl_ids, z_min, z_max_use, dz, overlap,
      search_radius, inlier_tol,
      hough_min_radius, hough_max_radius, hough_min_votes,
      max_radius, min_pts, has_stem, seed_ransac = TRUE,
      n_ransac, conf, inlier_frac, n_best, stem_col
    )
    seg_tbl$Fit_method <- "hough_seed_ransac"
  }

  # ---- Add Diameter column --------------------------------------------------
  seg_tbl$Diameter <- seg_tbl$Radius * 2.0

  # Reorder columns for readability
  col_order <- c("TreeID", "Segment", "Z_low", "Z_high", "Z_mid",
                 "X_ctr", "Y_ctr", "Radius", "Diameter", "RANSAC_err",
                 "N_pts", "N_inliers", "Inlier_frac", "Fit_method", "Valid")
  col_order <- col_order[col_order %in% names(seg_tbl)]
  seg_tbl   <- seg_tbl[, col_order, drop = FALSE]

  return(seg_tbl)
}


# ============================================================================
# Internal helpers for segment_bole
# ============================================================================

# Determine global z_max per tree or use a scalar lookup
.per_tree_zmax <- function(tl_ids, z_max_use, z_min) {
  vapply(tl_ids, function(tid) {
    v <- z_max_use[as.character(tid)]
    if (is.na(v) || v <= z_min) z_min + 1.0 else v
  }, numeric(1))
}

# Build a bole-point mask (logical vector, same length as las@data rows)
.bole_mask <- function(las, tx, ty, z_lo, z_hi, sr2, has_stem, stem_col = NULL) {
  X  <- las@data$X
  Y  <- las@data$Y
  Z  <- las@data$Z
  mask <- Z >= z_lo & Z <= z_hi & (X - tx)^2 + (Y - ty)^2 <= sr2
  if (has_stem && !is.null(stem_col))
    mask <- mask & !is.na(las@data[[stem_col]]) & las@data[[stem_col]]
  mask
}

# ---------------------------------------------------------------------------
# RANSAC path: delegates to C++ C_ransac_bole_slices for all trees at once
# ---------------------------------------------------------------------------
.segment_ransac <- function(las, tl_xy, tl_ids,
                            z_min, z_max_use,
                            dz, overlap, search_radius, inlier_tol,
                            n_ransac, conf, inlier_frac, n_best,
                            max_radius, min_pts, has_stem,
                            stem_col = NULL, ncpu = 1L) {

  # Subset to bole-candidate points only (respects Stem flag)
  X  <- las@data$X
  Y  <- las@data$Y
  Z  <- las@data$Z
  TID <- as.integer(las@data$treeID)

  # If using stem column, zero out non-stem tree IDs so C++ ignores them
  if (has_stem && !is.null(stem_col)) {
    is_stem <- !is.na(las@data[[stem_col]]) & las@data[[stem_col]]
    TID[!is_stem] <- 0L
  }

  # Per-tree z_max: use the maximum of the per-tree z_max values
  global_z_max <- max(unlist(z_max_use), na.rm = TRUE)
  global_z_max <- max(global_z_max, z_min + dz)

  result <- C_ransac_bole_slices(
    X          = X,
    Y          = Y,
    Z          = Z,
    treeIds    = TID,
    tree_X     = tl_xy[, 1],
    tree_Y     = tl_xy[, 2],
    tree_ID_vals = tl_ids,
    z_min      = z_min,
    z_max      = global_z_max,
    dz         = dz,
    overlap    = overlap,
    search_radius = search_radius,
    inlier_tol = inlier_tol,
    n_samples  = as.integer(n_ransac),
    confidence = conf,
    inlier_frac = inlier_frac,
    n_best     = as.integer(n_best),
    max_radius = max_radius,
    min_pts    = as.integer(min_pts),
    ncpu       = as.integer(ncpu)
  )

  # Filter out slices above each tree's own z_max
  z_max_lookup <- setNames(
    .per_tree_zmax(tl_ids, z_max_use, z_min),
    as.character(tl_ids)
  )
  keep <- vapply(seq_len(nrow(result)), function(i) {
    tmax <- z_max_lookup[as.character(result$TreeID[i])]
    if (is.na(tmax)) return(TRUE)
    result$Z_mid[i] <= tmax
  }, logical(1))

  result <- result[keep, , drop = FALSE]
  rownames(result) <- NULL
  return(result)
}


# ---------------------------------------------------------------------------
# Hough (+ optional RANSAC seed) path: per-tree, per-slice in R
# ---------------------------------------------------------------------------
.segment_hough <- function(las, tl_xy, tl_ids,
                           z_min, z_max_use, dz, overlap,
                           search_radius, inlier_tol,
                           hmin_r, hmax_r, hmin_v,
                           max_radius, min_pts, has_stem,
                           seed_ransac,
                           n_ransac, conf, inlier_frac, n_best,
                           stem_col = NULL) {

  X   <- las@data$X
  Y   <- las@data$Y
  Z   <- las@data$Z
  TID <- as.integer(las@data$treeID)
  sr2 <- search_radius^2

  step <- dz * max(0.01, 1.0 - overlap)

  out_rows <- list()

  for (ti in seq_along(tl_ids)) {
    tid  <- tl_ids[ti]
    tx   <- tl_xy[ti, 1]
    ty   <- tl_xy[ti, 2]
    zmax <- z_max_use[as.character(tid)]
    if (is.na(zmax) || zmax <= z_min) zmax <- z_min + 1.0

    n_slices <- max(1L, ceiling((zmax - z_min) / step))

    # Extract all bole-candidate points for this tree
    tree_mask <- TID == tid &
      Z >= z_min & Z <= zmax &
      (X - tx)^2 + (Y - ty)^2 <= sr2
    if (has_stem && !is.null(stem_col))
      tree_mask <- tree_mask & !is.na(las@data[[stem_col]]) & las@data[[stem_col]]

    tree_idx <- which(tree_mask)
    if (length(tree_idx) < min_pts) next

    tX <- X[tree_idx]; tY <- Y[tree_idx]; tZ <- Z[tree_idx]

    # Run Hough on the full bole column once per tree
    pts_mat <- cbind(tX, tY, tZ)
    hough_res <- tryCatch(
      houghPrimitives::detect_tree_boles(
        pts_mat,
        min_radius   = hmin_r,
        max_radius   = hmax_r,
        min_votes    = hmin_v,
        label_points = TRUE
      ),
      error = function(e) NULL
    )

    # Extract per-point shape_type for cylinder labelling
    hough_cyl <- rep(FALSE, length(tree_idx))
    if (!is.null(hough_res)) {
      st <- NULL
      if (is.data.frame(hough_res) && "shape_type" %in% names(hough_res)) {
        st <- hough_res$shape_type
      } else if (is.list(hough_res)) {
        if (!is.null(hough_res$labels) &&
            "shape_type" %in% names(hough_res$labels))
          st <- hough_res$labels$shape_type
        else if ("shape_type" %in% names(hough_res))
          st <- hough_res$shape_type
      }
      if (!is.null(st) && length(st) == length(tree_idx))
        hough_cyl <- st == 4L
    }

    # Per-slice fitting
    for (s in seq_len(n_slices)) {
      zlo  <- z_min + (s - 1L) * step
      zhi  <- zlo + dz
      zmid <- 0.5 * (zlo + zhi)
      if (zmid > zmax) next

      s_mask <- tZ >= zlo & tZ <= zhi
      npts   <- sum(s_mask)

      row <- list(
        TreeID = tid, Segment = s,
        Z_low = zlo, Z_high = zhi, Z_mid = zmid,
        X_ctr = NA_real_, Y_ctr = NA_real_, Radius = NA_real_,
        RANSAC_err = NA_real_, N_pts = npts,
        N_inliers = 0L, Inlier_frac = NA_real_, Valid = FALSE
      )

      if (npts < min_pts) {
        out_rows[[length(out_rows) + 1L]] <- row
        next
      }

      sX <- tX[s_mask]; sY <- tY[s_mask]; sZ <- tZ[s_mask]
      hc <- hough_cyl[s_mask]

      # --- Hough-based centre from cylinder-labelled points in this slice ---
      hough_cx <- NA_real_; hough_cy <- NA_real_; hough_r <- NA_real_
      if (sum(hc) >= 3L) {
        hough_cx <- mean(sX[hc])
        hough_cy <- mean(sY[hc])
        hough_r  <- median(sqrt((sX[hc] - hough_cx)^2 +
                                  (sY[hc] - hough_cy)^2))
      }

      # --- Fit circle to the slice -----------------------------------------
      if (seed_ransac && !is.na(hough_r) && hough_r > 0 &&
          hough_r <= max_radius) {
        # Narrow RANSAC to points near the Hough circle
        band    <- hough_r * 2.0
        close   <- sqrt((sX - hough_cx)^2 + (sY - hough_cy)^2) <= band
        fit_pts <- list(X = sX[close], Y = sY[close])
        use_ransac <- (sum(close) >= min_pts)
      } else {
        fit_pts    <- list(X = sX, Y = sY)
        use_ransac <- TRUE
      }

      if (use_ransac && length(fit_pts$X) >= min_pts) {
        # Use cylinderFit (irls/ransac method="bf" with max_angle=0 ≈ circle)
        pt_mat <- cbind(fit_pts$X, fit_pts$Y,
                        rep(zmid, length(fit_pts$X)))
        fit <- tryCatch(
          cylinderFit(lidR::LAS(data.table::data.table(
            X = fit_pts$X, Y = fit_pts$Y, Z = rep(zmid, length(fit_pts$X))
          )), method = "ransac", n = n_ransac, inliers = inlier_frac,
          conf = conf, n_best = n_best),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          # cylinderFit returns rho,theta,phi,alpha,radius,err,px,py,pz
          # For near-vertical cylinder px,py ≈ cx,cy
          cx <- if ("px" %in% names(fit)) fit$px else hough_cx
          cy <- if ("py" %in% names(fit)) fit$py else hough_cy
          r  <- if ("radius" %in% names(fit)) fit$radius else NA_real_
          er <- if ("err"    %in% names(fit)) fit$err    else NA_real_
          if (!is.na(r) && r > 0 && r <= max_radius) {
            n_in <- sum(abs(sqrt((sX - cx)^2 + (sY - cy)^2) - r) <= inlier_tol)
            row  <- list(
              TreeID = tid, Segment = s,
              Z_low = zlo, Z_high = zhi, Z_mid = zmid,
              X_ctr = cx, Y_ctr = cy, Radius = r,
              RANSAC_err = er, N_pts = npts,
              N_inliers = n_in, Inlier_frac = n_in / npts, Valid = TRUE
            )
          }
        }
      } else if (!is.na(hough_r) && hough_r > 0 && hough_r <= max_radius) {
        # Fall back to raw Hough estimate
        n_in <- sum(abs(sqrt((sX - hough_cx)^2 +
                               (sY - hough_cy)^2) - hough_r) <= inlier_tol)
        row  <- list(
          TreeID = tid, Segment = s,
          Z_low = zlo, Z_high = zhi, Z_mid = zmid,
          X_ctr = hough_cx, Y_ctr = hough_cy, Radius = hough_r,
          RANSAC_err = mean(abs(sqrt((sX - hough_cx)^2 +
                                      (sY - hough_cy)^2) - hough_r)),
          N_pts = npts, N_inliers = n_in,
          Inlier_frac = n_in / npts, Valid = TRUE
        )
      }

      out_rows[[length(out_rows) + 1L]] <- row
    }  # end slice loop
  }  # end tree loop

  if (length(out_rows) == 0L) {
    return(data.frame(
      TreeID = integer(), Segment = integer(),
      Z_low = numeric(), Z_high = numeric(), Z_mid = numeric(),
      X_ctr = numeric(), Y_ctr = numeric(), Radius = numeric(),
      RANSAC_err = numeric(), N_pts = integer(),
      N_inliers = integer(), Inlier_frac = numeric(), Valid = logical()
    ))
  }

  result <- do.call(rbind, lapply(out_rows, as.data.frame))
  rownames(result) <- NULL
  return(result)
}


# ----------------------------------------------------------------------------
#' Compute outside-bark bole volume from a segment table
#'
#' `compute_bole_volume` aggregates the per-slice circle fits returned by
#' [segment_bole()] into per-tree outside-bark stem volume estimates. Two
#' volume formulae are available: the simple cylinder formula and Smalian's
#' formula for tapering stems.
#'
#' @param seg_table `data.frame` returned by [segment_bole()].
#' @param method character. `"cylinder"` (default) uses `π r² h` per segment.
#'   `"smalian"` applies Smalian's formula between adjacent valid segments:
#'   `V = h × π × (r₁² + r₂²) / 2`.
#'
#' @return A `data.frame` with one row per tree and columns:
#'   \describe{
#'     \item{`TreeID`}{integer.}
#'     \item{`Volume_m3`}{numeric. Total outside-bark stem volume (m³).}
#'     \item{`N_segments`}{integer. Number of valid segments used.}
#'     \item{`Mean_radius`}{numeric. Mean fitted radius (m).}
#'     \item{`SD_radius`}{numeric. SD of fitted radii (m) — measure of taper.}
#'     \item{`Basal_area_m2`}{numeric. Basal area at the first valid segment.}
#'     \item{`Height_modeled`}{numeric. Vertical height span modeled (m).}
#'     \item{`Mean_RANSAC_err`}{numeric. Mean per-segment RANSAC residual (m).}
#'   }
#'
#' @seealso [segment_bole()], [assess_tree_quality()]
#'
#' @examples
#' \donttest{
#' # (continuing from segment_bole example)
#' # seg_tbl  <- segment_bole(las_bole, tree_locs)
#' # vol_tbl  <- compute_bole_volume(seg_tbl)
#' }
#'
#' @export
compute_bole_volume <- function(seg_table,
                                method = c("cylinder", "smalian")) {

  method <- match.arg(method)

  if (!is.data.frame(seg_table))
    stop("'seg_table' must be a data.frame returned by segment_bole().")

  required_cols <- c("TreeID", "Segment", "Z_low", "Z_high", "Radius", "Valid")
  missing_cols  <- setdiff(required_cols, names(seg_table))
  if (length(missing_cols) > 0)
    stop("'seg_table' is missing columns: ", paste(missing_cols, collapse = ", "),
         "\nMake sure it is the output of segment_bole().")

  tree_ids <- sort(unique(seg_table$TreeID))
  out_rows <- vector("list", length(tree_ids))

  for (i in seq_along(tree_ids)) {
    tid <- tree_ids[i]
    tr  <- seg_table[seg_table$TreeID == tid & !is.na(seg_table$Valid) &
                       seg_table$Valid, , drop = FALSE]

    n_seg <- nrow(tr)

    if (n_seg == 0L) {
      out_rows[[i]] <- data.frame(
        TreeID = tid, Volume_m3 = NA_real_, N_segments = 0L,
        Mean_radius = NA_real_, SD_radius = NA_real_,
        Basal_area_m2 = NA_real_, Height_modeled = NA_real_,
        Mean_RANSAC_err = NA_real_
      )
      next
    }

    # Sort by segment index
    tr <- tr[order(tr$Segment), ]
    r  <- tr$Radius
    h  <- tr$Z_high - tr$Z_low

    if (method == "cylinder") {
      vols <- pi * r^2 * h
      total_vol <- sum(vols, na.rm = TRUE)
    } else {
      # Smalian: between consecutive valid segments
      if (n_seg == 1L) {
        total_vol <- pi * r[1L]^2 * h[1L]
      } else {
        smalian_vols <- numeric(n_seg - 1L)
        for (k in seq_len(n_seg - 1L)) {
          smalian_vols[k] <-
            (tr$Z_high[k + 1L] - tr$Z_low[k]) *
            pi * (r[k]^2 + r[k + 1L]^2) / 2.0
        }
        # Add the top cone cap
        total_vol <- sum(smalian_vols, na.rm = TRUE) +
          pi * r[n_seg]^2 * h[n_seg] / 2.0
      }
    }

    ht_span <- max(tr$Z_high, na.rm = TRUE) - min(tr$Z_low, na.rm = TRUE)
    ba      <- pi * r[1L]^2

    mean_err <- if ("RANSAC_err" %in% names(tr))
      mean(tr$RANSAC_err, na.rm = TRUE) else NA_real_

    out_rows[[i]] <- data.frame(
      TreeID          = tid,
      Volume_m3       = total_vol,
      N_segments      = n_seg,
      Mean_radius     = mean(r, na.rm = TRUE),
      SD_radius       = if (n_seg > 1L) stats::sd(r, na.rm = TRUE) else 0.0,
      Basal_area_m2   = ba,
      Height_modeled  = ht_span,
      Mean_RANSAC_err = mean_err
    )
  }

  result <- do.call(rbind, out_rows)
  rownames(result) <- NULL
  return(result)
}
