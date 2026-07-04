# ============================================================================
# spanner – Stem_Classification.R
#
# Stem-point classification for individually segmented trees.
# This is the first step in the stem-segmentation pipeline that follows
# individual tree segmentation (segment_graph()).
#
# Public function:  classify_stem_points()
# Internal helpers: .stem_eigen(), .stem_lewos(), .stem_hough()
# ============================================================================


# ----------------------------------------------------------------------------
#' Classify stem points within a segmented TLS point cloud
#'
#' @description
#' **Deprecated in spanner v2.0 — use [stem_eigen()] or [stem_hough()] inside
#' [stem_points()] instead.**
#'
#' `classify_stem_points` labels every point in a segmented LAS object as
#' stem (wood stem) or non-stem and attaches two new columns—`Stem` (logical)
#' and `WoodProb` (numeric, 0–1)—to the returned LAS. It is designed to
#' operate on the output of [segment_graph()] and feeds directly into
#' [segment_stem()].
#'
#' @details
#' Four classification methods are available:
#' \describe{
#'   \item{`"eigen"`}{Computes local PCA eigenvalue metrics (Linearity and
#'     Verticality) via [eigen_metrics()] for each point. Points whose
#'     Linearity **≥** `linearity_threshold` **AND** Verticality **≥**
#'     `verticality_threshold` within `search_radius` of their tree's XY
#'     centroid are labelled as stem. Fastest option; works well for clean TLS
#'     data. Mirrors the `stm.eigen.knn` approach in TreeLS.}
#'   \item{`"lewos"`}{Like `"eigen"` but adds *k*-NN majority-vote label
#'     propagation (`n_propagation` iterations, `k_propagation` neighbours)
#'     inspired by the leaf–wood separation framework of Wang et al. (2020).
#'     Suitable for mixed leaf-wood scenes where isolated wood points may be
#'     missed by thresholding alone.}
#'   \item{`"hough"`}{Uses `houghPrimitives::detect_tree_boles()` to identify
#'     cylinder-classified points (shape type 4) per tree. Requires the
#'     **houghPrimitives** package. Robust to noise and partial occlusion.}
#'   \item{`"combined"`}{Runs `"eigen"` as a fast pre-filter then `"hough"` on
#'     the eigen-filtered candidates (logical AND). Requires **houghPrimitives**.
#'     Highest precision; recommended default for clean TLS data.}
#' }
#'
#' If `eigen_metrics` columns (`Linearity`, `Verticality`) are already present
#' in `las@data` they are reused without recomputation.
#'
#' @param las LAS object. Must contain a `treeID` column (output of
#'   [segment_graph()]).
#' @param tree_locs sf object of tree seed locations (output of
#'   [get_raster_eigen_treelocs()]).
#' @param method character. One of `"combined"` (default), `"eigen"`,
#'   `"lewos"`, or `"hough"`.
#' @param z_min numeric. Lower Z bound (m above ground) of the stem search
#'   window. Default `0.1`.
#' @param z_max numeric. Upper Z bound (m). Default `8.0`.
#' @param search_radius numeric. Maximum XY distance (m) from tree seed XY to
#'   include a point as a stem candidate. Defaults to `max(tree_locs$Radius) *
#'   3` when `NULL`.
#' @param neigh_k integer. Number of nearest neighbours (excluding the query
#'   point itself) for the local PCA neighbourhood used in `"eigen"` /
#'   `"lewos"`. Default `30L`.
#' @param linearity_threshold numeric \[0, 1\]. Minimum Linearity for a point to
#'   be considered wood. Default `0.65`.
#' @param verticality_threshold numeric \[0, 1\]. Minimum Verticality for a
#'   point to be considered wood. Default `0.70`.
#' @param n_propagation integer. Label-propagation iterations (`"lewos"`
#'   method). Default `3L`.
#' @param k_propagation integer. Number of neighbours for propagation (`"lewos"`
#'   method). Default `10L`.
#' @param hough_min_radius numeric. Minimum cylinder radius (m) for Hough
#'   detection. Default `0.03`.
#' @param hough_max_radius numeric. Maximum cylinder radius (m). Default `0.60`.
#' @param hough_min_votes integer. Minimum Hough accumulator votes to accept a
#'   detection. Default `20L`.
#' @param ncpu integer. CPU cores for parallel eigen computation. Default `4L`.
#' @param voxel_thin numeric or `NULL`. If non-NULL, the point cloud is
#'   voxel-thinned to this resolution (m) with [lidR::decimate_points()] and
#'   [lidR::random_per_voxel()] **before** [eigen_metrics()] is run. The
#'   resulting Linearity/Verticality labels are then projected back to the full
#'   cloud via 1-NN. Dramatically reduces computation time for very dense TLS.
#'   Default `NULL` (no thinning).
#'
#' @return A LAS object identical to `las` with two extra columns:
#'   \describe{
#'     \item{`Stem`}{logical — `TRUE` for stem-classified points.}
#'     \item{`WoodProb`}{numeric \[0, 1\] — classification confidence.}
#'   }
#'
#' @seealso [segment_graph()], [segment_stem()], [get_raster_eigen_treelocs()]
#'
#' @references
#' Wang, D., Momo Takoudjou, S., & Casella, E. (2020). LeWoS: A universal
#' leaf-wood classification method to facilitate the 3D modelling of large
#' tropical trees using terrestrial LiDAR. *Methods in Ecology and Evolution*,
#' 11(3), 376–389. \doi{10.1111/2041-210X.13342}
#'
#' de Conto, T., Olofsson, K., Gorgens, E. B., Rodriguez, L. C. E., &
#' Almeida, G. (2017). Performance of terrestrial and airborne LiDAR for
#' estimating stem volume and aboveground biomass of individual Pinus taeda
#' trees. *Computers and Electronics in Agriculture*, 140, 249–255.
#' \doi{10.1016/j.compag.2017.07.019}
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
#' las_seg <- segment_graph(las, tree_locs, k = 50,
#'                          distance.threshold = 0.5,
#'                          use.metabolic.scale = FALSE,
#'                          ptcloud_slice_min = 0.5, ptcloud_slice_max = 2.5,
#'                          subsample.graph = 0.1, return.dense = TRUE)
#'
#' # Classify stem points (eigen method — no houghPrimitives needed)
#' las_stem <- classify_stem_points(las_seg, tree_locs, method = "eigen")
#' plot(las_stem, color = "Stem")
#' }
#'
#' @export
classify_stem_points <- function(las,
                                 tree_locs,
                                 method              = c("combined", "eigen",
                                                         "lewos", "hough"),
                                 z_min               = 0.1,
                                 z_max               = 8.0,
                                 search_radius       = NULL,
                                 neigh_k             = 30L,
                                 linearity_threshold  = 0.65,
                                 verticality_threshold = 0.70,
                                 n_propagation       = 3L,
                                 k_propagation       = 10L,
                                 hough_min_radius    = 0.03,
                                 hough_max_radius    = 0.60,
                                 hough_min_votes     = 20L,
                                 ncpu                = 4L,
                                 voxel_thin          = NULL) {

  method <- match.arg(method)

  # ---- Validation -----------------------------------------------------------
  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object — run segment_graph() first.")
  if (!"treeID" %in% names(las@data))
    stop("'las' must contain a 'treeID' column — run segment_graph() first.")
  if (!inherits(tree_locs, "sf"))
    stop("'tree_locs' must be an sf object from get_raster_eigen_treelocs().")
  if (!"TreeID" %in% names(tree_locs))
    stop("'tree_locs' must have a 'TreeID' column.")

  if (method %in% c("hough", "combined")) {
    if (!requireNamespace("houghPrimitives", quietly = TRUE))
      stop("Method '", method, "' requires the houghPrimitives package.\n",
           "Install with: devtools::install_github('bi0m3trics/houghPrimitives')\n",
           "Or use method = 'eigen' or method = 'lewos'.")
  }

  # ---- Tree seed geometry ---------------------------------------------------
  tl_xy  <- sf::st_coordinates(tree_locs)   # n × 2 matrix (X, Y)
  tl_ids <- tree_locs$TreeID

  if (is.null(search_radius)) {
    search_radius <- max(tree_locs$Radius, na.rm = TRUE) * 3.0
    search_radius <- max(search_radius, 0.20)            # floor 20 cm
  }

  n_pts <- nrow(las@data)

  # ---- Method dispatch ------------------------------------------------------
  stem_flag <- rep(FALSE, n_pts)
  wood_prob  <- rep(0.0,   n_pts)

  if (method %in% c("eigen", "lewos", "combined")) {
    stem_flag <- .stem_eigen(
      las, tl_xy, tl_ids,
      z_min, z_max, search_radius,
      neigh_k, linearity_threshold, verticality_threshold, ncpu,
      propagate  = (method == "lewos"),
      n_prop     = as.integer(n_propagation),
      k_prop     = as.integer(k_propagation),
      voxel_thin = voxel_thin
    )
    wood_prob <- as.numeric(stem_flag)
  }

  if (method %in% c("hough", "combined")) {
    hough_flag <- .stem_hough(
      las, tl_xy, tl_ids,
      z_min, z_max, search_radius,
      hough_min_radius, hough_max_radius, as.integer(hough_min_votes)
    )
    if (method == "combined") {
      # Intersection of eigen + hough; confidence is average of both evidence
      wood_prob  <- (as.numeric(stem_flag) * 0.5 +
                       as.numeric(hough_flag) * 0.5)
      stem_flag  <- stem_flag & hough_flag
    } else {
      # Pure hough
      stem_flag <- hough_flag
      wood_prob  <- as.numeric(hough_flag)
    }
  }

  # ---- Attach columns to LAS ------------------------------------------------
  las@data[["Stem"]]     <- stem_flag
  las@data[["WoodProb"]] <- wood_prob

  return(las)
}


# ============================================================================
# Internal helpers
# ============================================================================

# ----------------------------------------------------------------------------
# .stem_eigen
# Linearity + Verticality threshold, optionally followed by kNN propagation.
# ----------------------------------------------------------------------------
.stem_eigen <- function(las, tl_xy, tl_ids,
                        z_min, z_max, search_radius,
                        neigh_k, lin_thr, vert_thr, ncpu,
                        propagate, n_prop, k_prop,
                        voxel_thin = NULL) {

  # Compute eigen metrics, optionally on a voxel-thinned cloud
  need_eigen <- !all(c("Linearity", "Verticality") %in% names(las@data))

  if (!is.null(voxel_thin) && voxel_thin > 0) {
    # Always recompute on thinned cloud (overrides cached columns)
    message("  [classify_stem_points] Voxel-thinning (res = ", voxel_thin,
            " m) before eigen computation ...")
    las_thin <- lidR::decimate_points(las,
                                      lidR::random_per_voxel(res = voxel_thin,
                                                             n   = 1L))
    em <- eigen_metrics(las_thin, k = neigh_k, ncpu = ncpu)
    # Transfer eigen labels from thinned to full cloud via 1-NN (cross-cloud)
    nn1 <- FNN::knnx.index(
      data  = cbind(las_thin@data$X, las_thin@data$Y, las_thin@data$Z),
      query = cbind(las@data$X,      las@data$Y,      las@data$Z),
      k     = 1L
    )[, 1L]
    lin_vals  <- em$Linearity[nn1]
    vert_vals <- em$Verticality[nn1]
    # Replace any NAs (edge points with too few thin-cloud neighbours) with 0
    lin_vals[is.na(lin_vals)]   <- 0
    vert_vals[is.na(vert_vals)] <- 0
    las@data$Linearity   <- lin_vals
    las@data$Verticality <- vert_vals
  } else if (need_eigen) {
    message("  [classify_stem_points] Computing eigen metrics",
            " (k = ", neigh_k, " neighbours) ...")
    em <- eigen_metrics(las, k = neigh_k, ncpu = ncpu)
    las@data$Linearity   <- em$Linearity
    las@data$Verticality <- em$Verticality
  }

  X      <- las@data$X
  Y      <- las@data$Y
  Z      <- las@data$Z
  tid_v  <- las@data$treeID
  lin    <- las@data$Linearity
  vert   <- las@data$Verticality
  n_pts  <- length(X)

  stem_flag <- rep(FALSE, n_pts)

  for (ti in seq_along(tl_ids)) {
    tid <- tl_ids[ti]
    tx  <- tl_xy[ti, 1]
    ty  <- tl_xy[ti, 2]
    sr2 <- search_radius^2

    idx <- which(
      tid_v == tid &
        Z >= z_min & Z <= z_max &
        (X - tx)^2 + (Y - ty)^2 <= sr2 &
        !is.na(lin) & !is.na(vert)
    )
    if (length(idx) == 0L) next

    wood <- lin[idx] >= lin_thr & vert[idx] >= vert_thr
    stem_flag[idx[wood]] <- TRUE
  }

  # ---- Optional kNN label propagation (LeWoS-inspired) --------------------
  if (propagate && any(stem_flag)) {
    stem_idx <- which(
      Z >= z_min & Z <= z_max &
        !is.na(lin) & !is.na(vert) &
        tid_v > 0
    )
    n_stem <- length(stem_idx)

    if (n_stem > k_prop + 1L) {
      coords     <- cbind(X[stem_idx], Y[stem_idx], Z[stem_idx])
      labels_sub <- stem_flag[stem_idx]

      # FNN::knn.index excludes self when test=NULL, so k = k_prop gives
      # exactly k true neighbours (same semantics as lidR::knn)
      nn_idx <- FNN::knn.index(coords, k = k_prop)

      for (iter in seq_len(n_prop)) {
        vote_mat   <- matrix(labels_sub[nn_idx], nrow = n_stem, ncol = k_prop)
        vote_count <- rowSums(vote_mat)
        labels_sub <- vote_count > (k_prop / 2L)
      }
      stem_flag[stem_idx] <- labels_sub
    }
  }

  return(stem_flag)
}


# ----------------------------------------------------------------------------
# .stem_hough
# Per-tree Hough Transform cylinder detection via houghPrimitives.
# ----------------------------------------------------------------------------
.stem_hough <- function(las, tl_xy, tl_ids,
                        z_min, z_max, search_radius,
                        min_r, max_r, min_votes) {

  X      <- las@data$X
  Y      <- las@data$Y
  Z      <- las@data$Z
  tid_v  <- las@data$treeID
  n_pts  <- length(X)
  sr2    <- search_radius^2

  stem_flag <- rep(FALSE, n_pts)

  for (ti in seq_along(tl_ids)) {
    tid <- tl_ids[ti]
    tx  <- tl_xy[ti, 1]
    ty  <- tl_xy[ti, 2]

    idx <- which(
      tid_v == tid &
        Z >= z_min & Z <= z_max &
        (X - tx)^2 + (Y - ty)^2 <= sr2
    )
    if (length(idx) < 20L) next  # need a minimum density for Hough

    pts_mat <- cbind(X[idx], Y[idx], Z[idx])

    result <- tryCatch(
      houghPrimitives::detect_tree_boles(
        pts_mat,
        min_radius   = min_r,
        max_radius   = max_r,
        min_votes    = min_votes,
        label_points = TRUE
      ),
      error = function(e) {
        message("  [classify_stem_points] Hough failed for TreeID ", tid,
                ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(result)) next

    # Robustly extract per-point shape_type regardless of return structure
    shape_type <- NULL
    if (is.data.frame(result) && "shape_type" %in% names(result)) {
      shape_type <- result$shape_type
    } else if (is.list(result)) {
      if (!is.null(result$labels) && "shape_type" %in% names(result$labels))
        shape_type <- result$labels$shape_type
      else if ("shape_type" %in% names(result))
        shape_type <- result$shape_type
    }

    if (!is.null(shape_type) && length(shape_type) == length(idx)) {
      stem_flag[idx[shape_type == 4L]] <- TRUE
    }
  }

  return(stem_flag)
}
