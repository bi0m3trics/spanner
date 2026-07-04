# ============================================================================
# spanner – LW_Classification.R
#
# Combined leaf–wood classification for individually segmented TLS trees.
# Computes eigen_metrics ONCE at a single neighbourhood size and reuses the
# result for both stem detection (Linearity + Verticality + LeWoS propagation)
# and branch detection (Curvature + Sphericity thresholds).
#
# Public function:  classify_lw_points()
# ============================================================================


# ----------------------------------------------------------------------------
#' Classify leaf–wood components in a segmented TLS point cloud
#'
#' `classify_lw_points` labels every point as **Stem**, **Branch**, or
#' **Other** and attaches five new columns — `Stem`, `StemProb`, `Branch`,
#' `BranchProb`, `Other` — to the returned LAS.
#'
#' A single call to [eigen_metrics()] is made at `neigh_k` neighbours.  The
#' resulting metrics drive **both** stem detection (Linearity ≥
#' `linearity_threshold` **AND** Verticality ≥ `verticality_threshold`,
#' followed by optional LeWoS kNN label propagation) **and** branch detection
#' (Curvature < `curvature_threshold` **AND** Sphericity <
#' `sphericity_threshold` on non-stem points).
#'
#' Stem detection logic mirrors [classify_stem_points()] with `method =
#' "lewos"`.
#'
#' Branch probability (`BranchProb`) is a continuous 0–1 affinity score for
#' all non-stem, non-ground points:
#' \deqn{\text{BranchProb} = \frac{1}{2}\left[\max\!\left(0,\,1 -
#'   \frac{\text{Curvature}}{\text{curvature\_threshold}}\right) +
#'   \max\!\left(0,\,1 -
#'   \frac{\text{Sphericity}}{\text{sphericity\_threshold}}\right)\right]}
#' Values near 1 indicate a smooth, elongated cylinder (large structural
#' branch); values near 0 indicate foliage or fine irregular structure.
#'
#' @param las LAS object.  Must contain a `treeID` column (output of
#'   [segment_graph()]).
#' @param tree_locs sf object of tree seed locations (output of
#'   [get_raster_eigen_treelocs()]).
#' @param method character.  `"lewos"` (default) adds kNN majority-vote label
#'   propagation after the initial Linearity/Verticality threshold filter.
#'   `"eigen"` applies the threshold only.
#' @param z_min numeric. Lower Z bound (m above ground) of the stem search
#'   window.  Default `0.1`.
#' @param z_max numeric. Upper Z bound (m).  Default `50.0`.
#' @param search_radius numeric. Maximum XY distance (m) from the tree seed XY
#'   centroid for stem candidacy.  `NULL` → `max(tree_locs$Radius) * 3`
#'   (floor 0.20 m).
#' @param axis_refine_stems logical. If `TRUE`, stem candidates are first found
#'   from Linearity/Verticality within each assigned `treeID`, then filtered
#'   around local per-tree stem axes inferred from those candidates. Default
#'   `TRUE`.
#' @param axis_search_radius numeric or `NULL`. Radius around inferred local
#'   stem axes. `NULL` uses `search_radius`. Default `NULL`.
#' @param axis_min_points integer. Minimum candidate points needed to infer
#'   local stem axes. Trees below this use the seed-centred fallback. Default
#'   `10L`.
#' @param neigh_k integer. Number of nearest neighbours for **all** eigen
#'   metric computations (stem and branch share this value).  Default `20L`.
#' @param linearity_threshold numeric \[0, 1\]. Minimum Linearity for stem
#'   candidacy.  Default `0.10`.
#' @param verticality_threshold numeric \[0, 1\]. Minimum Verticality for stem
#'   candidacy.  Default `0.75`.
#' @param n_propagation integer. LeWoS propagation iterations.  Default `3L`.
#' @param k_propagation integer. Neighbours for propagation voting.
#'   Default `15L`.
#' @param curvature_threshold numeric. Points with Curvature **below** this
#'   value are branch candidates.  Default `0.15` (derived from PCA/LDA
#'   analysis on a TLS snag).
#' @param sphericity_threshold numeric. Points with Sphericity **below** this
#'   value are branch candidates.  Default `0.29` (as above).
#' @param ncpu integer. CPU cores for parallel eigen computation.
#'   Default `4L`.
#' @param voxel_thin numeric or `NULL`.  If non-NULL, the cloud is voxel-
#'   thinned to this resolution (m) before [eigen_metrics()] and all four
#'   metrics are projected back to the full cloud via 1-NN.  Default `NULL`.
#' @param branch_r numeric or `NULL`.  If non-NULL, a **second**
#'   [eigen_metrics()] call is made at this fixed search radius (m) — e.g.
#'   `0.1` — and branch candidates are identified by high Linearity at that
#'   scale followed by [lidR::connected_components()] with radius `branch_r`.
#'   This is more physically grounded than the default Curvature/Sphericity
#'   thresholds: structural branches appear locally linear at 0.1 m while
#'   foliage does not.  If `NULL` (default) the Curvature + Sphericity
#'   threshold approach is used instead.
#' @param branch_lin_threshold numeric \[0, 1\].  Minimum Linearity at
#'   `branch_r` for a non-stem, non-ground point to be a branch candidate.
#'   Only used when `branch_r` is non-NULL.  Default `0.35`.
#' @param branch_cc_min_size integer.  Minimum connected-component size (number
#'   of points within `branch_r` of at least one neighbour) to accept a
#'   candidate cluster as a branch.  Smaller clusters are treated as noise and
#'   left as Other.  Only used when `branch_r` is non-NULL.  Default `10L`.
#' @param rescue_missing_stems logical. If `TRUE`, trees with too few stem
#'   points after the seed-centred pass are rechecked using their assigned
#'   `treeID` points. This helps MLS cases where the detected seed is offset
#'   from a visibly vertical stem. Default `FALSE`.
#' @param min_stem_points integer. Minimum stem points per tree before the
#'   rescue pass is skipped. Default `20L`.
#' @param rescue_search_radius numeric or `NULL`. Radius around the rescued
#'   local stem axis. `NULL` uses `search_radius`. Default `NULL`.
#' @param coarse_r numeric or `NULL`.  When non-NULL and `branch_r` is also
#'   non-NULL, activates **two-scale mode**.  [eigen_metrics()] is called at
#'   `branch_r` (fine scale) and `coarse_r` (coarse scale).  Stem candidates
#'   must be vertical at **both** scales (`Vert_fine >=
#'   verticality_threshold` **and** `Vert_coarse >= verticality_threshold`).
#'   Branch candidates must be linear at both scales (`Lin_fine >=
#'   branch_lin_threshold` **and** `Lin_coarse >= branch_lin_threshold`), with
#'   `BranchProb = Lin_fine * Lin_coarse * (1 - Vert_coarse)` peaking for
#'   horizontal structural branches.  Recommended: `coarse_r ~ 0.25` when
#'   `branch_r ~ 0.05`.  Default `NULL`.
#'
#' @return A LAS object identical to `las` with these extra columns:
#' \describe{
#'   \item{`Stem`}{logical — main stem points.}
#'   \item{`StemProb`}{numeric \[0, 1\] — stem confidence.  For the `"lewos"`
#'     method this is the fraction of k_propagation neighbours that voted stem
#'     in the final propagation round; for `"eigen"` it is 0 or 1.}
#'   \item{`Branch`}{logical — structural branch points (not stem, not ground).
#'     Two-scale: `Lin_fine >= branch_lin_threshold` **and** `Lin_coarse >=
#'     branch_lin_threshold`, filtered by connected components.  Single-scale:
#'     `Lin_fine >= branch_lin_threshold`, filtered by connected components.
#'     Default: Curvature < `curvature_threshold` **and** Sphericity <
#'     `sphericity_threshold`.}
#'   \item{`BranchProb`}{numeric \[0, 1\] — branch affinity score.  Two-scale:
#'     `Lin_fine * Lin_coarse * (1 - Vert_coarse)`.  Single-scale: `Lin_fine`.
#'     Default: mean of normalised Curvature and Sphericity scores.}
#'   \item{`Other`}{logical — points that are neither `Stem` nor `Branch`.}
#' }
#'
#' @seealso [classify_stem_points()], [eigen_metrics()], [segment_stem()]
#'
#' @references
#' Wang, D., Momo Takoudjou, S., & Casella, E. (2020). LeWoS: A universal
#' leaf-wood classification method to facilitate the 3D modelling of large
#' tropical trees using terrestrial LiDAR. *Methods in Ecology and Evolution*,
#' 11(3), 376–389. \doi{10.1111/2041-210X.13342}
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
#' las_lw <- classify_lw_points(las_seg, tree_locs)
#' plot(las_lw, color = "StemProb")
#' plot(las_lw, color = "BranchProb")
#' }
#'
#' @export
classify_lw_points <- function(las,
                               tree_locs,
                               method                = c("lewos", "eigen"),
                               z_min                 = 0.1,
                               z_max                 = 50.0,
                               search_radius         = NULL,
                               axis_refine_stems     = TRUE,
                               axis_search_radius    = NULL,
                               axis_min_points       = 10L,
                               neigh_k               = 20L,
                               linearity_threshold   = 0.10,
                               verticality_threshold = 0.75,
                               n_propagation         = 3L,
                               k_propagation         = 15L,
                               curvature_threshold   = 0.15,
                               sphericity_threshold  = 0.29,
                               ncpu                  = 4L,
                               voxel_thin            = NULL,
                               branch_r              = NULL,
                               coarse_r              = NULL,
                               branch_lin_threshold  = 0.35,
                               branch_cc_min_size    = 10L,
                               rescue_missing_stems  = FALSE,
                               min_stem_points       = 20L,
                               rescue_search_radius  = NULL) {

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

  # Drop the legacy Stem column from classify_stem_points() if present.
  if ("Stem" %in% names(las@data))
    las@data[, Stem := NULL]

  # ---- Tree seed geometry ---------------------------------------------------
  tl_xy  <- sf::st_coordinates(tree_locs)
  tl_ids <- tree_locs$TreeID

  if (is.null(search_radius)) {
    search_radius <- max(tree_locs$Radius, na.rm = TRUE) * 3.0
    search_radius <- max(search_radius, 0.20)
  }

  n_pts <- nrow(las@data)

  # ---- Two-scale mode flag -------------------------------------------------
  two_scale <- !is.null(branch_r) && branch_r > 0 &&
               !is.null(coarse_r) && coarse_r > 0

  # ---- Voxel-thinning (optional) -------------------------------------------
  las_thin <- NULL
  nn1      <- NULL

  if (!is.null(voxel_thin) && voxel_thin > 0) {
    message("  [classify_lw_points] Voxel-thinning (res = ", voxel_thin,
            " m) before eigen computation ...")
    las_thin <- lidR::decimate_points(las,
                                      lidR::random_per_voxel(res = voxel_thin,
                                                             n   = 1L))
    nn1 <- FNN::knnx.index(
      data  = cbind(las_thin@data$X, las_thin@data$Y, las_thin@data$Z),
      query = cbind(las@data$X,      las@data$Y,      las@data$Z),
      k     = 1L
    )[, 1L]
  }

  # Compute metrics on the thinned cloud (or full cloud if no thinning).
  # .proj() projects a thinned-cloud metric vector back to all original points.
  target <- if (!is.null(las_thin)) las_thin else las
  .proj  <- function(v) {
    if (!is.null(nn1)) { r <- v[nn1]; r[is.na(r)] <- 0; r }
    else               { v[is.na(v)] <- 0; v }
  }

  # ---- Eigen metrics -------------------------------------------------------
  lin_fine <- vert_fine <- lin_coarse <- vert_coarse <- NULL

  if (two_scale) {
    # Two radius-based passes: fine (branch_r) and coarse (coarse_r).
    # Stem = consistently vertical at both scales.
    # Branch = consistently linear at both scales, not vertical at coarse scale.
    message("  [classify_lw_points] Fine-scale eigen (r = ", branch_r, " m) ...")
    em_fine   <- eigen_metrics(target, r = branch_r,  ncpu = ncpu)
    message("  [classify_lw_points] Coarse-scale eigen (r = ", coarse_r, " m) ...")
    em_coarse <- eigen_metrics(target, r = coarse_r, ncpu = ncpu)
    lin_fine    <- .proj(em_fine$Linearity)
    vert_fine   <- .proj(em_fine$Verticality)
    lin_coarse  <- .proj(em_coarse$Linearity)
    vert_coarse <- .proj(em_coarse$Verticality)

  } else {
    # Single k-NN pass; store metrics in las@data for stem + branch use.
    need_eigen <- !all(c("Linearity", "Verticality",
                          "Curvature", "Sphericity") %in% names(las@data))
    if (need_eigen) {
      message("  [classify_lw_points] Computing eigen_metrics (k = ",
              neigh_k, ") ...")
      em <- eigen_metrics(target, k = as.integer(neigh_k), ncpu = ncpu)
      for (col in c("Linearity", "Verticality", "Curvature", "Sphericity"))
        las@data[[col]] <- .proj(em[[col]])
    }
  }

  # ---- Stem initial threshold -----------------------------------------------
  X     <- las@data$X
  Y     <- las@data$Y
  Z     <- las@data$Z
  tid_v <- las@data$treeID

  # In two-scale mode use scale-specific vectors; otherwise read from las@data.
  lin  <- if (two_scale) lin_fine    else las@data$Linearity
  vert <- if (two_scale) vert_coarse else las@data$Verticality

  ti_match <- match(tid_v, tl_ids)
  tx_v     <- tl_xy[ti_match, 1L]
  ty_v     <- tl_xy[ti_match, 2L]
  sr2      <- search_radius^2

  if (two_scale) {
    # Stem: high verticality at the COARSE scale only.
    # Fine-scale (branch_r) eigen is unstable on sparse voxelized surfaces —
    # a 5 cm sphere on a 5 cm voxel grid often has < 3 neighbours, producing
    # degenerate eigenvalues and Verticality = 0 for genuine stem points.
    # Spatial discrimination (branches vs stems) is handled by search_radius.
    # The fine-scale pass is still used exclusively for branch linearity below.
    stem_flag <- !is.na(ti_match) &
                 Z >= z_min & Z <= z_max &
                 (X - tx_v)^2 + (Y - ty_v)^2 <= sr2 &
                 vert_coarse >= verticality_threshold
  } else {
    stem_flag <- !is.na(ti_match) &
                 !is.na(lin) & !is.na(vert) &
                 Z >= z_min & Z <= z_max &
                 (X - tx_v)^2 + (Y - ty_v)^2 <= sr2 &
                 lin >= linearity_threshold &
                 vert >= verticality_threshold
  }

  metric_stem <- !is.na(ti_match) &
                 !is.na(vert) &
                 Z >= z_min & Z <= z_max &
                 vert >= verticality_threshold
  if (!two_scale) {
    metric_stem <- metric_stem & !is.na(lin) &
      lin >= linearity_threshold
  }

  if (isTRUE(axis_refine_stems)) {
    axis_radius <- if (is.null(axis_search_radius)) {
      search_radius
    } else {
      axis_search_radius
    }
    axis_radius <- max(as.numeric(axis_radius), 0.2)
    axis_min <- as.integer(max(1L, axis_min_points))
    axis_z_hi <- min(z_max, max(z_min + 3.0, 6.0))
    axis_stem <- rep(FALSE, n_pts)

    for (ti in seq_along(tl_ids)) {
      tid <- tl_ids[ti]
      cand_idx <- which(metric_stem & tid_v == tid)
      if (length(cand_idx) == 0L) next

      seed_fallback <- cand_idx[
        (X[cand_idx] - tl_xy[ti, 1])^2 + (Y[cand_idx] - tl_xy[ti, 2])^2 <= sr2
      ]

      if (length(cand_idx) < axis_min) {
        axis_stem[seed_fallback] <- TRUE
        next
      }

      axis_idx <- cand_idx[Z[cand_idx] <= axis_z_hi]
      if (length(axis_idx) < axis_min) axis_idx <- cand_idx

      axis_xy <- cbind(X[axis_idx], Y[axis_idx])
      axis_clusters <- tryCatch(
        dbscan::dbscan(axis_xy,
                       eps = max(0.15, min(axis_radius, search_radius) * 0.5),
                       minPts = max(3L, min(axis_min, 10L))),
        error = function(e) NULL
      )

      if (!is.null(axis_clusters) && any(axis_clusters$cluster > 0L)) {
        cl_ids <- sort(unique(axis_clusters$cluster[axis_clusters$cluster > 0L]))
        centers <- do.call(rbind, lapply(cl_ids, function(cl) {
          cl_idx <- axis_idx[axis_clusters$cluster == cl]
          c(stats::median(X[cl_idx], na.rm = TRUE),
            stats::median(Y[cl_idx], na.rm = TRUE))
        }))
      } else {
        centers <- matrix(
          c(stats::median(X[axis_idx], na.rm = TRUE),
            stats::median(Y[axis_idx], na.rm = TRUE)),
          ncol = 2L
        )
      }

      axis_keep_mask <- rep(FALSE, length(cand_idx))
      for (ci in seq_len(nrow(centers))) {
        axis_keep_mask <- axis_keep_mask |
          (X[cand_idx] - centers[ci, 1])^2 +
          (Y[cand_idx] - centers[ci, 2])^2 <= axis_radius^2
      }
      axis_keep <- cand_idx[axis_keep_mask]

      if (length(axis_keep) >= axis_min) {
        axis_stem[axis_keep] <- TRUE
      } else {
        axis_stem[seed_fallback] <- TRUE
      }
    }

    stem_flag <- axis_stem
  }

  # ---- LeWoS propagation + StemProb ----------------------------------------
  stem_prob <- rep(0.0, n_pts)

  if (method == "lewos" && any(stem_flag)) {
    # Save the strict two-scale seeds so propagation can only ADD stem labels,
    # never demote a point that passed the hard eigen threshold.
    initial_stem_flag <- stem_flag

    prop_idx <- which(
      Z >= z_min & Z <= z_max & !is.na(lin) & !is.na(vert) & tid_v > 0
    )
    n_sub <- length(prop_idx)

    if (n_sub > as.integer(k_propagation) + 1L) {
      coords <- cbind(X[prop_idx], Y[prop_idx], Z[prop_idx])
      nn_idx <- FNN::knn.index(coords, k = as.integer(k_propagation))

      # C++ kernel: iterates majority-vote in-place, returns converged labels
      # + final-round vote fractions (StemProb) in one call.  Eliminates the
      # n_propagation × (n_sub × k_prop) matrix allocations the pure-R loop
      # would make.
      res <- C_knn_majority_vote(nn_idx, stem_flag[prop_idx],
                                 as.integer(n_propagation), as.integer(ncpu))
      stem_flag[prop_idx] <- res$labels
      stem_prob[prop_idx] <- res$vote_frac
    }

    # Restore any seed that was demoted by the majority vote — a point that
    # cleared the strict two-scale threshold cannot be voted out by foliage
    # neighbours.  Force its probability to 1.
    demoted <- initial_stem_flag & !stem_flag
    stem_flag[demoted] <- TRUE
    stem_prob[demoted] <- 1.0

    # Seeds that survived propagation also get probability ≥ their vote
    # fraction, but we floor them at 1.0 so the viewer shows them as fully
    # confirmed stem.
    stem_prob[initial_stem_flag & stem_flag] <- 1.0

  } else {
    # "eigen" method: probability is 0 or 1
    stem_prob <- as.numeric(stem_flag)
  }

  # ---- Rescue trees whose seed-centred stem mask missed the visible stem ----
  # Hough/raster tree seeds can be offset from the actual stem in MLS. The
  # strict first pass uses distance to the seed, so a viable vertical stem can
  # be excluded before propagation. For trees with too few stem points, estimate
  # a local axis from assigned vertical/linear points and relabel candidates
  # close to that local axis.
  if (isTRUE(rescue_missing_stems)) {
    rescue_radius <- if (is.null(rescue_search_radius)) {
      search_radius
    } else {
      rescue_search_radius
    }
    rescue_radius <- max(as.numeric(rescue_radius), 0.2)
    rescue_z_hi <- min(z_max, max(z_min + 3.0, 6.0))
    min_rescue_pts <- as.integer(max(1L, min_stem_points))

    for (tid in tl_ids) {
      tree_idx <- which(tid_v == tid & Z >= z_min & Z <= z_max &
                          !is.na(lin) & !is.na(vert))
      if (length(tree_idx) == 0L) next

      n_existing <- sum(stem_flag[tree_idx], na.rm = TRUE)
      if (n_existing >= min_rescue_pts) next

      cand_idx <- tree_idx[
        lin[tree_idx] >= linearity_threshold &
          vert[tree_idx] >= verticality_threshold
      ]
      if (length(cand_idx) < min_rescue_pts) next

      axis_idx <- cand_idx[Z[cand_idx] <= rescue_z_hi]
      if (length(axis_idx) < min_rescue_pts) axis_idx <- cand_idx

      ax <- stats::median(X[axis_idx], na.rm = TRUE)
      ay <- stats::median(Y[axis_idx], na.rm = TRUE)
      near_axis <- (X[cand_idx] - ax)^2 + (Y[cand_idx] - ay)^2 <=
        rescue_radius^2
      rescued_idx <- cand_idx[near_axis]
      if (length(rescued_idx) < min_rescue_pts) next

      stem_flag[rescued_idx] <- TRUE
      stem_prob[rescued_idx] <- pmax(stem_prob[rescued_idx], 0.95)
    }
  }

  # ---- Branch classification + BranchProb ----------------------------------
  is_ground <- !is.na(las@data$Classification) &
               las@data$Classification == 2L

  if (two_scale) {
    # ---- Two-scale branch detection ----------------------------------------
    # Lin_fine x Lin_coarse: high only when the structure is linear at both
    # scales -> structural branches.  Weighted by (1 - Vert_coarse) to
    # suppress upright stem-like features that narrowly passed the stem filter.
    branch_prob <- rep(0.0, n_pts)
    non_stem_ng <- !stem_flag & !is_ground
    branch_prob[non_stem_ng] <-
      (lin_fine * lin_coarse * (1 - vert_coarse))[non_stem_ng]

    branch_cand <- lin_fine  >= branch_lin_threshold &
                   lin_coarse >= branch_lin_threshold &
                   !stem_flag & !is_ground
    branch_flag <- rep(FALSE, n_pts)
    if (any(branch_cand)) {
      cand_idx <- which(branch_cand)
      las_cand <- lidR::filter_poi(las, branch_cand)
      las_cand <- lidR::connected_components(las_cand,
                                             res     = branch_r,
                                             min_pts = as.integer(branch_cc_min_size),
                                             name    = "CompID")
      branch_flag[cand_idx] <- las_cand@data$CompID > 0L
    }

  } else if (!is.null(branch_r) && branch_r > 0) {
    # ---- Single-scale branch detection -------------------------------------
    message("  [classify_lw_points] Branch-scale eigen (r = ", branch_r, " m) ...")
    em_br  <- eigen_metrics(target, r = branch_r, ncpu = ncpu)
    lin_br <- .proj(em_br$Linearity)

    branch_prob <- rep(0.0, n_pts)
    non_stem_ng <- !stem_flag & !is_ground
    branch_prob[non_stem_ng] <- lin_br[non_stem_ng]

    branch_cand <- lin_br >= branch_lin_threshold & !stem_flag & !is_ground
    branch_flag <- rep(FALSE, n_pts)
    if (any(branch_cand)) {
      cand_idx <- which(branch_cand)
      las_cand <- lidR::filter_poi(las, branch_cand)
      las_cand <- lidR::connected_components(las_cand,
                                             res     = branch_r,
                                             min_pts = as.integer(branch_cc_min_size),
                                             name    = "CompID")
      branch_flag[cand_idx] <- las_cand@data$CompID > 0L
    }

  } else {
    # ---- Default: Curvature + Sphericity thresholds ------------------------
    crv <- las@data$Curvature
    sph <- las@data$Sphericity

    branch_flag <- !is.na(crv) & !is.na(sph) &
                   crv < curvature_threshold &
                   sph < sphericity_threshold &
                   !stem_flag & !is_ground

    crv_score   <- pmax(0.0, 1.0 - crv / curvature_threshold)
    sph_score   <- pmax(0.0, 1.0 - sph / sphericity_threshold)
    branch_prob <- (crv_score + sph_score) / 2.0
    branch_prob[stem_flag | is_ground | is.na(branch_prob)] <- 0.0
  }

  # ---- Other ---------------------------------------------------------------
  other_flag <- !stem_flag & !branch_flag

  # ---- Attach columns to LAS -----------------------------------------------
  las@data[["Stem"]]      <- stem_flag
  las@data[["StemProb"]]  <- stem_prob
  las@data[["Branch"]]    <- branch_flag
  las@data[["BranchProb"]] <- branch_prob
  las@data[["Other"]]     <- other_flag

  return(las)
}
