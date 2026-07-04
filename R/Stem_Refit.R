# ============================================================================
# spanner – Stem_Refit.R
#
# Iterative refit of "bad" or "marginal" tree stem models using a cascade of
# increasingly aggressive strategies:
#   1. relax_ransac  – loosen inlier_frac / increase n_samples
#   2. tighten_ransac – tighten search_radius / raise inlier_frac
#   3. dbscan_stem   – DBSCAN per-slice to isolate the main stem cluster
#   4. meanshift_stem – iterative mean-shift per-slice to converge on stem axis
#
# Public function: refit_trees()
# ============================================================================


# ----------------------------------------------------------------------------
#' Iteratively refit stem models for trees that failed quality assessment
#'
#' `refit_trees` identifies trees labelled `"bad"` (or optionally `"marginal"`)
#' in `qual_table` and attempts to produce improved segment tables by running
#' a cascade of refit strategies.  Each strategy is tried in order and the
#' first strategy that lifts a tree out of the `"bad"` class is adopted.  The
#' function returns an updated segment table, an updated quality table, and a
#' log of what was done per tree.
#'
#' @details
#' Strategies (in cascade order):
#' \describe{
#'   \item{`"relax_ransac"`}{Re-runs [segment_stem()] with `inlier_frac`
#'     reduced by `relax_delta` and `n_ransac` doubled.  Useful when the
#'     stem is well-separated but RANSAC is too strict.}
#'   \item{`"tighten_ransac"`}{Re-runs [segment_stem()] with `search_radius`
#'     halved and `inlier_frac` increased by `relax_delta`.  Useful when
#'     co-mingled understorey points pollute the slice.}
#'   \item{`"dbscan_stem"`}{For each height slice of a problem tree, runs
#'     `dbscan::dbscan()` on the XY coordinates and keeps only the cluster
#'     whose centroid is closest to the tree seed location.  The filtered
#'     points are then fed back into [segment_stem()] with the `"ransac"`
#'     algorithm.  Requires the **dbscan** package.}
#'   \item{`"meanshift_stem"`}{For each height slice, performs simple
#'     iterative mean-shift (bandwidth = `ms_bandwidth`) to locate the
#'     densest XY region and restricts fitting to points within that
#'     neighbourhood.  Does **not** require additional packages.}
#' }
#'
#' @param las LAS object (output of [segment_graph()]; optionally has `Stem`
#'   column from [classify_stem_points()]).
#' @param tree_locs sf object from [get_raster_eigen_treelocs()].
#' @param seg_table `data.frame` from [segment_stem()].
#' @param qual_table sf from [assess_tree_quality()].
#' @param strategies character vector. Ordered list of strategies to attempt.
#'   Default `c("relax_ransac", "tighten_ransac", "dbscan_stem",
#'   "meanshift_stem")`.
#' @param refit_marginal logical. Also attempt to refit `"marginal"` trees.
#'   Default `FALSE`.
#' @param tree_ids optional integer vector of tree IDs to consider for refit.
#'   `NULL` uses all trees in the requested quality classes. Default `NULL`.
#' @param max_iterations integer. Maximum number of strategy passes per tree.
#'   Default `3L`.
#' @param z_min numeric. Lower stem height (m). Default `0.1`.
#' @param z_max numeric. Upper stem height (m). `NULL` = per-tree max.
#' @param dz numeric. Slice height (m). Default `0.5`.
#' @param overlap numeric. Slice overlap fraction. Default `0.1`.
#' @param search_radius numeric. Initial XY search radius (m). `NULL` =
#'   auto-detect.
#' @param inlier_tol numeric. RANSAC inlier tolerance (m). Default `0.03`.
#' @param n_ransac integer. Initial RANSAC samples. Default `20L`.
#' @param conf numeric. RANSAC confidence. Default `0.99`.
#' @param inlier_frac numeric. Initial inlier fraction. Default `0.85`.
#' @param n_best integer. RANSAC best-of aggregation. Default `25L`.
#' @param max_ransac_iterations integer. Maximum allowed RANSAC iterations per
#'   slice in refit strategies. If a strategy's `inlier_frac` and sample count
#'   would exceed this budget, the sample count is reduced automatically.
#'   Default `1000L`.
#' @param max_radius numeric. Maximum valid radius (m). Default `1.0`.
#' @param min_pts integer. Minimum slice points. Default `10L`.
#' @param relax_delta numeric. Inlier-fraction change per relaxation step.
#'   Default `0.15`.
#' @param dbscan_eps_factor numeric. DBSCAN `eps = mean_radius * dbscan_eps_factor`.
#'   Default `2.0`.
#' @param dbscan_min_pts integer. DBSCAN `minPts`. Default `5L`.
#' @param ms_bandwidth numeric. Mean-shift kernel bandwidth (m). `NULL` =
#'   `search_radius / 3`.
#' @param ms_max_iter integer. Max mean-shift iterations per slice. Default `30L`.
#' @param ms_tol numeric. Mean-shift convergence tolerance (m). Default `0.001`.
#' @param ncpu integer. OpenMP threads passed to the internal RANSAC C++ routine
#'   (parallelises across height slices within each per-tree call). Default `1L`.
#' @param quality_args list. Additional arguments forwarded to
#'   [assess_tree_quality()] for re-evaluation after refit. Default `list()`.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{`seg_table`}{Updated segment table (all trees, bad trees replaced).}
#'     \item{`qual_table`}{Updated quality table (sf).}
#'     \item{`refit_log`}{`data.frame` with columns `TreeID`, `original_class`,
#'       `strategy_used`, `final_class`, `n_attempts`.}
#'   }
#'
#' @seealso [segment_stem()], [assess_tree_quality()], [classify_stem_points()]
#'
#' @examples
#' \donttest{
#' # (continuing from assess_tree_quality example)
#' # refit_result <- refit_trees(las_seg, tree_locs, seg_tbl, qual_tbl)
#' # final_segs   <- refit_result$seg_table
#' # final_qual   <- refit_result$qual_table
#' # print(refit_result$refit_log)
#' }
#'
#' @export
refit_trees <- function(las,
                        tree_locs,
                        seg_table,
                        qual_table,
                        strategies    = c("relax_ransac", "tighten_ransac",
                                          "dbscan_stem", "meanshift_stem"),
                        refit_marginal = FALSE,
                        tree_ids       = NULL,
                        max_iterations = 3L,
                        z_min          = 0.1,
                        z_max          = NULL,
                        dz             = 0.5,
                        overlap        = 0.1,
                        search_radius  = NULL,
                        inlier_tol     = 0.03,
                        n_ransac       = 20L,
                        conf           = 0.99,
                        inlier_frac    = 0.85,
                        n_best         = 25L,
                        max_ransac_iterations = 1000L,
                        max_radius     = 1.0,
                        min_pts        = 10L,
                        relax_delta    = 0.15,
                        dbscan_eps_factor = 2.0,
                        dbscan_min_pts = 5L,
                        ms_bandwidth   = NULL,
                        ms_max_iter    = 30L,
                        ms_tol         = 0.001,
                        ncpu           = 1L,
                        quality_args   = list()) {

  # ---- Validation -----------------------------------------------------------
  strategies <- match.arg(strategies,
                          choices = c("relax_ransac", "tighten_ransac",
                                      "dbscan_stem", "meanshift_stem"),
                          several.ok = TRUE)

  if ("dbscan_stem" %in% strategies &&
      !requireNamespace("dbscan", quietly = TRUE))
    stop("Strategy 'dbscan_stem' requires the dbscan package. ",
         "Install with: install.packages('dbscan')")

  if (!inherits(qual_table, "sf"))
    stop("'qual_table' must be the sf from assess_tree_quality().")
  if (!"quality_class" %in% names(qual_table))
    stop("'qual_table' must contain a 'quality_class' column.")

  # ---- Identify problem trees -----------------------------------------------
  target_classes <- if (refit_marginal) c("bad", "marginal") else "bad"
  target_ids     <- qual_table$TreeID[
    qual_table$quality_class %in% target_classes
  ]
  if (!is.null(tree_ids)) {
    target_ids <- intersect(target_ids, as.integer(tree_ids))
  }

  if (length(target_ids) == 0L) {
    message("[refit_trees] No trees require refitting.")
    return(list(
      seg_table  = seg_table,
      qual_table = qual_table,
      refit_log  = data.frame(TreeID = integer(), original_class = character(),
                              strategy_used = character(),
                              final_class = character(),
                              n_attempts = integer(),
                              stringsAsFactors = FALSE)
    ))
  }

  message("[refit_trees] Attempting refit for ", length(target_ids),
          " tree(s): ", paste(head(target_ids, 10), collapse = ", "),
          if (length(target_ids) > 10) " ..." else "")

  # ---- Set up geometry ------------------------------------------------------
  tl_xy  <- sf::st_coordinates(tree_locs)
  tl_ids <- as.integer(tree_locs$TreeID)

  if (is.null(search_radius)) {
    search_radius <- max(tree_locs$Radius, na.rm = TRUE) * 3.0
    search_radius <- max(search_radius, 0.20)
  }
  if (is.null(ms_bandwidth)) ms_bandwidth <- search_radius / 3.0

  # Persistent copies that we update iteratively
  cur_seg   <- seg_table
  cur_qual  <- qual_table

  log_rows  <- vector("list", length(target_ids))

  # ---- Per-tree cascade -----------------------------------------------------
  for (ii in seq_along(target_ids)) {
    tid <- target_ids[ii]
    orig_class <- as.character(
      qual_table$quality_class[qual_table$TreeID == tid])

    ti <- which(tl_ids == tid)
    if (length(ti) == 0L) next
    tx <- tl_xy[ti[1L], 1L]
    ty <- tl_xy[ti[1L], 2L]

    # Extract this tree's LAS data (guard against NA treeIDs in unclassified pts)
    tree_mask <- !is.na(las@data$treeID) & las@data$treeID == tid
    if (sum(tree_mask) < min_pts) next
    tlas_data <- las@data[tree_mask, ]

    # Determine per-tree z_max
    zmax_tree <- if (!is.null(z_max)) z_max else
      max(tlas_data$Z, na.rm = TRUE)
    zmax_tree <- max(zmax_tree, z_min + dz)

    strategy_used <- NA_character_
    n_attempts    <- 0L
    final_class   <- orig_class

    for (strategy in strategies) {
      if (n_attempts >= max_iterations) break
      n_attempts <- n_attempts + 1L

      new_segs <- tryCatch(
        .refit_strategy(
          strategy      = strategy,
          tlas_data     = tlas_data,
          tid           = tid,
          tx            = tx, ty = ty,
          z_min         = z_min, z_max = zmax_tree,
          dz            = dz, overlap = overlap,
          search_radius = search_radius,
          inlier_tol    = inlier_tol,
          n_ransac      = n_ransac, conf = conf,
          inlier_frac   = inlier_frac, n_best = n_best,
          max_ransac_iterations = max_ransac_iterations,
          max_radius    = max_radius, min_pts = min_pts,
          relax_delta   = relax_delta,
          dbscan_eps_factor = dbscan_eps_factor,
          dbscan_min_pts = dbscan_min_pts,
          ms_bandwidth  = ms_bandwidth,
          ms_max_iter   = ms_max_iter, ms_tol = ms_tol,
          ncpu          = ncpu
        ),
        error = function(e) {
          message("  [refit_trees] Strategy '", strategy,
                  "' failed for tree ", tid, ": ", conditionMessage(e))
          NULL
        }
      )

      if (is.null(new_segs) || nrow(new_segs) == 0L) next

      # Build a temporary qual table for this tree to assess improvement
      aligned <- .align_refit_segments(
        cur_seg[cur_seg$TreeID != tid, , drop = FALSE],
        new_segs
      )
      tmp_seg <- rbind(
        aligned$existing,
        aligned$replacement
      )

      qa_args <- c(
        list(tree_locs = tree_locs[tree_locs$TreeID == tid, ],
             seg_table = new_segs,
             vol_table = NULL),
        quality_args
      )
      tmp_qual <- tryCatch(do.call(assess_tree_quality, qa_args),
                           error = function(e) NULL)

      if (is.null(tmp_qual)) next

      new_class <- as.character(
        tmp_qual$quality_class[tmp_qual$TreeID == tid])

      if (length(new_class) == 0L) next

      strategy_used <- strategy
      final_class   <- new_class

      # Update persistent tables regardless of improvement (take the best we got)
      cur_seg <- tmp_seg

      # Update quality row for this tree in cur_qual
      qual_cols <- setdiff(names(tmp_qual), attr(tmp_qual, "sf_column"))
      for (qc in qual_cols) {
        if (qc %in% names(cur_qual)) {
          cur_qual[[qc]][cur_qual$TreeID == tid] <- tmp_qual[[qc]][tmp_qual$TreeID == tid]
        }
      }

      # Stop cascade if no longer "bad" (or "marginal" if that was also targeted)
      if (!final_class %in% target_classes) break
    }

    log_rows[[ii]] <- data.frame(
      TreeID        = tid,
      original_class = orig_class,
      strategy_used  = if (is.na(strategy_used)) "none" else strategy_used,
      final_class    = final_class,
      n_attempts     = n_attempts,
      stringsAsFactors = FALSE
    )
  }

  refit_log <- do.call(rbind, log_rows[!sapply(log_rows, is.null)])
  if (is.null(refit_log))
    refit_log <- data.frame(TreeID = integer(), original_class = character(),
                            strategy_used = character(),
                            final_class = character(),
                            n_attempts = integer(),
                            stringsAsFactors = FALSE)
  rownames(refit_log) <- NULL

  return(list(
    seg_table  = cur_seg,
    qual_table = cur_qual,
    refit_log  = refit_log
  ))
}


# ============================================================================
# Internal: dispatch a single refit strategy for one tree
# ============================================================================
.refit_strategy <- function(strategy, tlas_data, tid, tx, ty,
                            z_min, z_max, dz, overlap,
                            search_radius, inlier_tol,
                            n_ransac, conf, inlier_frac, n_best,
                            max_ransac_iterations,
                            max_radius, min_pts,
                            relax_delta, dbscan_eps_factor, dbscan_min_pts,
                            ms_bandwidth, ms_max_iter, ms_tol,
                            ncpu = 1L) {

  X   <- tlas_data$X
  Y   <- tlas_data$Y
  Z   <- tlas_data$Z
  sr2 <- search_radius^2

  # Stem-candidate mask (within search cylinder + z range)
  base_mask <- Z >= z_min & Z <= z_max &
    (X - tx)^2 + (Y - ty)^2 <= sr2
  stem_col <- "Stem"
  if (stem_col %in% names(tlas_data) && any(!is.na(tlas_data[[stem_col]])))
    base_mask <- base_mask & !is.na(tlas_data[[stem_col]]) & tlas_data[[stem_col]]

  step     <- dz * max(0.01, 1.0 - overlap)
  n_slices <- max(1L, ceiling((z_max - z_min) / step))

  # ---- Strategy-specific logic returns a filtered point set per slice -------
  # We'll build slice data and call RANSAC on each

  if (strategy == "relax_ransac") {
    new_if  <- max(0.40, inlier_frac - relax_delta)
    new_ns  <- .bounded_ransac_samples(
      requested = as.integer(n_ransac),
      inlier_frac = new_if,
      confidence = conf,
      n_best = n_best,
      max_iterations = max_ransac_iterations
    )
    return(.slice_and_ransac(X, Y, Z, base_mask, tid, tx, ty,
                             z_min, z_max, dz, overlap, search_radius,
                             inlier_tol, new_ns, conf, new_if, n_best,
                             max_radius, min_pts, ncpu))
  }

  if (strategy == "tighten_ransac") {
    new_sr  <- search_radius * 0.5
    new_if  <- min(0.97, inlier_frac + relax_delta)
    new_ns  <- .bounded_ransac_samples(
      requested = as.integer(n_ransac),
      inlier_frac = new_if,
      confidence = conf,
      n_best = n_best,
      max_iterations = max_ransac_iterations
    )
    tight_mask <- Z >= z_min & Z <= z_max &
      (X - tx)^2 + (Y - ty)^2 <= new_sr^2
    if (stem_col %in% names(tlas_data) && any(!is.na(tlas_data[[stem_col]])))
      tight_mask <- tight_mask & !is.na(tlas_data[[stem_col]]) & tlas_data[[stem_col]]
    return(.slice_and_ransac(X, Y, Z, tight_mask, tid, tx, ty,
                             z_min, z_max, dz, overlap, new_sr,
                             inlier_tol, new_ns, conf, new_if, n_best,
                             max_radius, min_pts, ncpu))
  }

  if (strategy == "dbscan_stem") {
    # Estimate expected radius from the seg_table if available
    # (We don't have it here, use search_radius / 3 as heuristic)
    eps <- (search_radius / 3.0) * dbscan_eps_factor
    return(.slice_dbscan_ransac(X, Y, Z, base_mask, tid, tx, ty,
                                z_min, z_max, dz, overlap, search_radius,
                                inlier_tol,
                                .bounded_ransac_samples(
                                  requested = as.integer(n_ransac),
                                  inlier_frac = inlier_frac,
                                  confidence = conf,
                                  n_best = n_best,
                                  max_iterations = max_ransac_iterations
                                ),
                                conf, inlier_frac,
                                n_best, max_radius, min_pts,
                                eps, as.integer(dbscan_min_pts), ncpu))
  }

  if (strategy == "meanshift_stem") {
    return(.slice_meanshift_ransac(X, Y, Z, base_mask, tid, tx, ty,
                                   z_min, z_max, dz, overlap, search_radius,
                                   inlier_tol,
                                   .bounded_ransac_samples(
                                     requested = as.integer(n_ransac),
                                     inlier_frac = inlier_frac,
                                     confidence = conf,
                                     n_best = n_best,
                                     max_iterations = max_ransac_iterations
                                   ),
                                   conf, inlier_frac,
                                   n_best, max_radius, min_pts,
                                   ms_bandwidth, as.integer(ms_max_iter), ms_tol, ncpu))
  }

  stop("Unknown strategy: ", strategy)
}


# ---------------------------------------------------------------------------
# Helper: align legacy refit rows with stem_segments() / segment_stem() schemas
# ---------------------------------------------------------------------------
.align_refit_segments <- function(existing, replacement) {
  all_cols <- unique(c(names(existing), names(replacement)))

  add_missing <- function(x) {
    missing_cols <- setdiff(all_cols, names(x))
    for (nm in missing_cols) x[[nm]] <- NA
    x[, all_cols, drop = FALSE]
  }

  list(
    existing    = add_missing(existing),
    replacement = add_missing(replacement)
  )
}


# ---------------------------------------------------------------------------
# Helper: keep the C++ RANSAC iteration count within a sane per-slice budget
# ---------------------------------------------------------------------------
.bounded_ransac_samples <- function(requested, inlier_frac, confidence, n_best,
                                    max_iterations) {
  requested <- max(3L, as.integer(requested))
  max_iterations <- max(1L, as.integer(max_iterations))

  iterations_for <- function(n_samples) {
    p_success <- inlier_frac^n_samples
    if (!is.finite(p_success) || p_success <= 0 || p_success >= 1)
      return(Inf)

    k_base <- ceiling(log(1 - confidence) / log(1 - p_success))
    (as.integer(n_best) + 5L) * k_base
  }

  if (iterations_for(requested) <= max_iterations)
    return(requested)

  candidates <- requested:3L
  costs <- vapply(candidates, iterations_for, numeric(1L))
  ok <- candidates[costs <= max_iterations]

  if (length(ok) > 0L) max(ok) else 3L
}


# ---------------------------------------------------------------------------
# Helper: slice the stem + call C_ransac_stem_slices for a single tree
# ---------------------------------------------------------------------------
.slice_and_ransac <- function(X, Y, Z, mask, tid, tx, ty,
                              z_min, z_max, dz, overlap, search_radius,
                              inlier_tol, n_ransac, conf, inlier_frac,
                              n_best, max_radius, min_pts, ncpu = 1L) {

  X_s  <- X[mask]; Y_s <- Y[mask]; Z_s <- Z[mask]
  if (length(X_s) < min_pts) return(NULL)

  result <- C_ransac_stem_slices(
    X            = X_s, Y = Y_s, Z = Z_s,
    treeIds      = rep(tid, length(X_s)),
    tree_X       = tx, tree_Y = ty,
    tree_ID_vals = as.integer(tid),
    z_min        = z_min, z_max = z_max,
    dz           = dz, overlap = overlap,
    search_radius = search_radius,
    inlier_tol   = inlier_tol,
    n_samples    = as.integer(n_ransac),
    confidence   = conf,
    inlier_frac  = inlier_frac,
    n_best       = as.integer(n_best),
    max_radius   = max_radius,
    min_pts      = as.integer(min_pts),
    ncpu         = as.integer(ncpu)
  )
  result$Diameter   <- result$Radius * 2.0
  result$Fit_method <- "ransac_refit"
  result
}


# ---------------------------------------------------------------------------
# Helper: DBSCAN per slice → keep stem cluster → single batch RANSAC call
# ---------------------------------------------------------------------------
.slice_dbscan_ransac <- function(X, Y, Z, base_mask, tid, tx, ty,
                                 z_min, z_max, dz, overlap, search_radius,
                                 inlier_tol, n_ransac, conf, inlier_frac,
                                 n_best, max_radius, min_pts, eps, dbs_min,
                                 ncpu = 1L) {

  X_b <- X[base_mask]; Y_b <- Y[base_mask]; Z_b <- Z[base_mask]
  if (length(X_b) < min_pts) return(NULL)

  step      <- dz * max(0.01, 1.0 - overlap)
  n_slices  <- max(1L, ceiling((z_max - z_min) / step))
  stem_flag <- rep(FALSE, length(X_b))

  # Pre-filter: DBSCAN per slice, flag points in the nearest stem cluster
  for (s in seq_len(n_slices)) {
    zlo   <- z_min + (s - 1L) * step
    zhi   <- zlo + dz
    s_idx <- which(Z_b >= zlo & Z_b <= zhi)
    if (length(s_idx) < min_pts) next

    pts2d  <- cbind(X_b[s_idx], Y_b[s_idx])
    db_res <- tryCatch(
      dbscan::dbscan(pts2d, eps = eps, minPts = dbs_min),
      error = function(e) NULL
    )
    if (is.null(db_res) || all(db_res$cluster == 0L)) next

    cl_ids  <- sort(unique(db_res$cluster))
    cl_ids  <- cl_ids[cl_ids > 0L]
    best_cl <- cl_ids[which.min(vapply(cl_ids, function(cl) {
      m  <- db_res$cluster == cl
      dx <- mean(X_b[s_idx][m]) - tx
      dy <- mean(Y_b[s_idx][m]) - ty
      dx^2 + dy^2
    }, numeric(1L)))]

    stem_flag[s_idx[db_res$cluster == best_cl]] <- TRUE
  }

  keep <- which(stem_flag)
  if (length(keep) < min_pts) return(NULL)

  # Single batch RANSAC call on all pre-filtered stem points
  result <- C_ransac_stem_slices(
    X            = X_b[keep], Y = Y_b[keep], Z = Z_b[keep],
    treeIds      = rep(tid, length(keep)),
    tree_X       = tx, tree_Y = ty,
    tree_ID_vals = as.integer(tid),
    z_min        = z_min, z_max = z_max,
    dz           = dz, overlap = overlap,
    search_radius = max_radius * 2.0,
    inlier_tol   = inlier_tol,
    n_samples    = as.integer(n_ransac),
    confidence   = conf,
    inlier_frac  = inlier_frac,
    n_best       = as.integer(n_best),
    max_radius   = max_radius,
    min_pts      = as.integer(min_pts),
    ncpu         = as.integer(ncpu)
  )

  if (is.null(result) || nrow(result) == 0L) return(NULL)
  result$Diameter   <- result$Radius * 2.0
  result$Fit_method <- "dbscan_ransac_refit"
  rownames(result)  <- NULL
  result
}


# ---------------------------------------------------------------------------
# Helper: mean-shift per slice → single batch RANSAC on converged points
# ---------------------------------------------------------------------------
.slice_meanshift_ransac <- function(X, Y, Z, base_mask, tid, tx, ty,
                                    z_min, z_max, dz, overlap, search_radius,
                                    inlier_tol, n_ransac, conf, inlier_frac,
                                    n_best, max_radius, min_pts,
                                    bw, ms_max_iter, ms_tol, ncpu = 1L) {

  X_b <- X[base_mask]; Y_b <- Y[base_mask]; Z_b <- Z[base_mask]
  if (length(X_b) < min_pts) return(NULL)

  step      <- dz * max(0.01, 1.0 - overlap)
  n_slices  <- max(1L, ceiling((z_max - z_min) / step))
  bw2       <- bw^2
  stem_flag <- rep(FALSE, length(X_b))

  # Pre-filter: mean-shift per slice, flag converged-neighbourhood points
  for (s in seq_len(n_slices)) {
    zlo   <- z_min + (s - 1L) * step
    zhi   <- zlo + dz
    s_idx <- which(Z_b >= zlo & Z_b <= zhi)
    if (length(s_idx) < min_pts) next

    sX <- X_b[s_idx]; sY <- Y_b[s_idx]
    cx <- tx; cy <- ty

    for (iter in seq_len(ms_max_iter)) {
      d2   <- (sX - cx)^2 + (sY - cy)^2
      inbd <- d2 <= bw2
      if (sum(inbd) < 3L) break
      new_cx <- mean(sX[inbd])
      new_cy <- mean(sY[inbd])
      if ((new_cx - cx)^2 + (new_cy - cy)^2 < ms_tol^2) {
        cx <- new_cx; cy <- new_cy; break
      }
      cx <- new_cx; cy <- new_cy
    }

    stem_flag[s_idx[(sX - cx)^2 + (sY - cy)^2 <= bw2]] <- TRUE
  }

  keep <- which(stem_flag)
  if (length(keep) < min_pts) return(NULL)

  # Single batch RANSAC call on all pre-filtered stem points
  result <- C_ransac_stem_slices(
    X            = X_b[keep], Y = Y_b[keep], Z = Z_b[keep],
    treeIds      = rep(tid, length(keep)),
    tree_X       = tx, tree_Y = ty,
    tree_ID_vals = as.integer(tid),
    z_min        = z_min, z_max = z_max,
    dz           = dz, overlap = overlap,
    search_radius = bw * 2.0,
    inlier_tol   = inlier_tol,
    n_samples    = as.integer(n_ransac),
    confidence   = conf,
    inlier_frac  = inlier_frac,
    n_best       = as.integer(n_best),
    max_radius   = max_radius,
    min_pts      = as.integer(min_pts),
    ncpu         = as.integer(ncpu)
  )

  if (is.null(result) || nrow(result) == 0L) return(NULL)
  result$Diameter   <- result$Radius * 2.0
  result$Fit_method <- "meanshift_ransac_refit"
  rownames(result)  <- NULL
  result
}
