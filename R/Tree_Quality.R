# ============================================================================
# spanner – Tree_Quality.R
#
# Assess the quality of per-tree bole models and assign quality classes.
#
# Public function: assess_tree_quality()
# ============================================================================


# ----------------------------------------------------------------------------
#' Assess bole-model quality and flag trees for iterative refit
#'
#' `assess_tree_quality` joins the per-segment fit statistics from
#' [segment_bole()] with per-tree volume statistics from
#' [compute_bole_volume()] and allometric information from `tree_locs` to
#' produce a composite fitness score and quality class for each tree.
#'
#' @details
#' Each tree receives a `fit_score` in \[0, 1\] assembled from five equally
#' weighted components (each component is clamped to \[0, 1\]):
#'
#' \describe{
#'   \item{*Residual quality* (30 %)}{`1 − clamp(mean_RANSAC_err / max_mean_error)`}
#'   \item{*Inlier fraction* (20 %)}{Mean per-segment `Inlier_frac` across valid
#'     segments.}
#'   \item{*Bole coverage* (20 %)}{Fraction of the expected bole height that
#'     contains at least one valid segment.}
#'   \item{*Radius consistency* (15 %)}{`1 − clamp(sd(Radius) / max_radius_cv_abs)`,
#'     where `max_radius_cv_abs = mean(Radius) × max_radius_cv`.}
#'   \item{*Allometric plausibility* (15 %)}{1 if the DBH/height ratio is
#'     within `dbh_ht_ratio_range`, linearly decaying outside.}
#' }
#'
#' `quality_class`:
#' \describe{
#'   \item{`"good"`}{`fit_score >= 0.70` **AND** all hard thresholds satisfied.}
#'   \item{`"marginal"`}{`fit_score` in \[0.40, 0.70) or one hard threshold
#'     violated.}
#'   \item{`"bad"`}{`fit_score < 0.40` or multiple hard thresholds violated or
#'     fewer than `min_segments` valid slices.}
#' }
#'
#' @param tree_locs sf object from [get_raster_eigen_treelocs()]. Must have a
#'   `TreeID` column. Optionally has `Height` (m) and `Radius` (m).
#' @param seg_table `data.frame` from [segment_bole()].
#' @param vol_table `data.frame` from [compute_bole_volume()]. Pass `NULL` to
#'   recompute internally (uses cylinder method).
#' @param max_mean_error numeric. Per-tree mean RANSAC residual threshold for
#'   a "good" fit (m). Default `0.02`.
#' @param max_max_error numeric. Per-tree max RANSAC residual hard cutoff (m).
#'   Default `0.05`.
#' @param min_inlier_frac numeric. Minimum mean inlier fraction across valid
#'   segments. Default `0.70`.
#' @param min_bole_coverage numeric. Minimum fraction of expected bole height
#'   with valid segments. Default `0.50`.
#' @param dbh_ht_ratio_range numeric vector of length 2. Plausible
#'   \[DBH(m) / TreeHeight(m)\] range. Default `c(0.004, 0.05)`.
#' @param max_radius_cv numeric. Maximum coefficient of variation for per-slice
#'   radii (sd/mean). Default `0.40`.
#' @param min_segments integer. Minimum valid segments for a non-"bad"
#'   classification. Default `3L`.
#'
#' @return An sf object identical to `tree_locs` with the following extra
#'   columns:
#'   \describe{
#'     \item{`N_segments`}{integer. Valid segment count.}
#'     \item{`Mean_radius`}{numeric. Mean fitted radius (m).}
#'     \item{`SD_radius`}{numeric. SD of fitted radii (m).}
#'     \item{`DBH_m`}{numeric. `2 × Mean_radius` at the lowest valid segment.}
#'     \item{`Mean_RANSAC_err`}{numeric. Mean per-segment RANSAC residual.}
#'     \item{`Max_RANSAC_err`}{numeric. Maximum per-segment RANSAC residual.}
#'     \item{`Mean_inlier_frac`}{numeric. Mean inlier fraction.}
#'     \item{`Bole_coverage`}{numeric. Fraction of expected bole height modeled.}
#'     \item{`DBH_ht_ratio`}{numeric. DBH (m) / tree height (m).}
#'     \item{`Radius_CV`}{numeric. Coefficient of variation of fitted radii.}
#'     \item{`fit_score`}{numeric \[0, 1\]. Composite quality score.}
#'     \item{`quality_class`}{factor. `"good"`, `"marginal"`, or `"bad"`.}
#'   }
#'
#' @seealso [segment_bole()], [compute_bole_volume()], [refit_trees()]
#'
#' @examples
#' \donttest{
#' # (continuing from segment_bole / compute_bole_volume examples)
#' # qual_tbl <- assess_tree_quality(tree_locs, seg_tbl, vol_tbl)
#' # table(qual_tbl$quality_class)
#' }
#'
#' @export
assess_tree_quality <- function(tree_locs,
                                seg_table,
                                vol_table         = NULL,
                                max_mean_error    = 0.02,
                                max_max_error     = 0.05,
                                min_inlier_frac   = 0.70,
                                min_bole_coverage = 0.50,
                                dbh_ht_ratio_range = c(0.004, 0.05),
                                max_radius_cv     = 0.40,
                                min_segments      = 3L) {

  # ---- Validation -----------------------------------------------------------
  if (!inherits(tree_locs, "sf"))
    stop("'tree_locs' must be an sf object from get_raster_eigen_treelocs().")
  if (!"TreeID" %in% names(tree_locs))
    stop("'tree_locs' must have a 'TreeID' column.")
  if (!is.data.frame(seg_table))
    stop("'seg_table' must be a data.frame from segment_bole().")

  required_seg <- c("TreeID", "Segment", "Z_low", "Z_high",
                    "Radius", "Valid", "Inlier_frac", "RANSAC_err")
  miss <- setdiff(required_seg, names(seg_table))
  if (length(miss) > 0)
    stop("'seg_table' is missing columns: ", paste(miss, collapse = ", "))

  # ---- Compute vol_table if not supplied ------------------------------------
  if (is.null(vol_table))
    vol_table <- compute_bole_volume(seg_table, method = "cylinder")

  # ---- Determine expected bole height for each tree -------------------------
  # Prefer a Height column in tree_locs; fall back to Z_high max from seg_table
  has_ht <- "Height" %in% names(tree_locs)
  tl_ht  <- if (has_ht) {
    setNames(tree_locs$Height, as.character(tree_locs$TreeID))
  } else {
    tapply(seg_table$Z_high, seg_table$TreeID, max, na.rm = TRUE)
  }

  tree_ids <- sort(unique(tree_locs$TreeID))
  out_list <- vector("list", length(tree_ids))

  # Pre-split to avoid O(n_trees × n_rows) repeated subsetting in the loop
  valid_mask        <- !is.na(seg_table$Valid) & seg_table$Valid
  segs_all_split    <- split(seg_table, seg_table$TreeID)
  segs_valid_split  <- split(seg_table[valid_mask, , drop = FALSE],
                             seg_table$TreeID[valid_mask])
  vol_split         <- split(vol_table, vol_table$TreeID)

  for (i in seq_along(tree_ids)) {
    tid <- tree_ids[i]
    tid_ch <- as.character(tid)

    segs_all   <- segs_all_split[[tid_ch]]
    if (is.null(segs_all)) segs_all <- seg_table[0L, , drop = FALSE]
    segs_valid <- segs_valid_split[[tid_ch]]
    if (is.null(segs_valid)) segs_valid <- seg_table[0L, , drop = FALSE]
    vol_row    <- vol_split[[tid_ch]]
    if (is.null(vol_row)) vol_row <- vol_table[0L, , drop = FALSE]

    n_seg <- nrow(segs_valid)

    # Expected bole height
    exp_ht <- as.numeric(tl_ht[tid_ch])
    if (is.na(exp_ht) || exp_ht <= 0) exp_ht <- NA_real_

    # Basic metrics
    rads    <- segs_valid$Radius
    errs    <- segs_valid$RANSAC_err
    infrac  <- segs_valid$Inlier_frac

    mean_r   <- mean(rads, na.rm = TRUE)
    sd_r     <- if (n_seg > 1L) stats::sd(rads, na.rm = TRUE) else 0.0
    mean_err <- mean(errs, na.rm = TRUE)
    max_err  <- max(errs,  na.rm = TRUE)
    mean_inf <- mean(infrac, na.rm = TRUE)

    # DBH = 2 × radius of the lowest valid segment
    dbh_m <- if (n_seg >= 1L) {
      lowest <- segs_valid[which.min(segs_valid$Z_mid), , drop = FALSE]
      2.0 * lowest$Radius
    } else NA_real_

    # Bole coverage fraction
    bole_cov <- if (!is.na(exp_ht) && exp_ht > 0 && n_seg >= 1L) {
      ht_span <- max(segs_valid$Z_high, na.rm = TRUE) -
        min(segs_valid$Z_low, na.rm = TRUE)
      min(1.0, ht_span / exp_ht)
    } else NA_real_

    # DBH/height allometric ratio
    dbh_ht <- if (!is.na(dbh_m) && !is.na(exp_ht) && exp_ht > 0)
      dbh_m / exp_ht else NA_real_

    # Radius coefficient of variation
    r_cv <- if (!is.na(mean_r) && mean_r > 0 && !is.na(sd_r))
      sd_r / mean_r else NA_real_

    # ---- Component scores [0, 1] -------------------------------------------
    clamp <- function(x) min(1.0, max(0.0, x))

    sc_resid  <- if (!is.na(mean_err))
      clamp(1.0 - mean_err / max(max_mean_error, 1e-9)) else 0.0
    sc_inlier <- if (!is.na(mean_inf)) clamp(mean_inf)             else 0.0
    sc_cov    <- if (!is.na(bole_cov)) clamp(bole_cov)             else 0.0
    sc_cv     <- if (!is.na(r_cv) && !is.na(mean_r)) {
      max_abs_cv <- mean_r * max_radius_cv
      clamp(1.0 - sd_r / max(max_abs_cv, 1e-9))
    } else 0.0
    sc_allo   <- if (!is.na(dbh_ht)) {
      lo <- dbh_ht_ratio_range[1L]; hi <- dbh_ht_ratio_range[2L]
      if (dbh_ht >= lo && dbh_ht <= hi) {
        1.0
      } else if (dbh_ht < lo) {
        max(0.0, 1.0 - (lo - dbh_ht) / lo)
      } else {
        max(0.0, 1.0 - (dbh_ht - hi) / hi)
      }
    } else 0.5  # neutral when no height available

    fit_score <- 0.30 * sc_resid +
      0.20 * sc_inlier +
      0.20 * sc_cov    +
      0.15 * sc_cv     +
      0.15 * sc_allo

    # ---- Hard threshold violations ----------------------------------------
    hard_fail_count <- 0L
    if (!is.na(mean_err) && mean_err > max_mean_error)  hard_fail_count <- hard_fail_count + 1L
    if (!is.na(max_err)  && max_err  > max_max_error)   hard_fail_count <- hard_fail_count + 1L
    if (!is.na(mean_inf) && mean_inf < min_inlier_frac) hard_fail_count <- hard_fail_count + 1L
    if (!is.na(bole_cov) && bole_cov < min_bole_coverage) hard_fail_count <- hard_fail_count + 1L
    if (!is.na(dbh_ht) &&
        (dbh_ht < dbh_ht_ratio_range[1L] || dbh_ht > dbh_ht_ratio_range[2L]))
      hard_fail_count <- hard_fail_count + 1L

    # ---- Quality class -----------------------------------------------------
    q_class <- if (n_seg < min_segments || hard_fail_count >= 3L ||
                   (!is.na(fit_score) && fit_score < 0.40)) {
      "bad"
    } else if (fit_score >= 0.70 && hard_fail_count == 0L) {
      "good"
    } else {
      "marginal"
    }

    out_list[[i]] <- data.frame(
      TreeID          = tid,
      N_segments      = n_seg,
      Mean_radius     = mean_r,
      SD_radius       = sd_r,
      DBH_m           = dbh_m,
      Mean_RANSAC_err = mean_err,
      Max_RANSAC_err  = max_err,
      Mean_inlier_frac = mean_inf,
      Bole_coverage   = bole_cov,
      DBH_ht_ratio    = dbh_ht,
      Radius_CV       = r_cv,
      fit_score       = fit_score,
      quality_class   = q_class,
      stringsAsFactors = FALSE
    )
  }

  # ---- Join quality metrics back to tree_locs sf ---------------------------
  qual_df <- do.call(rbind, out_list)
  qual_df$quality_class <- factor(qual_df$quality_class,
                                  levels = c("good", "marginal", "bad"))

  # Drop any pre-existing quality columns to avoid duplicate merges
  drop_cols <- intersect(names(qual_df)[-1L], names(tree_locs))
  if (length(drop_cols) > 0)
    tree_locs <- tree_locs[, !names(tree_locs) %in% drop_cols, drop = FALSE]

  result <- merge(tree_locs, qual_df, by = "TreeID", all.x = TRUE,
                  sort = FALSE)

  return(result)
}
