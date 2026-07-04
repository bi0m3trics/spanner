# ============================================================================
# spanner v2.0 – Stem_Points_Methods.R
#
# Method factories for the stem_points() stage.
# Each factory returns a closure of class "spanner_stm_mtd" that accepts
# (las, map) and returns a LAS with per-point stem classification columns.
# ============================================================================


# ----------------------------------------------------------------------------
#' Stem-point method: eigenvalue filtering (+ optional LeWoS propagation)
#'
#' Wraps [classify_lw_points()] as a composable method object.  Uses local
#' PCA eigenvalue metrics — primarily Linearity and Verticality — computed
#' over a k-NN neighbourhood to distinguish stem points from branches and
#' other returns.
#'
#' Internally this calls `classify_lw_points(method = "lewos")` by default,
#' which adds a majority-vote label-propagation pass after the initial
#' threshold filter (Wang et al. 2020 LeWoS framework).  Set `method =
#' "eigen"` for a pure threshold-only classification.
#'
#' @param ... Named arguments forwarded verbatim to [classify_lw_points()].
#'   Key parameters: `method` (`"lewos"` or `"eigen"`), `neigh_k`,
#'   `linearity_threshold`, `verticality_threshold`, `z_min`, `z_max`,
#'   `n_propagation`, `k_propagation`, `curvature_threshold`,
#'   `sphericity_threshold`, `branch_r`, `coarse_r`, `ncpu`, `voxel_thin`.
#'
#' @return A closure of class `"spanner_stm_mtd"` / `"function"` with
#'   attribute `method_name = "eigen"`.  When called it returns the input LAS
#'   with columns `Stem`, `StemProb`, `Branch`, `BranchProb`, `Other` added.
#'
#' @seealso [stem_points()], [stem_hough()], [classify_lw_points()]
#'
#' @examples
#' \donttest{
#' # (continuing from tree_points example)
#' las_s <- stem_points(las_v, map, method = stem_eigen(neigh_k = 20L))
#' }
#'
#' @export
stem_eigen <- function(...) {
  args <- list(...)
  func <- function(las, map) {
    do.call(classify_lw_points, c(list(las = las, tree_locs = map), args))
  }
  structure(func,
            class       = c("spanner_stm_mtd", "function"),
            method_name = "eigen")
}


# ----------------------------------------------------------------------------
#' Stem-point method: Hough slice-stack detector
#'
#' Uses `C_hough_stem_plot()` to detect stem points by tracking the best
#' circular cross-section through stacked horizontal slices, per tree.  A
#' point is labelled as a stem point if it lies within the tracked circle's
#' radius (plus `pixel_size` tolerance) at its corresponding height band.
#'
#' This method is complementary to `stm_eigen()`:
#' \itemize{
#'   \item More robust to mixed-return and high-noise conditions where
#'     eigenvalue thresholds become unreliable.
#'   \item Directly comparable to `stm.hough()` from **TreeLS**.
#'   \item Returns `Stem`, `Segment`, `Radius`, and `Votes` columns rather
#'     than the Stem/Branch/Other columns of `stm_eigen()`.
#' }
#'
#' @param h_base_lo  numeric.  Lower bound of the stem seeding zone (m above
#'   ground).  Default `1.0`.
#' @param h_base_hi  numeric.  Upper bound of the seeding zone (m).  Must be
#'   greater than `h_base_lo`.  Default `2.5`.
#' @param h_step     numeric.  Slice thickness used for tracking (m).
#'   Default `0.5`.
#' @param max_d      numeric.  Maximum stem diameter (m) — sets the Hough
#'   search radius as `max_d / 2`.  Default `0.5`.
#' @param pixel_size numeric.  Hough accumulator resolution (m).  Default
#'   `0.025`.
#' @param min_density numeric \[0, 1\].  Fraction of slice max-count required
#'   for a pixel to vote.  Default `0.1`.
#' @param min_votes  integer.  Minimum Hough votes to accept a circle centre.
#'   Default `3L`.
#'
#' @return A closure of class `"spanner_stm_mtd"` / `"function"` with
#'   attribute `method_name = "hough"`.  When called it returns the input LAS
#'   with columns `Stem` (logical), `Segment` (integer), `Radius` (numeric),
#'   `Votes` (integer) added.
#'
#' @seealso [stem_points()], [stem_eigen()], [classify_lw_points()]
#'
#' @examples
#' \donttest{
#' las_h <- stem_points(las_v, map,
#'            method = stem_hough(h_base_lo = 1, h_base_hi = 2.5,
#'                               max_d = 0.5, min_votes = 5L))
#' }
#'
#' @export
stem_hough <- function(h_base_lo   = 1.0,
                      h_base_hi   = 2.5,
                      h_step      = 0.5,
                      max_d       = 0.5,
                      pixel_size  = 0.025,
                      min_density = 0.1,
                      min_votes   = 3L) {

  stopifnot(
    is.numeric(h_base_lo), is.numeric(h_base_hi), h_base_hi > h_base_lo,
    is.numeric(h_step),     h_step     > 0,
    is.numeric(max_d),      max_d      > 0,
    is.numeric(pixel_size), pixel_size > 0,
    is.numeric(min_density), min_density > 0, min_density <= 1,
    length(min_votes) == 1, min_votes >= 1L
  )

  func <- function(las, map) {
    if (!inherits(las, "LAS"))
      stop("stm_hough: 'las' must be a LAS object", call. = FALSE)
    if (!"treeID" %in% names(las@data))
      stop("stm_hough: 'las' must have a 'treeID' column — run tree_points() first.",
           call. = FALSE)

    # Exclude ground points from stem detection
    has_class <- "Classification" %in% names(las@data)
    survey <- if (has_class) las@data$Classification != 2L
              else           rep(TRUE, nrow(las@data))

    xyz <- as.matrix(las@data[survey, c("X", "Y", "Z")])
    ids <- as.integer(las@data$treeID[survey])

    results <- C_hough_stem_plot(
      las         = xyz,
      tree_ids    = ids,
      h_base_lo   = h_base_lo,
      h_base_hi   = h_base_hi,
      h_step      = h_step,
      max_radius  = max_d / 2,
      pixel_size  = pixel_size,
      min_density = min_density,
      min_votes   = as.integer(min_votes)
    )

    # Initialise output columns
    n <- nrow(las@data)
    las@data$Stem    <- FALSE
    las@data$Segment <- 0L
    las@data$Radius  <- 0.0
    las@data$Votes   <- 0L

    # Write results back to the rows that were surveyed
    las@data$Stem[survey]    <- results$Stem
    las@data$Segment[survey] <- as.integer(results$Segment)
    las@data$Radius[survey]  <- results$Radius
    las@data$Votes[survey]   <- as.integer(results$Votes)

    las
  }

  structure(func,
            class       = c("spanner_stm_mtd", "function"),
            method_name = "hough")
}
