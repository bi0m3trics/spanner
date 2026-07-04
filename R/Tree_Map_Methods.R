# ============================================================================
# spanner v2.0 – Tree_Map_Methods.R
#
# Method factories for the find_trees() stage.
# Each factory returns a closure of class "spanner_map_mtd" that accepts
# a single LAS argument and returns an sf POINT object compatible with
# get_raster_eigen_treelocs().
# ============================================================================


# ----------------------------------------------------------------------------
#' Tree-map method: raster-eigen pipeline
#'
#' Wraps [get_raster_eigen_treelocs()] as a composable method object for use
#' inside [find_trees()].  All parameters are forwarded unchanged so existing
#' call-sites can be migrated without any parameter changes.
#'
#' @param ... Named arguments forwarded verbatim to
#'   [get_raster_eigen_treelocs()].  Any parameter that function accepts can
#'   be supplied here (e.g. `res`, `pt_spacing`, `dens_threshold`,
#'   `neigh_sizes`, `eigen_threshold`, `grid_slice_min`, `grid_slice_max`,
#'   `minimum_polygon_area`, `cylinder_fit_type`, `max_dia`, `SDvert`,
#'   `n_pts`, `n_best`, `inliers`, `conf`, `max_angle`).
#'
#' @return A closure of class `"spanner_map_mtd"` / `"function"` with
#'   attribute `method_name = "raster_eigen"`.
#'
#' @seealso [find_trees()], [use_hough()], [get_raster_eigen_treelocs()]
#'
#' @examples
#' \donttest{
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#'
#' map <- find_trees(las, method = use_raster_eigen(
#'   res = 0.025, pt_spacing = 0.0254, dens_threshold = 0.25,
#'   neigh_sizes = c(0.25, 0.15, 0.66), eigen_threshold = 0.75,
#'   grid_slice_min = 1, grid_slice_max = 2,
#'   minimum_polygon_area = 0.005, cylinder_fit_type = "ransac",
#'   max_dia = 1, SDvert = 0.33, n_pts = 20, n_best = 25,
#'   inliers = 0.9, conf = 0.99, max_angle = 20
#' ))
#' }
#'
#' @export
use_raster_eigen <- function(...) {
  args <- list(...)
  func <- function(las) {
    do.call(get_raster_eigen_treelocs, c(list(las = las), args))
  }
  structure(func,
            class       = c("spanner_map_mtd", "function"),
            method_name = "raster_eigen")
}


# ----------------------------------------------------------------------------
#' Tree-map method: Hough slice-stack
#'
#' Detects tree positions by scanning stacked horizontal slices with the
#' circular Hough Transform.  Finds vertically persistent circular
#' cross-sections and returns one representative tree centre per detected stem.
#'
#' The algorithm mirrors `map.hough()` from the **TreeLS** package but runs
#' entirely on spanner's own compiled C++ back end.  Key design choices:
#' \itemize{
#'   \item Height slices of thickness `h_step` are scanned between `min_h`
#'     and `max_h`.
#'   \item The Hough accumulator casts votes at radius increments of
#'     `pixel_size` up to `max_d / 2`.
#'   \item A circle centre is accepted if it accumulates at least `min_votes`
#'     and the local pixel count exceeds `min_density * max_count`.
#'   \item A candidate is promoted to a tree only when it appears in at least
#'     75 % of the height slices (vertical persistence filter).
#' }
#'
#' @param min_h      numeric.  Lower height bound for slice scanning (m).
#'   Default `1.0`.
#' @param max_h      numeric.  Upper height bound (m).  Must be > `min_h`.
#'   Default `3.0`.
#' @param h_step     numeric.  Slice thickness (m).  Default `0.5`.
#' @param pixel_size numeric.  Hough accumulator grid resolution (m).
#'   Smaller values give finer radius precision but increase memory.
#'   Default `0.025`.
#' @param max_d      numeric.  Maximum stem diameter searched (m).  This sets
#'   the maximum circle radius as `max_d / 2`.  Default `0.5`.
#' @param min_density numeric \[0, 1\].  Minimum fraction of the local
#'   maximum point count for a pixel to contribute a Hough vote.
#'   Default `0.1`.
#' @param min_votes  integer.  Minimum accumulated votes to accept a circle
#'   centre.  Lower values find more (possibly false) stems; higher values
#'   are more conservative.  Default `3L`.
#'
#' @return A closure of class `"spanner_map_mtd"` / `"function"` with
#'   attribute `method_name = "hough"`.  When called with a LAS object it
#'   returns an sf POINT object with columns `TreeID`, `Radius`, `Error`
#'   identical in format to [get_raster_eigen_treelocs()] output.
#'
#' @seealso [find_trees()], [use_raster_eigen()]
#'
#' @examples
#' \donttest{
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#'
#' map <- find_trees(las, method = use_hough(
#'   min_h = 1, max_h = 3, h_step = 0.5, max_d = 0.5,
#'   min_density = 0.1, min_votes = 3
#' ))
#' }
#'
#' @export
use_hough <- function(min_h       = 1.0,
                      max_h       = 3.0,
                      h_step      = 0.5,
                      pixel_size  = 0.025,
                      max_d       = 0.5,
                      min_density = 0.1,
                      min_votes   = 3L) {

  stopifnot(
    is.numeric(min_h), is.numeric(max_h), max_h > min_h,
    is.numeric(h_step),     h_step     > 0,
    is.numeric(pixel_size), pixel_size > 0,
    is.numeric(max_d),      max_d      > 0,
    is.numeric(min_density), min_density > 0, min_density <= 1,
    length(min_votes) == 1, min_votes >= 1L
  )

  func <- function(las) {
    if (!inherits(las, "LAS"))
      stop("use_hough: 'las' must be a LAS object")

    xyz <- as.matrix(las@data[, c("X", "Y", "Z")])

    raw <- C_stack_map(
      las         = xyz,
      h_min       = min_h,
      h_max       = max_h,
      h_step      = h_step,
      pixel_size  = pixel_size,
      max_radius  = max_d / 2,
      min_density = min_density,
      min_votes   = as.integer(min_votes)
    )

    if (length(raw) == 0 || length(raw[["X"]]) == 0)
      stop("use_hough: no trees found — try relaxing min_density / min_votes.",
           call. = FALSE)

    map_df <- as.data.frame(raw, stringsAsFactors = FALSE)

    # TreePosition == TRUE rows are the mean tree centre positions
    tree_pos <- map_df[map_df$TreePosition == TRUE, , drop = FALSE]
    if (nrow(tree_pos) == 0)
      stop("use_hough: could not determine tree positions.", call. = FALSE)

    tree_pos$TreeID <- as.integer(tree_pos$TreeID)

    # Mean radius from the keypoint disc candidates for each tree
    key_rows <- map_df[map_df$Keypoint == TRUE & !map_df$TreePosition, ]
    if (nrow(key_rows) > 0) {
      mean_r_by_tree <- tapply(key_rows$Radius, key_rows$TreeID,
                               mean, na.rm = TRUE)
      tree_pos$Radius <- as.numeric(
        mean_r_by_tree[as.character(tree_pos$TreeID)]
      )
      tree_pos$Radius[is.na(tree_pos$Radius)] <- 0.05
    } else {
      tree_pos$Radius <- 0.05
    }

    tree_pos$Error <- NA_real_

    sf_out <- sf::st_as_sf(tree_pos,
                            coords = c("X", "Y"),
                            crs    = sf::st_crs(las))
    sf_out[, c("TreeID", "Radius", "Error")]
  }

  structure(func,
            class       = c("spanner_map_mtd", "function"),
            method_name = "hough")
}
