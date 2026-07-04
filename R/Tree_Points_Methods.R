# ============================================================================
# spanner v2.0 – Tree_Points_Methods.R
#
# Method factories for the tree_points() stage.
# Each factory returns a closure of class "spanner_tpt_mtd" that accepts
# (las, map) and returns a LAS with a treeID integer column.
# ============================================================================


# ----------------------------------------------------------------------------
#' Tree-points method: Voronoi nearest-tree partition
#'
#' Assigns every point to the closest tree centre in 2-D Euclidean space.
#' This is a complete spatial partition — every point receives a tree ID
#' unless `max_dist` is set.
#'
#' Implemented via `C_tree_ids_voronoi()`.  For typical TLS/MLS scenes
#' (< a few hundred trees) the brute-force nearest-neighbour search is fast
#' enough.  For very large mobile-mapping datasets tile the cloud first.
#'
#' @param max_dist numeric.  Points farther than this distance (m) from every
#'   tree centre are left with `treeID = 0`.  Set to `-1` (default) to assign
#'   **all** points — a pure Voronoi partition.
#'
#' @return A closure of class `"spanner_tpt_mtd"` / `"function"` with
#'   attribute `method_name = "voronoi"`.
#'
#' @seealso [tree_points()], [assign_crop()], [assign_graph()]
#'
#' @examples
#' \donttest{
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#'
#' map   <- find_trees(las, use_raster_eigen())
#' las_v <- tree_points(las, map, method = assign_voronoi())
#' }
#'
#' @export
assign_voronoi <- function(max_dist = -1.0) {
  stopifnot(is.numeric(max_dist), length(max_dist) == 1)

  func <- function(las, map) {
    tl_xy  <- sf::st_coordinates(map)
    tl_ids <- as.integer(map$TreeID)

    ids <- C_tree_ids_voronoi(
      pt_x     = las@data$X,
      pt_y     = las@data$Y,
      tree_x   = tl_xy[, 1],
      tree_y   = tl_xy[, 2],
      tree_ids = tl_ids,
      max_dist = max_dist
    )

    las@data$treeID <- ids
    las
  }

  structure(func,
            class       = c("spanner_tpt_mtd", "function"),
            method_name = "voronoi")
}


# ----------------------------------------------------------------------------
#' Tree-points method: fixed-radius crop
#'
#' Assigns points within a fixed circle (or square) around each tree centre
#' to that tree.  Points outside every crop receive `treeID = 0`.  Where
#' crops overlap, priority is given to the tree that appears first when
#' iterating in `map` row order.
#'
#' Implemented via `C_tree_ids_crop()`.  Useful for conservative stem-only
#' point extraction where crown points should be excluded.
#'
#' @param l      numeric.  Circle radius (metres) when `circle = TRUE`, or
#'   half side-length when `circle = FALSE`.  Default `1.0`.
#' @param circle logical.  `TRUE` = circular crop (default); `FALSE` = square
#'   crop.
#'
#' @return A closure of class `"spanner_tpt_mtd"` / `"function"` with
#'   attribute `method_name = "crop"`.
#'
#' @seealso [tree_points()], [assign_voronoi()], [assign_graph()]
#'
#' @examples
#' \donttest{
#' map   <- find_trees(las, use_raster_eigen())
#' las_c <- tree_points(las, map, method = assign_crop(l = 1.5))
#' }
#'
#' @export
assign_crop <- function(l = 1.0, circle = TRUE) {
  stopifnot(is.numeric(l), length(l) == 1, l > 0,
            is.logical(circle), length(circle) == 1)

  func <- function(las, map) {
    tl_xy  <- sf::st_coordinates(map)
    tl_ids <- as.integer(map$TreeID)

    ids <- C_tree_ids_crop(
      pt_x     = las@data$X,
      pt_y     = las@data$Y,
      tree_x   = tl_xy[, 1],
      tree_y   = tl_xy[, 2],
      tree_ids = tl_ids,
      length   = l,
      circle   = circle
    )

    las@data$treeID <- ids
    las
  }

  structure(func,
            class       = c("spanner_tpt_mtd", "function"),
            method_name = "crop")
}


# ----------------------------------------------------------------------------
#' Tree-points method: graph propagation
#'
#' Wraps [segment_graph()] as a composable method object for [tree_points()].
#' Builds a k-nearest-neighbour graph over the point cloud and propagates
#' tree IDs from seed locations via shortest paths — the most spatially
#' accurate but slowest assignment method.
#'
#' All arguments are forwarded verbatim to [segment_graph()].  Refer to that
#' function's documentation for the full parameter list including `k`,
#' `distance.threshold`, `use.metabolic.scale`, `subsample.graph`,
#' `return.dense`, etc.
#'
#' @param ... Named arguments forwarded to [segment_graph()].
#'
#' @return A closure of class `"spanner_tpt_mtd"` / `"function"` with
#'   attribute `method_name = "graph"`.
#'
#' @seealso [tree_points()], [assign_voronoi()], [assign_crop()],
#'   [segment_graph()]
#'
#' @examples
#' \donttest{
#' map   <- find_trees(las, use_raster_eigen())
#' las_g <- tree_points(las, map,
#'            method = assign_graph(k = 50,
#'                                  distance.threshold = 0.5,
#'                                  subsample.graph = 0.1,
#'                                  return.dense = TRUE))
#' }
#'
#' @export
assign_graph <- function(...) {
  args <- list(...)
  func <- function(las, map) {
    do.call(segment_graph, c(list(las = las, tree.locations = map), args))
  }
  structure(func,
            class       = c("spanner_tpt_mtd", "function"),
            method_name = "graph")
}
