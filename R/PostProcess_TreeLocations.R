#' Obtain tree information by processing point cloud data
#'
#' `process_tree_data` processes the output of `get_raster_eigen_treelocs` and `segment_graph` to add information
#' about the height, crown area, and diameter for each unique TreeID. It also has an optional parameter to return
#' an `sf` object representing the convex hulls for each tree.
#'
#' For terrestrial and mobile lidar datasets, tree locations and estimates of DBH are provided
#' by rasterizing individual point cloud values of relative neighborhood density (at 0.3 and
#' 1 m radius) and verticality within a slice of the normalized point cloud around breast height
#' (1.34 m). The algorithm then uses defined threshold values to classify the resulting rasters
#' and create unique polygons from the resulting classified raster. These point-density and
#' verticality polygons were selected by their intersection with one another, resulting in a
#' final set of polygons which were used to clip out regions of the point cloud that were most
#' likely to represent tree boles. A RANSAC cylinder fitting algorithm was then used to estimate
#' the fit of a cylinder to individual bole points. Cylinder centers and radius were used as inputs
#' to an individual tree segmentation.
#'
#' @param treeData An `sf` object containing the following tree information: `TreeID`,
#' `X`, `Y`, `Z`, `Radius`, and `Error`, output from the `get_raster_eigen_treelocs` function.
#' @param segmentedLAS A LAS object that is the output from `segment_graph`.
#' @param return_sf logical: If TRUE, returns an `sf` object representing the convex hulls for each tree.
#' @return sf object An updated `sf` object with the original columns plus:
#'   \describe{
#'     \item{height}{numeric: Height of the highest point for each TreeID.}
#'     \item{crown_area}{numeric: Area of the convex hull for each TreeID.}
#'     \item{diameter}{numeric: Diameter of the tree, calculated as twice the Radius.}
#'   }
#'   If `return_sf` is TRUE, also returns an `sf` object representing the convex hulls for each tree.
#'
#' @examples
#'
#' \dontrun{
#' # Set the number of threads to use in lidR
#' set_lidr_threads(8)
#'
#' # Download and read an example laz
#' getExampleData("DensePatchA")
#' LASfile = system.file("extdata", "DensePatchA.laz", package="spanner")
#' las = readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' # Don't forget to make sure the las object has a projection
#' projection(las) = sp::CRS("+init=epsg:26912")
#'
#' # Pre-process the example lidar dataset by classifying the ground points
#' # using lidR::csf(), normalizing it, and removing outlier points
#' # using lidR::ivf()
#' las = classify_ground(las, csf(sloop_smooth = FALSE,
#'                                 class_threshold = 0.5,
#'                                 cloth_resolution = 0.5, rigidness = 1L,
#'                                 iterations = 500L, time_step = 0.65))
#' las = normalize_height(las, tin())
#' las = classify_noise(las, ivf(0.25, 3))
#' las = filter_poi(las, Classification != LASNOISE)
#'
#' # Plot the non-ground points, colored by height
#' plot(filter_poi(las, Classification != 2), color = "Z")
#'
#' # Perform a deep inspection of the las object. If you see any
#' # red text, you may have issues!
#' las_check(las)
#'
#' # Find individual tree locations and attribute data
#' myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05,
#'                                        pt_spacing = 0.0254,
#'                                        dens_threshold = 0.2,
#'                                        neigh_sizes = c(0.333, 0.166, 0.5),
#'                                        eigen_threshold = 0.5,
#'                                        grid_slice_min = 0.6666,
#'                                        grid_slice_max = 2.0,
#'                                        minimum_polygon_area = 0.025,
#'                                        cylinder_fit_type = "ransac",
#'                                        max_dia = 0.5,
#'                                        SDvert = 0.25,
#'                                        n_pts = 20,
#'                                        n_best = 25)
#'
#' # Plot the tree information over a CHM
#' plot(lidR::grid_canopy(las, res = 0.2, p2r()))
#' points(myTreeLocs$X, myTreeLocs$Y, col = "black", pch = 16,
#'        cex = myTreeLocs$Radius^2 * 10, asp = 1)
#'
#' # Segment the point cloud
#' myTreeGraph = segment_graph(las = las, tree.locations = myTreeLocs, k = 50,
#'                              distance.threshold = 0.5,
#'                              use.metabolic.scale = FALSE,
#'                              ptcloud_slice_min = 0.6666,
#'                              ptcloud_slice_max = 2.0,
#'                              subsample.graph = 0.1,
#'                              return.dense = FALSE)
#'
#' # Plot it in 3D colored by treeID
#' plot(myTreeGraph, color = "treeID", pal=spanner_pal())
#'
#' # Process the data
#' processed_data <- process_tree_data(myTreeLocs, myTreeGraph, return_sf = TRUE)
#'
#' # Print the processed data
#' print(processed_data$data)
#' # Print the sf object if return_sf is TRUE
#' if (!is.null(processed_data$sf)) {
#'   print(processed_data$sf)
#' }
#' }
#'
#' @export
process_tree_data <- function(treeData, segmentedLAS, return_sf = FALSE) {
  # Extract unique tree IDs from segmentedLAS
  segmented_tree_ids <- unique(segmentedLAS$treeID)
  
  # Check if the same value and number of treeIDs exist in both treeData and segmentedLAS
  if (!all(treeData$TreeID %in% segmented_tree_ids) || length(unique(treeData$TreeID)) != length(segmented_tree_ids)) {
    stop("TreeIDs do not match between treeData and segmentedLAS.")
  }

  unique_tree_ids <- unique(treeData$TreeID)

  highest_points <- data.frame(TreeID = integer(), X = numeric(), Y = numeric(), Z = numeric(), Radius = numeric(), Error = numeric())
  convex_hulls <- list()
  crown_areas <- numeric()

  for (tree_id in unique_tree_ids) {
    tree_data <- subset(treeData, TreeID == tree_id)
    tree_las <- lidR::filter_poi(segmentedLAS, treeID == tree_id)

    # Find the highest point
    highest_point <- tree_las[which.max(tree_las$Z), ]
    highest_point_df <- data.frame(TreeID = tree_id, X = highest_point$X, Y = highest_point$Y, Z = highest_point$Z, Radius = NA, Error = NA)
    highest_points <- rbind(highest_points, highest_point_df)

    # Compute the convex hull
    tree_coords <- data.frame(X = tree_las$X, Y = tree_las$Y)
    if (nrow(tree_coords) > 2) {
      convex_hull <- sf::st_convex_hull(sf::st_union(sf::st_sfc(sf::st_multipoint(as.matrix(tree_coords)))))
      convex_hulls[[as.character(tree_id)]] <- convex_hull

      # Calculate the crown area
      crown_area <- sf::st_area(convex_hull)
      crown_areas <- c(crown_areas, crown_area)
    } else {
      convex_hulls[[as.character(tree_id)]] <- NULL
      crown_areas <- c(crown_areas, NA)
    }
  }

  # Add height, crown area, and diameter to the original data
  treeData$height <- NA
  treeData$crown_area <- NA
  treeData$diameter <- treeData$Radius * 2

  for (tree_id in unique_tree_ids) {
    treeData[treeData$TreeID == tree_id, "height"] <- highest_points[highest_points$TreeID == tree_id, "Z"]
    treeData[treeData$TreeID == tree_id, "crown_area"] <- crown_areas[unique_tree_ids == tree_id]
  }

  result <- list(data = treeData)

  if (return_sf) {
    hulls_sf <- do.call(rbind, lapply(names(convex_hulls), function(tree_id) {
      if (!is.null(convex_hulls[[tree_id]])) {
        sf::st_sf(TreeID = as.integer(tree_id), geometry = convex_hulls[[tree_id]])
      }
    }))
    result$sf <- hulls_sf
  }

  return(result)
}