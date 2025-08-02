#' Obtain tree information by rasterizing point cloud values of relative neighborhood density and
#' verticality within a slice of a normalized point cloud
#'
#' `get_raster_eigen_treelocs` returns a data.frame containing `TreeID`, `X`, `Y`, `Z`, `Radius`
#' and `Error` in the same units as the .las
#'
#' For terrestrial and mobile lidar datasets, tree locations and estimates of DBH are provided
#' by rasterizing individual point cloud values of relative neighborhood density (at 0.3 and
#' 1 m radius) and verticality within a slice of the normalized point cloud around breast height
#' (1.34 m). The algorithim then uses defined threshold values to classify the resulting rasters
#' and create unique polygons from the resulting classified raster. These point-density and
#' verticality polygons were selected by their intersection with one another, resulting in a
#' final set of polygons which were used to clip out regions of the point cloud that were most
#' likely to represent tree boles. A RANSAC cylinder fitting algorithm was then used to estimate
#' the fit of a cylinder to individual bole points. Cylinder centers and radius were used as inputs
#' to an individual tree segmentation
#'
#' @param las LAS Normalized las object.
#' @param res numeric Pixel width of rasterized point cloud metrics.
#' @param pt_spacing numeric Subsample spacing for graph connections.
#' @param dens_threshold numeric Minimum point density in raster cell to be considered as potential tree bole.
#' @param neigh_size numeric Vector for verticality and relative density (small and large neighborhoods) calculations
#' @param eigen_threshold numeric Minimum average verticality in raster cell to be considered as potential tree bole.
#' @param grid_slice_min numeric Lower bound of point cloud slice in normalized point cloud.
#' @param grid_slice_max numeric Upper bound of point cloud slice in normalized point cloud.
#' @param minimum_polygon_area numeric Smallest allowable polygon area of potential tree boles.
#' @param cylinder_fit_type  character Choose "ransac" or "irls" cylinder fitting.
#' @param max_dia numeric The max diameter (in m) of a resulting tree (use to eliminate commission errors).
#' @param SDvert numeric The standard deviation threshold below whihc polygons will be considered as tree boles.
#' @param n_best integer number of "best" ransac fits to keep when evaluating the best fit.
#' @param n_pts integer number of point to be selected per ransac iteraiton for fitting.
#' @param inliers integer expected proportion of inliers among cylinder points
#' @param conf numeric confidence level
#' @param max_angle numeric maximum tolerated deviation, in degrees, from vertical.
#' @param run_full_pipeline logical If TRUE (default), runs the complete tree detection pipeline. If FALSE, only performs initial processing and verticality calculation for testing purposes.
#' #' @return sf A sf object containing the following tree seed information: `TreeID`,
#' `Radius`, and `Error` in the same units as the .las, as well as the point geometry
#'
#' @examples
#'
#' \dontrun{
#' # set the number of threads to use in lidR
#' set_lidr_threads(8)
#'
#' # read the las (which must be downloaded with getExampleData())
#' LASfile <- system.file("extdata", "DensePatchA.laz", package="spanner")
#' las = readLAS(LASfile, select = "xyzc")
#'
#' # plot(las, color="Z", backend="lidRviewer", trim=30)
#'
#' # find tree locations and attribute data (optimized for dense patches)
#' myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.03, pt_spacing = 0.0254,
#'                                        dens_threshold = 0.1,
#'                                        neigh_sizes = c(0.2, 0.1, 0.3),
#'                                        eigen_threshold = 0.5,
#'                                        grid_slice_min = 0.5,
#'                                        grid_slice_max = 2.5,
#'                                        minimum_polygon_area = 0.01,
#'                                        cylinder_fit_type = "ransac",
#'                                        max_dia = 0.8,
#'                                        SDvert = 0.15,
#'                                        n_pts = 20,
#'                                        n_best = 25,
#'                                        inliers = 0.9,
#'                                        conf = 0.99,
#'                                        max_angle = 20)
#'
#' # For testing/debugging: run only initial processing (faster)
#' initial_results = get_raster_eigen_treelocs(las = las, res = 0.03,
#'                                            run_full_pipeline = FALSE)
#' print(initial_results)  # Shows verticality calculation summary
#'
#' plot(lidR::grid_canopy(las, res = 0.2, p2r()))
#' symbols(st_coordinates(myTreeLocs)[,1], st_coordinates(myTreeLocs)[,2],
#' circles = myTreeLocs$Radius^2, inches = FALSE, add = TRUE, bg = 'black')
#' }
#'
#' @export
get_raster_eigen_treelocs <- function(las = las, res = 0.03, pt_spacing = 0.0254, dens_threshold = 0.1,
                                      neigh_sizes=c(0.2, 0.1, 0.3), eigen_threshold = 0.5,
                                      grid_slice_min = 0.5, grid_slice_max = 2.5,
                                      minimum_polygon_area = 0.01, cylinder_fit_type = "ransac",
                                      output_location = getwd(), max_dia=0.8, SDvert = 0.15, n_best=25, n_pts=20,
                                      inliers = 0.9, conf = 0.99, max_angle = 20, run_full_pipeline = TRUE)
{

  ##---------------------- Preprocesssing -------------------------------------
  ## Subsample using systematic voxel grid to 1in
  message("Downsampling the scan... (step 1 of 14)\n")
  las <- lidR::decimate_points(las, random_per_voxel(res = pt_spacing, n = 1))
  if(is.null(las)) stop("No points in the las object after downsampling! Try increasing the point spacing.", call. = FALSE)

  ## Create the processing slice based on user's grid slice min/max
  message("Creating processing slice... (2/14)\n")
  slice_extra <- lidR::filter_poi(las, Z >= grid_slice_min, Z <= grid_slice_max)
  if(is.null(slice_extra)) stop("No points in the las object after processing slice! Try increasing the grid slice min/max.", call. = FALSE)

  message("Calculating verticality... (3/14)\n")
  # Use optimized geometric features function with automatic chunking for large datasets
  n_points <- nrow(slice_extra@data)
  cat("Processing", n_points, "points for verticality calculation...\n")

  vert_temp <- spanner:::C_geometric_features_simple(slice_extra, radius=neigh_sizes[1], max_neighbors=50)

  # Convert to data.table
  vert_temp <- data.table::as.data.table(vert_temp)
  if(nrow(vert_temp) == 0) stop("Problem calculating verticality...!", call. = FALSE)

  # vert_temp <- spanner::C_vert_in_sphere(slice_extra, radius = neigh_sizes[1], ncpu = lidR::get_lidr_threads())
  slice_extra <- lidR::add_lasattribute(slice_extra, as.numeric(vert_temp$verticality),
                                        name = "verticality", desc = "verticality")
  slice_extra@data$verticality[is.na(slice_extra@data$verticality)] <- 0.5

  # Check if user wants to run full pipeline or stop after initial processing
  if (!run_full_pipeline) {
    message("Stopping after verticality calculation (run_full_pipeline = FALSE)\n")
    message("Verticality calculation completed successfully!\n")
    message("Points processed: ", nrow(slice_extra@data), "\n")
    message("Verticality range: ", round(min(slice_extra@data$verticality, na.rm = TRUE), 3), " - ",
            round(max(slice_extra@data$verticality, na.rm = TRUE), 3), "\n")

    # Return a summary data frame instead of sf object
    summary_result <- data.frame(
      stage_completed = "verticality_calculation",
      points_processed = nrow(slice_extra@data),
      min_verticality = min(slice_extra@data$verticality, na.rm = TRUE),
      max_verticality = max(slice_extra@data$verticality, na.rm = TRUE),
      mean_verticality = mean(slice_extra@data$verticality, na.rm = TRUE),
      processing_time = as.character(Sys.time()),
      message = "Initial processing completed. Set run_full_pipeline = TRUE to continue with tree detection."
    )

    return(summary_result)
  }

  message("Creating canopy cover raster... (4/14)\n")
  cancov <- lidR::grid_canopy(las, 0.1, p2r(0.1))

  message("Calculating relative density... (5/14)\n")
  cnt_local <- spanner:::C_count_in_disc(X = slice_extra@data$X, Y = slice_extra@data$Y,
                                         x = slice_extra@data$X, y = slice_extra@data$Y,
                                         radius = neigh_sizes[2], ncpu = lidR::get_lidr_threads())
  cnt_large <- spanner:::C_count_in_disc(X = slice_extra@data$X, Y = slice_extra@data$Y,
                                         x = slice_extra@data$X, y = slice_extra@data$Y,
                                         radius = neigh_sizes[3], ncpu = lidR::get_lidr_threads())
  slice_extra <- lidR::add_lasattribute(slice_extra, as.numeric(cnt_local/cnt_large),
                                        name = "relative_density", desc = "relative density")

  ##---------------------- Eigen & Raster Processing -------------------------------------
  message("Gridding relative density... (6/14)\n")

  density_grid <- terra::rast(lidR::grid_metrics(slice_extra, ~stats::median(relative_density, na.rm = T), res = res,
                                                 start = c(min(slice_extra@data$X), min(slice_extra@data$Y))))
  density_polygon <- terra::as.polygons(terra::classify(density_grid, rbind(c(0,dens_threshold,0),
                                                                            c(dens_threshold,1,1))), dissolve = T)
  density_polygon <- terra::disagg(density_polygon[density_polygon$V1 > 0,])
  density_polygon <- sf::st_as_sf(density_polygon)
  density_polygon$area = as.numeric(sf::st_area(density_polygon))
  density_polygon <- density_polygon[density_polygon$area > minimum_polygon_area,]
  density_polygon <- sfheaders::sf_remove_holes(density_polygon)
  density_polygon <- sf::st_buffer(sf::st_buffer(density_polygon, 1.5), -1.475) ## smoothing out the resulting polygons
  if(nrow(density_polygon) == 0 | is.null(density_polygon)) stop("No density polygons were created from the rasterized point cloud metrics! Try adjusting the threshold values.", call. = FALSE)

  message("Gridding verticality... (7/14)\n")

  verticality_grid <- terra::rast(lidR::grid_metrics(slice_extra, ~stats::quantile(verticality, 0.5, na.rm = T), res = res))
  verticality_polygon <- terra::as.polygons(terra::classify(verticality_grid, rbind(c(0,eigen_threshold,0),
                                                                                    c(eigen_threshold,1,1))), dissolve = T)
  verticality_polygon <- terra::disagg(verticality_polygon[verticality_polygon$V1 > 0,])
  verticality_polygon <- sf::st_as_sf(verticality_polygon)
  verticality_polygon$area = as.numeric(sf::st_area(verticality_polygon))
  verticality_polygon <- verticality_polygon[verticality_polygon$area > minimum_polygon_area,]
  verticality_polygon <- sfheaders::sf_remove_holes(verticality_polygon)
  verticality_polygon <- sf::st_buffer(sf::st_buffer(verticality_polygon, 1.5), -1.475) ## smoothing out the resulting polygons
  if(nrow(verticality_polygon) == 0 | is.null(verticality_polygon)) stop("No verticality polygons were created from the rasterized point cloud metrics! Try adjusting the threshold values.", call. = FALSE)

  message("Merging relative density and verticality... (8/14)\n")
  merged <- sf::st_intersection(verticality_polygon, density_polygon)
  merged$area <- as.numeric(sf::st_area(merged))
  sf::st_crs(merged) = lidR::crs(las)
  merged <- merged[merged$area > minimum_polygon_area,]
  if(nrow(merged) == 0 | is.null(merged)) stop("No merged polygons were created from the rasterized point cloud metrics! Try adjusting the threshold values.", call. = FALSE)

  message("Filtering merged polygon... (9/14)\n")
  ## sd verticality
  verticalitySD_grid <- terra::rast(grid_metrics(slice_extra, ~stats::sd(verticality), res = res))
  merged$SDvert <- round(terra::extract(verticalitySD_grid, terra::vect(sf::as_Spatial(merged)),
                                        fun = function(x){stats::median(x, na.rm = T)})[,2],2)
  ## some simple filtering: lower canopy
  Z01 <- terra::rast(grid_metrics(slice_extra, res = res, ~stats::quantile(Z, 0.01)))
  merged$lowHGT <- round(terra::extract(Z01, terra::vect(sf::as_Spatial(merged)),
                                        fun = function(x){stats::quantile(x, 0.05, na.rm = T)})[,2],2)
  Z99 <- terra::rast( grid_metrics(slice_extra, res = res, ~stats::quantile(Z, 0.99)) )
  merged$highHGT <- round(terra::extract(Z99, terra::vect(sf::as_Spatial(merged)),
                                         fun = function(x){stats::quantile(x, 0.99, na.rm = T)})[,2],2)
  Z50 <- terra::rast(grid_metrics(slice_extra, res = res, ~stats::quantile(Z, 0.5)))
  merged$medHGT <- round(terra::extract(Z50, terra::vect(sf::as_Spatial(merged)),
                                        fun = function(x){stats::median(x, na.rm = T)})[,2],2)
  merged$diffHGT <- round(merged$highHGT - merged$lowHGT,2)
  # merged <- merged[merged$diffHGT > ((grid_slice_max - grid_slice_min) / 2), ]
  merged$maxHGT <- round(terra::extract(terra::rast(cancov), terra::vect(sf::as_Spatial(merged)), fun = function(x) max(x, na.rm = TRUE))[,2], 2)

  merged <- merged[merged$SDvert < SDvert, ]
  merged <- merged[merged$diffHGT > (grid_slice_max-grid_slice_min)*0.5, ]
  # ## low and high within the slice
  merged <- merged[merged$lowHGT < (grid_slice_min*1.25), ]
  #merged <- merged[merged$medHGT <= grid_slice_max*0.75, ]
  merged <- merged[merged$maxHGT >= grid_slice_max, ]

  message("Obtaining polygon attributes...(10/14)\n")
  merged$TreeID = sample(1:nrow(merged), nrow(merged), replace=F)

  message("Processing the resulting clipped slice...(11/14)\n")
  coords <- data.frame(sf::st_coordinates(na.omit(merged)))
  circles <- list()
  for(id in unique(coords$L2)){
    circles[[id]] <- data.frame(conicfit::CircleFitByLandau(coords[coords$L2 == id, c("X","Y")]))
    names(circles[[id]]) <- c("X","Y","R")
  }
  circles <- dplyr::bind_rows((circles))
  circles_sf <- sf::st_sf(sf::st_buffer(sf::st_cast(sf::st_sfc(sf::st_multipoint(as.matrix(circles)[,1:2])),
                                                    to = "POINT"), circles$R))
  circles_sf$R <- circles$R
  circles_sf$TreeID <- sample(1:nrow(circles_sf), nrow(circles_sf), replace=F)
  circles_sf <- sf::st_buffer(circles_sf, 0.075)
  sf::st_crs(circles_sf) = sf::st_crs(slice_extra)
  slice_clip <- lidR::merge_spatial(las = lidR::clip_roi(slice_extra, sf::st_sf(sf::st_union(circles_sf))),
                                    source = circles_sf, attribute = "TreeID")

  slice_clip <- lidR::filter_poi(slice_clip, verticality >= eigen_threshold)
  if(is.null(slice_clip)) stop("No points in the las object after processing the resulting clipped slice! Try adjusting the threshold values.", call. = FALSE)

  ##---------------------- Identify Trees -------------------------------------
  message("Fitting nested height (length) cylinders...(12/14)\n")
  cyl_fit <- list()
  for(t in 1:length(sort(unique(slice_clip$TreeID))))
  {
    min = grid_slice_min; max = grid_slice_max
    if(cylinder_fit_type == "ransac"){
      n_pts = n_pts
      n_best = n_best
    } else if(cylinder_fit_type == "irls"){
      n_pts = min(data.frame(lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data %>%
                               dplyr::group_by(TreeID) %>%
                               dplyr::summarize(count = length(X)))$count)/2
      n_best = 100
    } else {
      stop("You will need to choose either 'ransac' or 'irls' cylinder fitting method...", call. = FALSE)
    }
    fit <- spanner:::cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t]),
                                 method = cylinder_fit_type, n = n_pts, inliers = inliers,
                                 conf = conf, max_angle = max_angle, n_best = n_best)
    fit$TreeID <- sort(unique(slice_clip$TreeID))[t]
    fit$dbh_width <- max-min

    cyl_fit[[t]] <- fit

    min = 1.1 - ((1.1 - grid_slice_min)/2); max =  1.6 + ((grid_slice_max - 1.6)/2)
    if(cylinder_fit_type == "ransac"){
      n_pts = n_pts
      n_best = n_best
    } else if(cylinder_fit_type == "irls"){
      n_pts = min(data.frame(lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data %>%
                               dplyr::group_by(TreeID) %>%
                               dplyr::summarize(count = length(X)))$count)/2
      n_best = 100
    } else {
      stop("You will need to choose either 'ransac' or 'irls' cylinder fitting method...", call. = FALSE)
    }
    fit <- spanner:::cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t],
                                                  Z >= min, Z <= max),
                                 method = cylinder_fit_type, n = n_pts, inliers = inliers,
                                 conf = conf, max_angle = max_angle, n_best = n_best)
    fit$TreeID <- sort(unique(slice_clip$TreeID))[t]
    fit$dbh_width <- max-min

    cyl_fit[[t+length(unique(slice_clip$TreeID))]] <- fit

    min = 1.1; max = 1.6
    if(cylinder_fit_type == "ransac"){
      n_pts = n_pts
      n_best = n_best
    } else if(cylinder_fit_type == "irls"){
      n_pts = min(data.frame(lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data %>%
                               dplyr::group_by(TreeID) %>%
                               dplyr::summarize(count = length(X)))$count)/2
      n_best = 100
    } else {
      stop("You will need to choose either 'ransac' or 'irls' cylinder fitting method...", call. = FALSE)
    }
    fit <- spanner:::cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t],
                                                  Z >= min, Z <= max),
                                 method = cylinder_fit_type, n = n_pts, inliers = inliers,
                                 conf = conf, max_angle = max_angle, n_best = n_best)
    fit$TreeID <- sort(unique(slice_clip$TreeID))[t]
    fit$dbh_width <- max-min

    cyl_fit[[t + ((length(unique(slice_clip$TreeID))) * 2)]] <- fit
  }
  (cyl_fit <- dplyr::bind_rows(cyl_fit))
  cyl_fit <- cyl_fit[cyl_fit$radius <= max_dia, ]
  cyl_fit <- cyl_fit[complete.cases(cyl_fit),]
  if(nrow(cyl_fit) == 0 | is.null(cyl_fit)) stop("No cylinders were fit to the point cloud! Try adjusting the threshold values.", call. = FALSE)
  summary_cyl_fit <- cyl_fit %>% dplyr::group_by(TreeID) %>% dplyr::summarize_all(mean)
  if(nrow(summary_cyl_fit) == 0 | is.null(summary_cyl_fit)) stop("No summary data was created from the fitted cylinders! Try adjusting the threshold values.", call. = FALSE)

  ##---------------------- Cleaning up the Output -------------------------------------
  message("Successfully obtained the cylinder summaries... (13/14)\n")
  output<-summary_cyl_fit[,c("TreeID","px","py","pz","radius","err")]
  colnames(output)<-c("TreeID","X","Y","Z","Radius","Error")
  if(nrow(output) == 0 | is.null(output)) stop("No output data was created from the fitted cylinders! Try adjusting the threshold values.", call. = FALSE)
  return(st_as_sf(output, coords = c("X", "Y", "Z"), crs = lidR::st_crs(las)))
  message("Done! (14/14)\n")

  ##---------------------- Secondary Segmentation for High Error Trees ------
  message("Checking for under-segmented trees using RANSAC error... (13.5/14)\n")

  # Define error threshold for detecting potentially multi-tree segments
  error_threshold <- quantile(summary_cyl_fit$err, 0.75, na.rm = TRUE)  # Top 25% of errors
  high_error_trees <- summary_cyl_fit[summary_cyl_fit$err > error_threshold, ]

  if(nrow(high_error_trees) > 0) {
    message("Found ", nrow(high_error_trees), " trees with high RANSAC error - attempting secondary segmentation...")

    # Store original results
    original_output <- summary_cyl_fit
    additional_trees <- list()

    for(i in 1:nrow(high_error_trees)) {
      tree_id <- high_error_trees$TreeID[i]
      message("  Analyzing tree ", tree_id, " (error: ", round(high_error_trees$err[i], 4), ")")

      # Extract points for this high-error tree
      tree_points <- lidR::filter_poi(slice_clip, TreeID == tree_id)

      if(nrow(tree_points@data) > 50) {  # Need sufficient points for secondary segmentation

        # Try more aggressive clustering using smaller resolution
        tree_coords <- tree_points@data[, c("X", "Y", "Z")]

        # Create a finer grid for this specific tree area
        tree_bbox <- list(
          xmin = min(tree_coords$X) - 0.1,
          xmax = max(tree_coords$X) + 0.1,
          ymin = min(tree_coords$Y) - 0.1,
          ymax = max(tree_coords$Y) + 0.1
        )

        # Use 2/3 of the original resolution for finer detection
        fine_res <- res * 0.67

        # Create density raster for just this tree area
        tree_las_temp <- tree_points
        tree_density <- grid_density(tree_las_temp, res = fine_res)

        # Lower threshold for secondary detection
        secondary_dens_threshold <- dens_threshold * 0.7

        # Find high-density cells
        high_dens_cells <- terra::as.data.frame(tree_density, xy = TRUE)
        high_dens_cells <- high_dens_cells[high_dens_cells$density >= secondary_dens_threshold, ]

        if(nrow(high_dens_cells) > 2) {
          # Apply clustering to find potential tree centers
          if(nrow(high_dens_cells) > 4) {
            # Use k-means with k=2 or k=3 based on point distribution
            max_clusters <- min(3, floor(nrow(high_dens_cells) / 2))

            if(max_clusters >= 2) {
              # Try clustering
              clusters <- kmeans(high_dens_cells[, c("x", "y")],
                                centers = max_clusters,
                                nstart = 10)

              # Check if clusters are spatially separated
              cluster_centers <- clusters$centers
              min_distance <- min(dist(cluster_centers))

              # Only proceed if clusters are separated by at least 0.3m
              if(min_distance > 0.3) {
                message("    Found ", max_clusters, " potential sub-trees")

                # Assign points to clusters based on proximity
                tree_coords$cluster <- clusters$cluster[
                  apply(tree_coords[, c("X", "Y")], 1, function(point) {
                    distances <- apply(high_dens_cells[, c("x", "y")], 1, function(cell) {
                      sqrt(sum((point - cell)^2))
                    })
                    which.min(distances)
                  })
                ]

                tree_coords$cluster <- clusters$cluster[
                  apply(tree_coords[, c("X", "Y")], 1, function(point) {
                    distances <- apply(high_dens_cells[, c("x", "y")], 1, function(cell) {
                      sqrt(sum((point - cell)^2))
                    })
                    if(min(distances) < fine_res * 2) {
                      high_dens_cells$cluster[which.min(distances)]
                    } else {
                      NA  # Point too far from any high-density cell
                    }
                  })
                ]

                # Remove points not assigned to any cluster
                tree_coords <- tree_coords[!is.na(tree_coords$cluster), ]

                # Fit cylinders to each cluster
                for(cluster_id in unique(tree_coords$cluster)) {
                  cluster_points <- tree_coords[tree_coords$cluster == cluster_id, ]

                  if(nrow(cluster_points) >= 20) {  # Need minimum points for cylinder fitting
                    # Create temporary LAS object for this cluster
                    cluster_las <- LAS(cluster_points[, c("X", "Y", "Z")])

                    # Fit cylinder using same parameters as main function
                    cluster_fit <- spanner:::cylinderFit(
                      cluster_las,
                      method = cylinder_fit_type,
                      n = n_pts,
                      inliers = inliers,
                      conf = conf,
                      max_angle = max_angle,
                      n_best = n_best
                    )

                    if(!is.null(cluster_fit) && !is.na(cluster_fit$err)) {
                      # Check if this sub-tree has better (lower) error than original
                      if(cluster_fit$err < high_error_trees$err[i] * 0.8) {  # 20% improvement
                        cluster_fit$TreeID <- max(c(summary_cyl_fit$TreeID,
                                                   sapply(additional_trees, function(x) x$TreeID)),
                                                 na.rm = TRUE) + 1
                        cluster_fit$dbh_width <- grid_slice_max - grid_slice_min

                        # Check size constraints
                        if(cluster_fit$radius <= max_dia && cluster_fit$radius > 0.05) {
                          additional_trees[[length(additional_trees) + 1]] <- cluster_fit
                          message("      Added sub-tree with radius ", round(cluster_fit$radius, 3),
                                 " and error ", round(cluster_fit$err, 4))
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    # Combine original results with additional trees
    if(length(additional_trees) > 0) {
      additional_df <- dplyr::bind_rows(additional_trees)

      # Remove original high-error trees that were successfully split
      successfully_split <- sapply(additional_trees, function(x) {
        # Check if any additional tree is close to the original high-error tree
        orig_tree <- high_error_trees[high_error_trees$TreeID %in%
                                     summary_cyl_fit$TreeID, ]
        if(nrow(orig_tree) > 0) {
          distance <- sqrt((x$px - orig_tree$px[1])^2 + (x$py - orig_tree$py[1])^2)
          return(distance < 1.0)  # Within 1 meter
        }
        return(FALSE)
      })

      if(any(successfully_split)) {
        # Remove trees that were successfully split
        split_tree_ids <- unique(sapply(additional_trees[successfully_split],
                                       function(x) high_error_trees$TreeID[
                                         which.min(sqrt((x$px - high_error_trees$px)^2 +
                                                        (x$py - high_error_trees$py)^2))
                                       ]))

        summary_cyl_fit <- summary_cyl_fit[!summary_cyl_fit$TreeID %in% split_tree_ids, ]
        message("  Removed ", length(split_tree_ids), " original high-error trees")
      }

      # Add new trees
      summary_cyl_fit <- dplyr::bind_rows(summary_cyl_fit, additional_df)
      message("  Added ", nrow(additional_df), " new trees from secondary segmentation")
    }
  }
}
