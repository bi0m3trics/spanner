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
#' @param neigh_sizes numeric Vector for verticality and relative density (small and large neighborhoods) calculations
#' @param eigen_threshold numeric Minimum average verticality in raster cell to be considered as potential tree bole.
#' @param grid_slice_min numeric Lower bound of point cloud slice in normalized point cloud.
#' @param grid_slice_max numeric Upper bound of point cloud slice in normalized point cloud.
#' @param minimum_polygon_area numeric Smallest allowable polygon area of potential tree boles.
#' @param cylinder_fit_type  character Choose "ransac" or "irls" cylinder fitting.
#' @param max_dia numeric The max diameter (in m) of a resulting tree (use to eliminate commission errors).
#' @param SDvert numeric The standard deviation threshold below which polygons will be considered as tree boles.
#' @param n_best integer number of "best" ransac fits to keep when evaluating the best fit.
#' @param n_pts integer number of point to be selected per ransac iteraiton for fitting.
#' @param inliers integer expected proportion of inliers among cylinder points
#' @param conf numeric confidence level
#' @param max_angle numeric maximum tolerated deviation, in degrees, from vertical.
#' #' @return sf A sf object containing the following tree seed information: `TreeID`,
#' `Radius`, and `Error` in the same units as the .las, as well as the point geometry
#'
#' @examples
#'
#' \dontrun{
#' # Set the number of threads to use in lidR
#' set_lidr_threads(8)
#'
#' LASfile = system.file("extdata", "TLS_Clip.laz", package="spanner")
#' las = readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
#' # Don't forget to make sure the las object has a projection
#' sf::st_crs(las) <- 26912
#'
#' # Pre-process the example lidar dataset by classifying the ground  and noise points
#' # using lidR::csf(), normalizing it, and removing outlier points
#' # using lidR::ivf()
#' # las = classify_ground(las, csf(sloop_smooth = FALSE,
#' #                                 class_threshold = 0.5,
#' #                                cloth_resolution = 0.5, rigidness = 1L,
#' #                                 iterations = 500L, time_step = 0.65))
#' # las = normalize_height(las, tin())
#' # las = classify_noise(las, ivf(0.25, 3))
#' # las = filter_poi(las, Classification != LASNOISE)
#'
#' # Plot the non-ground points, colored by height
#' # plot(filter_poi(las, Classification != 2), color = "Z")
#'
#' # find tree locations and attribute data
#' myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.025, pt_spacing = 0.0254,
#'                                        dens_threshold = 0.25,
#'                                        neigh_sizes = c(0.25, 0.15, 0.66),
#'                                        eigen_threshold = 0.75,
#'                                        grid_slice_min = 1,
#'                                        grid_slice_max = 2,
#'                                        minimum_polygon_area = 0.005,
#'                                        cylinder_fit_type = "ransac",
#'                                        max_dia = 1,
#'                                        SDvert = 0.33,
#'                                        n_pts = 20,
#'                                        n_best = 25,
#'                                        inliers = 0.9,
#'                                        conf = 0.99,
#'                                        max_angle = 20)
#'
#' # Plot results if trees were found
#' if (!is.null(myTreeLocs) && nrow(myTreeLocs) > 0) {
#'   plot(lidR::rasterize_canopy(las, res = 0.2, p2r()))
#'   symbols(sf::st_coordinates(myTreeLocs)[,1], sf::st_coordinates(myTreeLocs)[,2],
#'           circles = myTreeLocs$Radius^2*3.14, inches = FALSE, add = TRUE, bg = 'black')
#' } else {
#'   message("No tree locations were found. Try adjusting the parameters.")
#' }
#' }
#'
#' @export
get_raster_eigen_treelocs <- function(las = las, res = 0.05, pt_spacing = 0.0254, dens_threshold = 0.2,
                                      neigh_sizes=c(0.333, 0.166, 0.5), eigen_threshold = 0.6666,
                                      grid_slice_min = 0.6666, grid_slice_max = 2.0,
                                      minimum_polygon_area = 0.025, cylinder_fit_type = "ransac",
                                      max_dia=0.5, SDvert = 0.25, n_best=25, n_pts=20, inliers = 0.9,
                                      conf = 0.99, max_angle = 20)
{

  ##---------------------- Preprocesssing -------------------------------------
  ## Subsample using systematic voxel grid to 1in
  message("Downsampling the scan... (step 1 of 14)\n")
  las <- lidR::decimate_points(las, random_per_voxel(res = pt_spacing, n = 1))
  if(lidR::is.empty(las)) stop("No points in the las object after downsampling! Try increasing the point spacing.", call. = FALSE)

  ## Create the processing slice based on user's grid slice min/max
  message("Creating processing slice... (2/14)\n")
  slice_extra <- lidR::filter_poi(las, Z >= grid_slice_min, Z <= grid_slice_max)
  if(lidR::is.empty(slice_extra)) stop("No points in the las object after processing slice! Try increasing the slice min/max.", call. = FALSE)

  message("Calculating verticality... (3/14)\n")
  vert_temp <- eigen_metrics(slice_extra, radius=neigh_sizes[1], ncpu=lidR::get_lidr_threads())
  if(nrow(vert_temp) == 0) stop("Problemn calculating verticality...!", call. = FALSE)

  # vert_temp <- spanner::C_vert_in_sphere(slice_extra, radius = neigh_sizes[1], ncpu = lidR::get_lidr_threads())
  slice_extra <- lidR::add_lasattribute(slice_extra, as.numeric(vert_temp$Verticality),
                                        name = "verticality", desc = "verticality")
  slice_extra@data$verticality[is.na(slice_extra@data$verticality)] <- 0.5

  message("Creating canopy height raster... (4/14)\n")
  cancov <- lidR::pixel_metrics(las, ~max(Z), res = 0.1)

  message("Calculating relative density... (5/14)\n")
  cnt_local <- C_count_in_disc(X = slice_extra@data$X, Y = slice_extra@data$Y,
                               x = slice_extra@data$X, y = slice_extra@data$Y,
                               radius = neigh_sizes[2], ncpu = lidR::get_lidr_threads())
  cnt_large <- C_count_in_disc(X = slice_extra@data$X, Y = slice_extra@data$Y,
                               x = slice_extra@data$X, y = slice_extra@data$Y,
                               radius = neigh_sizes[3], ncpu = lidR::get_lidr_threads())
  slice_extra <- lidR::add_lasattribute(slice_extra, as.numeric(cnt_local/cnt_large),
                                        name = "relative_density", desc = "relative density")

  ##---------------------- Eigen & Raster Processing -------------------------------------
  message("Gridding relative density... (6/14)\n")

  density_grid <- lidR::pixel_metrics(slice_extra, ~stats::median(relative_density, na.rm = T), res = res,
                                                 start = c(min(slice_extra@data$X), min(slice_extra@data$Y)))
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

  verticality_grid <- lidR::pixel_metrics(slice_extra, ~stats::quantile(verticality, 0.5, na.rm = T), res = res)
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
  sf::st_crs(merged) = sf::st_crs(las)
  merged <- merged[merged$area > minimum_polygon_area,]
  if(nrow(merged) == 0 | is.null(merged)) stop("No merged polygons were created from the rasterized point cloud metrics! Try adjusting the threshold values.", call. = FALSE)

  message("Filtering merged polygon... (9/14)\n")
  ## sd verticality
  verticalitySD_grid <- lidR::pixel_metrics(slice_extra, ~stats::sd(verticality), res = res)
  merged$SDvert <- round(terra::extract(verticalitySD_grid, terra::vect(merged),
                                        fun = function(x){stats::median(x, na.rm = T)})[,2],2)
  ## some simple filtering: lower canopy
  Z01 <- lidR::pixel_metrics(slice_extra, res = res, ~stats::quantile(Z, 0.01))
  merged$lowHGT <- round(terra::extract(Z01, terra::vect(merged),
                                        fun = function(x){stats::quantile(x, 0.05, na.rm = T)})[,2],2)
  Z99 <- lidR::pixel_metrics(slice_extra, res = res, ~stats::quantile(Z, 0.99))
  merged$highHGT <- round(terra::extract(Z99, terra::vect(merged),
                                         fun = function(x){stats::quantile(x, 0.99, na.rm = T)})[,2],2)
  Z50 <- lidR::pixel_metrics(slice_extra, res = res, ~stats::quantile(Z, 0.5))
  merged$medHGT <- round(terra::extract(Z50, terra::vect(merged),
                                        fun = function(x){stats::median(x, na.rm = T)})[,2],2)
  merged$diffHGT <- round(merged$highHGT - merged$lowHGT,2)
  # merged <- merged[merged$diffHGT > ((grid_slice_max - grid_slice_min) / 2), ]
  merged$maxHGT <- round(terra::extract(cancov, terra::vect(merged), fun = function(x) max(x, na.rm = TRUE))[,2], 2)

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
  if (length(circles) == 0) {
    message(" No stem detection for this slice.")
    return(NULL) # stop function and return NULL
  }
  circles <- dplyr::bind_rows((circles))
  circles_sf <- sf::st_sf(sf::st_buffer(sf::st_cast(sf::st_sfc(sf::st_multipoint(as.matrix(circles)[,1:2, drop = FALSE])),
                                                    to = "POINT"), circles$R))
  circles_sf$R <- circles$R
  circles_sf$TreeID <- sample(1:nrow(circles_sf), nrow(circles_sf), replace=F)
  circles_sf <- sf::st_buffer(circles_sf, 0.075)
  sf::st_crs(circles_sf) = sf::st_crs(slice_extra)
  slice_clip <- lidR::merge_spatial(las = lidR::clip_roi(slice_extra, sf::st_sf(sf::st_union(circles_sf))),
                                    source = circles_sf, attribute = "TreeID")

  slice_clip <- lidR::filter_poi(slice_clip, verticality >= eigen_threshold)
  if(lidR::is.empty(slice_clip)) stop("No points in the las object after processing the resulting clipped slice! Try adjusting the threshold values.", call. = FALSE)

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
      tmp_dat <- lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data
      tmp_dat_grouped <- dplyr::group_by(tmp_dat, TreeID)
      tmp_dat_summary <- dplyr::summarize(tmp_dat_grouped, count = length(X))
      n_pts = min(tmp_dat_summary$count)/2
      n_best = 100
    } else {
      stop("You will need to choose either 'ransac' or 'irls' cylinder fitting method...", call. = FALSE)
    }
    fit <- cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t]),
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
      tmp_dat <- lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data
      tmp_dat_grouped <- dplyr::group_by(tmp_dat, TreeID)
      tmp_dat_summary <- dplyr::summarize(tmp_dat_grouped, count = length(X))
      n_pts = min(tmp_dat_summary$count)/2
      n_best = 100
    } else {
      stop("You will need to choose either 'ransac' or 'irls' cylinder fitting method...", call. = FALSE)
    }
    fit <- cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t],
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
      tmp_dat <- lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data
      tmp_dat_grouped <- dplyr::group_by(tmp_dat, TreeID)
      tmp_dat_summary <- dplyr::summarize(tmp_dat_grouped, count = length(X))
      n_pts = min(tmp_dat_summary$count)/2
      n_best = 100
    } else {
      stop("You will need to choose either 'ransac' or 'irls' cylinder fitting method...", call. = FALSE)
    }
    fit <- cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t],
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
  summary_cyl_fit <- dplyr::summarize_all(dplyr::group_by(cyl_fit, TreeID), mean)
  if(nrow(summary_cyl_fit) == 0 | is.null(summary_cyl_fit)) stop("No summary data was created from the fitted cylinders! Try adjusting the threshold values.", call. = FALSE)

  ##---------------------- Cleaning up the Output -------------------------------------
  message("Successfully obtained the cylinder summaries... (13/14)\n")
  output<-summary_cyl_fit[,c("TreeID","px","py","pz","radius","err")]
  colnames(output)<-c("TreeID","X","Y","Z","Radius","Error")
  if(nrow(output) == 0 | is.null(output)) stop("No output data was created from the fitted cylinders! Try adjusting the threshold values.", call. = FALSE)
  message("Done! (14/14)\n")
  return(st_as_sf(output, coords = c("X", "Y", "Z"), crs = lidR::st_crs(las)))
}
