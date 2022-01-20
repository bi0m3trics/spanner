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
#' @param output_location character Where to write intermediary and output data layers.
#' @param max_dai numeric The max diameter (in m) of a resulting tree (use to eliminate commission errors).
#' @param SDvert numeric The standard deviation threshold below whihc polygons will be considered as tree boles
#'
#' @return data.frame A data.frame containing the following seed information: `TreeID`,
#' `X`, `Y`, `Z`, and `Radius` in the same units as the .las
#'
#' @examples
#'
#'
#' \dontrun{
#' # set the number of threads to use in lidR
#' set_lidr_threads(8)
#'
#' # read the las (which must be downloaded with getExampleData())
#' LASfile <- system.file("extdata", "Pine_Example.laz", package="spanner")
#' las = readLAS(LASfile, select = "xyzc")
#'
#' # plot(las, color="Z", backend="lidRviewer", trim=30)
#'
#' # find tree locations and attribute data
#' myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05, pt_spacing = 0.0254,
#'                                        dens_threshold = 0.25,
#'                                        neigh_sizes = c(0.333, 0.166, 0.5),
#'                                        eigen_threshold = 0.5,
#'                                        grid_slice_min = 0.666,
#'                                        grid_slice_max = 2,
#'                                        minimum_polygon_area = 0.01,
#'                                        cylinder_fit_type = "ransac",
#'                                        output_location = getwd(),
#'                                        max_dia = 1,
#'                                        SDvert = 0.25)
# plot(lidR::grid_canopy(las, res = 0.2, p2r()))
# points(myTreeLocs$X, myTreeLocs$Y, col = "black", pch=16, cex = myTreeLocs$Radius^2*10)
#' }
#'
#' @export
get_raster_eigen_treelocs <- function(las = las, res = 0.05, pt_spacing = 0.0254, dens_threshold = 0.2,
                                      neigh_sizes=c(0.333, 0.166, 0.5), eigen_threshold = 0.6666,
                                      grid_slice_min = 0.6666, grid_slice_max = 2.0,
                                      minimum_polygon_area = 0.025, cylinder_fit_type = "ransac",
                                      output_location = getwd(), max_dia=0.5, SDvert = 0.25)
{

  ## Validate output directory
  if (!dir.exists(output_location))
  {
    stop("Output directory does not exist!")
  }
  if (!(file.access(output_location, mode = 0) == 0))
  {
    stop("Output directory not writeable!")
  }
  ##---------------------- Preprocesssing -------------------------------------
  ## Subsample using systematic voxel grid to 1in
  message("Downsampling the scan... (step 1 of 14)\n")
  las <- TreeLS::tlsSample(las, TreeLS::smp.voxelize(pt_spacing))

  ## Create the processing slice based on user's grid slice min/max
  message("Creating processing slice... (2/14)\n")
  slice_extra <- lidR::filter_poi(las, Z >= grid_slice_min, Z <= grid_slice_max)

  message("Calculating verticality... (3/14)\n")
  vert_temp <- spanner::eigen_metrics(slice_extra, radius=neigh_sizes[1], ncpu=lidR::get_lidr_threads())
  # vert_temp <- spanner::C_vert_in_sphere(slice_extra, radius = neigh_sizes[1], ncpu = lidR::get_lidr_threads())
  slice_extra <- lidR::add_lasattribute(slice_extra, as.numeric(vert_temp$Verticality),
                                        name = "verticality", desc = "verticality")
  slice_extra@data$verticality[is.na(slice_extra@data$verticality)] <- 0.5

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
  filename = paste0(output_location, paste0("temp_RelDens_",grid_slice_min,"-",grid_slice_max,"_",
                                                                 dens_threshold,"threshold_res-",res,".shp"))

  density_grid <- terra::rast(lidR::grid_metrics(slice_extra, ~median(relative_density, na.rm = T), res = res,
                                                 start = c(min(slice_extra@data$X), min(slice_extra@data$Y))))
  density_polygon <- terra::as.polygons(terra::classify(density_grid, rbind(c(0,dens_threshold,0),
                                                              c(dens_threshold,1,1))), dissolve = T)
  terra::writeVector(terra::disagg(density_polygon[density_polygon$V1 > 0,]), filename = filename, overwrite=TRUE)

  sf::sf_use_s2(FALSE) ## Turn off the s2 processing

  density_polygon <- sf::read_sf(filename)
  density_polygon$area = as.numeric(sf::st_area(density_polygon))
  density_polygon <- density_polygon[density_polygon$area > minimum_polygon_area,]
  density_polygon <- sfheaders::sf_remove_holes(density_polygon)
  density_polygon <- sf::st_buffer(sf::st_buffer(density_polygon, 1.5), -1.475) ## smoothing out the resulting polygons

  message("Gridding verticality... (7/14)\n")
  filename = paste0(output_location, paste0("temp_vertical_",grid_slice_min,"-",grid_slice_max,"_",
                                                                 eigen_threshold,"threshold_res-",res,".shp"))

  verticality_grid <- terra::rast(lidR::grid_metrics(slice_extra, ~quantile(verticality, 0.5, na.rm = T), res = res))
  verticality_polygon <- terra::as.polygons(terra::classify(verticality_grid, rbind(c(0,eigen_threshold,0),
                                                                      c(eigen_threshold,1,1))), dissolve = T)
  terra::writeVector(terra::disagg(verticality_polygon[verticality_polygon$V1 > 0,]), filename = filename, overwrite=TRUE)

  verticality_polygon <- sf::read_sf(filename)
  verticality_polygon$area = as.numeric(sf::st_area(verticality_polygon))
  verticality_polygon <- verticality_polygon[verticality_polygon$area > minimum_polygon_area,]
  verticality_polygon <- sfheaders::sf_remove_holes(verticality_polygon)
  verticality_polygon <- sf::st_buffer(sf::st_buffer(verticality_polygon, 1.5), -1.475) ## smoothing out the resulting polygons

  message("Merging relative density and verticality... (8/14)\n")
  merged <- sf::st_intersection(verticality_polygon,
                            density_polygon)
  merged$area <- as.numeric(sf::st_area(merged))
  merged <- merged[merged$area > minimum_polygon_area,]

  message("Filtering merged polygon... (9/14)\n")
  ## sd verticality
  verticalitySD_grid <- terra::rast(grid_metrics(slice_extra, ~sd(verticality), res = res))
  merged$SDvert <- round(terra::extract(verticalitySD_grid, terra::vect(sf::as_Spatial(merged)),
                                        fun = function(x){median(x, na.rm = T)})[,2],2)
  ## some simple filtering: lower canopy
  Z01 <- terra::rast(grid_metrics(slice_extra, res = res, ~quantile(Z, 0.01)))
  merged$lowHGT <- round(terra::extract(Z01, terra::vect(sf::as_Spatial(merged)),
                                        fun = function(x){quantile(x, 0.05, na.rm = T)})[,2],2)
  Z99 <- terra::rast( grid_metrics(slice_extra, res = res, ~quantile(Z, 0.99)) )
  merged$highHGT <- round(terra::extract(Z99, terra::vect(sf::as_Spatial(merged)),
                                         fun = function(x){quantile(x, 0.99, na.rm = T)})[,2],2)
  Z50 <- terra::rast(grid_metrics(slice_extra, res = res, ~quantile(Z, 0.5)))
  merged$medHGT <- round(terra::extract(Z50, terra::vect(sf::as_Spatial(merged)),
                                        fun = function(x){median(x, na.rm = T)})[,2],2)
  merged$diffHGT <- round(merged$highHGT - merged$lowHGT,2)
  # merged <- merged[merged$diffHGT > ((grid_slice_max - grid_slice_min) / 2), ]
  merged$maxHGT <- round(terra::extract(terra::rast(cancov), terra::vect(sf::as_Spatial(merged)),
                                        fun = max)[,2], 2)

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
  slice_clip <- lidR::merge_spatial(las = lidR::clip_roi(slice_extra, sf::st_sf(sf::st_union(circles_sf))),
                              source = circles_sf, attribute = "TreeID")

  slice_clip <- lidR::filter_poi(slice_clip, verticality >= 0.5)

  ##---------------------- Identify Trees -------------------------------------
  message("Fitting nested height (length) cylinders...(12/14)\n")
  cyl_fit <- list()
  for(t in 1:length(sort(unique(slice_clip$TreeID))))
    {
    min = grid_slice_min; max = grid_slice_max
    if(cylinder_fit_type == "ransac"){
      n_pts = 20
      n_best = 25
    } else if(cylinder_fit_type == "irls"){
      n_pts = min(data.frame(lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data %>%
                               dplyr::group_by(TreeID) %>%
                               dplyr::summarize(count = length(X)))$count)/2
      n_best = 100
    } else {
      print("You will need to choose either 'ransac' or 'irls' cylinder fitting method...")
      stop()
    }
    fit <- TreeLS:::cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t]),
                                method = cylinder_fit_type, n = n_pts, inliers = 0.9,
                                conf = 0.99, max_angle = 20, n_best = n_best)
    fit$TreeID <- sort(unique(slice_clip$TreeID))[t]
    fit$dbh_width <- max-min

    cyl_fit[[t]] <- fit

    min = 1.1 - ((1.1 - grid_slice_min)/2); max =  1.6 + ((grid_slice_max - 1.6)/2)
    if(cylinder_fit_type == "ransac"){
      n_pts = 20
      n_best = 25
    } else if(cylinder_fit_type == "irls"){
      n_pts = min(data.frame(lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data %>%
                               dplyr::group_by(TreeID) %>%
                               dplyr::summarize(count = length(X)))$count)/2
      n_best = 100
    } else {
      print("You will need to choose either 'ransac' or 'irls' cylinder fitting method...")
      stop()
    }
    fit <- TreeLS:::cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t],
                                           Z >= min, Z <= max),
                                method = cylinder_fit_type, n = n_pts, inliers = 0.9,
                                conf = 0.99, max_angle = 20, n_best = n_best)
    fit$TreeID <- sort(unique(slice_clip$TreeID))[t]
    fit$dbh_width <- max-min

    cyl_fit[[t+length(unique(slice_clip$TreeID))]] <- fit

    min = 1.1; max = 1.6
    if(cylinder_fit_type == "ransac"){
      n_pts = 20
      n_best = 25
    } else if(cylinder_fit_type == "irls"){
      n_pts = min(data.frame(lidR::filter_poi(slice_clip, Z <= max, Z >= min)@data %>%
                               dplyr::group_by(TreeID) %>%
                               dplyr::summarize(count = length(X)))$count)/2
      n_best = 100
    } else {
      print("You will need to choose either 'ransac' or 'irls' cylinder fitting method...")
      stop()
    }
    fit <- TreeLS:::cylinderFit(lidR::filter_poi(slice_clip, TreeID == sort(unique(slice_clip$TreeID))[t],
                                           Z >= min, Z <= max),
                                method = cylinder_fit_type, n = n_pts, inliers = 0.9,
                                conf = 0.99, max_angle = 20, n_best = n_best)
    fit$TreeID <- sort(unique(slice_clip$TreeID))[t]
    fit$dbh_width <- max-min

    cyl_fit[[t + ((length(unique(slice_clip$TreeID))) * 2)]] <- fit
  }
  (cyl_fit <- dplyr::bind_rows(cyl_fit))
  cyl_fit <- cyl_fit[cyl_fit$radius <= max_dia, ]
  cyl_fit <- cyl_fit[complete.cases(cyl_fit),]
  summary_cyl_fit <- cyl_fit %>% dplyr::group_by(TreeID) %>% dplyr::summarize_all(mean)

  ##---------------------- Cleaning up the Output -------------------------------------
  message("Successfully obtained the cylinder summaries... (13/14)\n")
  output<-summary_cyl_fit[,c("TreeID","px","py","pz","radius","err")]
  colnames(output)<-c("TreeID","X","Y","Z","Radius","Error")

  return(output)
  message("Done! (14/14)\n")
}
