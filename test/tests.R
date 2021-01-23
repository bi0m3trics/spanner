set_lidr_threads(10)

LASfile = system.file("extdata", "MixedConifer.laz", package="spanner")
las = readTLSLAS(LASfile, select = "xyz0", filter = "-drop_z_below 0")

temp1 = spanner::C_count_in_disc(las@data$X, las@data$Y, las@data$X, las@data$Y, radius = 3, ncpu = get_lidr_threads())
temp2 = spanner::C_count_in_sphere(las, radius = 3, ncpu = get_lidr_threads())
temp3 = spanner::C_vert_in_sphere(las, radius = 4, ncpu = get_lidr_threads())

las = add_attribute(las, temp3, name="verticality")
plot(las, color="verticality")

las = readTLSLAS("C:/Dropbox/GitHubRepos/Rudd_CO_SW_zebcam_ALSregistered.laz")
xcenter = ((lidR::extent(las)[2]-lidR::extent(las)[1])/2)+lidR::extent(las)[1]
ycenter = ((lidR::extent(las)[4]-lidR::extent(las)[3])/2)+lidR::extent(las)[3]
las = lidR::clip_rectangle(las, xleft=xcenter-15, ybottom=ycenter-15, xright=xcenter+15, ytop=ycenter+15)
las = lidR::classify_ground(las, lidR::csf(sloop_smooth = FALSE, class_threshold = 0.05, cloth_resolution = 0.05,
                                           rigidness = 1L, iterations = 500L, time_step = 0.65), last_returns = FALSE)
las = lidR::normalize_height(las, lidR::knnidw(k = 9, p = 2), use_class = 2L)
las = lidR::classify_noise(las, lidR::ivf(0.5,5))
las = lidR::filter_poi(las, Classification != LASNOISE)
# plot(las, color="Classification", backend="pcv")

t1 = Sys.time()
myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05, pt_spacing = 0.01, dens_threshold = 0.25, neigh_sizes=c(0.1, 0.1, 0.25),
                                       eigen_threshold = 0.25, grid_slice_min = 0.666, grid_slice_max = 2.0,
                                       minimum_polygon_area = 0.01, cylinder_fit_type = "ransac", 
                                       output_location = "C:/Dropbox/GitHubRepos/")
t2 = Sys.time()
t2-t1

# plot(lidR::grid_canopy(las, res = 0.5, pitfree(c(0,2,5,10,15), c(0, 1.5))))
# points(myTreeLocs$X, myTreeLocs$Y, col = "black", pch=16, cex = myTreeLocs$Radius^2*pi*10)


t1 = Sys.time()
myTreeGraph <- segment_graph(las = las, tree.locations = myTreeLocs, k = 50, distance.threshold = 0.38,
                             use.metabolic.scale = FALSE, subsample.graph = 0.2, return.dense = FALSE,
                             output_location = "C:/Dropbox/GitHubRepos/")
t2 = Sys.time()
t2-t1

# plot(myTreeGraph, color="treeID", backend="pcv")
