library(spanner)

# set the number of threads to use in lidR
set_lidr_threads(8)

# download and read the laz
getExampleData()
LASfile = system.file("extdata", "Pine_Example.laz", package="spanner")
las = readLAS(LASfile, select = "xyzc")

# plot(las, color="Z", backend="lidRviewer", trim=30)

# find tree locations and attribute data
myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05, pt_spacing = 0.02, 
                                       dens_threshold = 0.25, 
                                       neigh_sizes=c(0.1, 0.1, 0.25), 
                                       eigen_threshold = 0.25, 
                                       grid_slice_min = 0.666, 
                                       grid_slice_max = 2,
                                       minimum_polygon_area = 0.01, 
                                       cylinder_fit_type = "ransac", 
                                       output_location = "E:/GithubRepos/", 
                                       max_dia = 0.5, 
                                       SDvert = 0.25)
# plot(lidR::grid_canopy(las, res = 0.2, p2r()))
# points(myTreeLocs$X, myTreeLocs$Y, col = "black", pch=16, cex = myTreeLocs$Radius^2*10, asp=1)

# segment the point cloud 
myTreeGraph = segment_graph(las = las, tree.locations = myTreeLocs, k = 50, 
                             distance.threshold = 0.5,
                             use.metabolic.scale = FALSE, subsample.graph = 0.1, 
                             return.dense = FALSE,
                             output_location = "E:/GithubRepos/")

# plot(myTreeGraph, color="treeID",  backend="lidRviewer")





