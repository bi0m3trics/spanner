% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Segment_Graph.R
\name{segment_graph}
\alias{segment_graph}
\title{Segment a terrestrial point cloud using graph theory.}
\usage{
segment_graph(
  las,
  tree.locations,
  k = 50,
  distance.threshold = 0.38,
  use.metabolic.scale = FALSE,
  ptcloud_slice_min = 0.5,
  ptcloud_slice_max = 2,
  metabolic.scale.function = NULL,
  subsample.graph = 0.1,
  return.dense = FALSE
)
}
\arguments{
\item{las}{LAS normalized las object.}

\item{tree.locations}{sf object \code{sf} object containing the following tree information: \code{TreeID},
\code{X}, \code{Y}, \code{Z}, \code{Radius}, and \code{Error}, output from the get_raster_eigen_treelocs function.}

\item{k}{integer Number of nearest neighbors to be used in processing (k >= 50 suggested)}

\item{distance.threshold}{numeric Maximum distance (in the same units as the .las) under which
two points are connected in the graph object (greater than point spacing). Two points with a greater
distance than this threshold are not connected in the graph for processing.}

\item{use.metabolic.scale}{bool Use of weights in the assignment of points to a given treeID. Useful
when interlocking crowns are present and trees are of different sizes.}

\item{ptcloud_slice_min}{numeric Lower bound of point cloud slice in normalized point cloud used for treeID matching.}

\item{ptcloud_slice_max}{numeric Upper bound of point cloud slice in normalized point cloud used for treeID matching.}

\item{metabolic.scale.function}{string Supply your own function for defining segmentation weights
based on a function of estimated tree diameter (e.g. metabolic.scale.function = 'x/2'). use.metabolic.scale
must be set to TRUE. If not supplied, defaults to metabolic scale function from Tao et al., 2015.}

\item{subsample.graph}{numeric The subsampled point spacing to use during processing. Note: processing
time increases quickly with smaller point spacing with negligible returns in accuracy.}

\item{return.dense}{bool Decision to return the subsampled point cloud or assign treeIDs back to
points in the input dense point cloud.}

\item{plot_radius}{numeric Radius (in the same units as the .las) used to define a region of
interest for processing}
}
\value{
a sparse/dense normalized .las with the column treeID added.
}
\description{
\code{segment_graph} returns a .las object with a new column "treeID".
}
\details{
Preforms Individual tree segmentation following ecological principles for “growing” trees
based on these input locations in a graph-theory approach inspired by work of Tao and
others (2015). Point coordinates are linked together based on proximity and turned into
a connected graph object, using the estimated tree bole locations as origin points, connecting individual
points back to those tree bole origins based on shortest paths within the graph network, and finally
assigning those points a unique tree identification based on the bole coordinate for which
they are connected. Input point cloud is subsampled to a lower resolution before processing to
increase processing efficiency. However, graph objects can still get large quite rapidly. Take
this into consideration when choosing the extent of the input las object.
}
\examples{
\dontrun{
\dontrun{
# Set the number of threads to use in lidR
set_lidr_threads(8)

# Download and read an example laz
getExampleData("DensePatchA")
LASfile = system.file("extdata", "DensePatchA.laz", package="spanner")
las = readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.01")
# Don't forget to make sure the las object has a projection
projection(las) = sp::CRS("+init=epsg:26912")

# Pre-process the example lidar dataset by classifying the ground points
# using lidR::csf(), normalizing it, and removing outlier points
# using lidR::ivf()
las = classify_ground(las, csf(sloop_smooth = FALSE,
                                class_threshold = 0.5,
                                cloth_resolution = 0.5, rigidness = 1L,
                                iterations = 500L, time_step = 0.65))
las = normalize_height(las, tin())
las = classify_noise(las, ivf(0.25, 3))
las = filter_poi(las, Classification != LASNOISE)

# Plot the non-ground points, colored by height
plot(filter_poi(las, Classification != 2), color = "Z")

# Perform a deep inspection of the las object. If you see any
# red text, you may have issues!
las_check(las)

# Find individual tree locations and attribute data
myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05,
                                       pt_spacing = 0.0254,
                                       dens_threshold = 0.2,
                                       neigh_sizes = c(0.333, 0.166, 0.5),
                                       eigen_threshold = 0.5,
                                       grid_slice_min = 0.6666,
                                       grid_slice_max = 2.0,
                                       minimum_polygon_area = 0.025,
                                       cylinder_fit_type = "ransac",
                                       max_dia = 0.5,
                                       SDvert = 0.25,
                                       n_pts = 20,
                                       n_best = 25)

# Plot the tree information over a CHM
plot(lidR::grid_canopy(las, res = 0.2, p2r()))
points(myTreeLocs$X, myTreeLocs$Y, col = "black", pch = 16,
       cex = myTreeLocs$Radius^2 * 10, asp = 1)

# Segment the point cloud
myTreeGraph = segment_graph(las = las, tree.locations = myTreeLocs, k = 50,
                             distance.threshold = 0.5,
                             use.metabolic.scale = FALSE,
                             ptcloud_slice_min = 0.6666,
                             ptcloud_slice_max = 2.0,
                             subsample.graph = 0.1,
                             return.dense = FALSE)

# Plot it in 3D colored by treeID
plot(myTreeGraph, color = "treeID", pal=spanner_pal())

}
}

}
\references{
Tao, S., Wu, F., Guo, Q., Wang, Y., Li, W., Xue, B., ... & Fang, J. (2015). Segmenting tree
crowns from terrestrial and mobile LiDAR data by exploring ecological theories. ISPRS Journal of
Photogrammetry and Remote Sensing, 110, 66-76.
}
