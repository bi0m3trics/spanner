% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Segment_Graph.R
\name{segment_graph}
\alias{segment_graph}
\title{Segment a terrestrial point cloud using graph theory.}
\usage{
segment_graph(
  las = las,
  tree.locations = NULL,
  k = NULL,
  distance.threshold = 0.38,
  use.metabolic.scale = FALSE,
  ptcloud_slice_min = 0.5,
  ptcloud_slice_max = 2,
  metabolic.scale.function = NULL,
  subsample.graph = 0.1,
  return.dense = FALSE,
  output_location = getwd()
)
}
\arguments{
\item{las}{LAS normalized las object.}

\item{tree.locations}{data.frame A data.frame contained the following seed information: \code{TreeID},
\code{X}, \code{Y}, \code{Z}, and \code{Radius} in the same units as the .las, output from the raster_eigen function.}

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

\item{output_location}{character Where to save processing outputs (rasters).}

\item{plot_radius}{numeric Radius (in the same units as the .las) used to define a region of
interest for processing}
}
\value{
a .las with the column treeID added.
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
# set the number of threads to use in lidR
set_lidr_threads(8)

# read the las (which must be downloaded with getExampleData())
LASfile <- system.file("extdata", "Pine_Example.laz", package="spanner")
las = readLAS(LASfile, select = "xyzc")

# plot(las, color="Z", backend="lidRviewer", trim=30)

# find tree locations and attribute data
myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05, pt_spacing = 0.02,
                                       dens_threshold = 0.25,
                                       neigh_sizes=c(0.333, 0.166, 0.5),
                                       eigen_threshold = 0.5,
                                       grid_slice_min = 0.666,
                                       grid_slice_max = 2,
                                       minimum_polygon_area = 0.01,
                                       cylinder_fit_type = "ransac",
                                       output_location = getwd(),
                                       max_dia = 0.5,
                                       SDvert = 0.25)

# segment the point cloud
myTreeGraph <- segment_graph(las = las, tree.locations = myTreeLocs, k = 50,
                             distance.threshold = 0.5,
                             use.metabolic.scale = FALSE, ptcloud_slice_min = 0.6666,
                             ptcloud_slice_max = 2.0, subsample.graph = 0.1,
                             return.dense = FALSE,
                             output_location = getwd())

plot(myTreeGraph, color="treeID",  backend="lidRviewer")
}

}
\references{
Tao, S., Wu, F., Guo, Q., Wang, Y., Li, W., Xue, B., ... & Fang, J. (2015). Segmenting tree
crowns from terrestrial and mobile LiDAR data by exploring ecological theories. ISPRS Journal of
Photogrammetry and Remote Sensing, 110, 66-76.
}
