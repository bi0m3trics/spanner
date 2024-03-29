% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Raster_Eigen_TreeLocations.R
\name{get_raster_eigen_treelocs}
\alias{get_raster_eigen_treelocs}
\title{Obtain tree information by rasterizing point cloud values of relative neighborhood density and
verticality within a slice of a normalized point cloud}
\usage{
get_raster_eigen_treelocs(
  las = las,
  res = 0.05,
  pt_spacing = 0.0254,
  dens_threshold = 0.2,
  neigh_sizes = c(0.333, 0.166, 0.5),
  eigen_threshold = 0.6666,
  grid_slice_min = 0.6666,
  grid_slice_max = 2,
  minimum_polygon_area = 0.025,
  cylinder_fit_type = "ransac",
  output_location = getwd(),
  max_dia = 0.5,
  SDvert = 0.25
)
}
\arguments{
\item{las}{LAS Normalized las object.}

\item{res}{numeric Pixel width of rasterized point cloud metrics.}

\item{pt_spacing}{numeric Subsample spacing for graph connections.}

\item{dens_threshold}{numeric Minimum point density in raster cell to be considered as potential tree bole.}

\item{eigen_threshold}{numeric Minimum average verticality in raster cell to be considered as potential tree bole.}

\item{grid_slice_min}{numeric Lower bound of point cloud slice in normalized point cloud.}

\item{grid_slice_max}{numeric Upper bound of point cloud slice in normalized point cloud.}

\item{minimum_polygon_area}{numeric Smallest allowable polygon area of potential tree boles.}

\item{cylinder_fit_type}{character Choose "ransac" or "irls" cylinder fitting.}

\item{output_location}{character Where to write intermediary and output data layers.}

\item{SDvert}{numeric The standard deviation threshold below whihc polygons will be considered as tree boles}

\item{neigh_size}{numeric Vector for verticality and relative density (small and large neighborhoods) calculations}

\item{max_dai}{numeric The max diameter (in m) of a resulting tree (use to eliminate commission errors).}
}
\value{
data.frame A data.frame containing the following seed information: \code{TreeID},
\code{X}, \code{Y}, \code{Z}, and \code{Radius} in the same units as the .las
}
\description{
\code{get_raster_eigen_treelocs} returns a data.frame containing \code{TreeID}, \code{X}, \code{Y}, \code{Z}, \code{Radius}
and \code{Error} in the same units as the .las
}
\details{
For terrestrial and mobile lidar datasets, tree locations and estimates of DBH are provided
by rasterizing individual point cloud values of relative neighborhood density (at 0.3 and
1 m radius) and verticality within a slice of the normalized point cloud around breast height
(1.34 m). The algorithim then uses defined threshold values to classify the resulting rasters
and create unique polygons from the resulting classified raster. These point-density and
verticality polygons were selected by their intersection with one another, resulting in a
final set of polygons which were used to clip out regions of the point cloud that were most
likely to represent tree boles. A RANSAC cylinder fitting algorithm was then used to estimate
the fit of a cylinder to individual bole points. Cylinder centers and radius were used as inputs
to an individual tree segmentation
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
myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05, pt_spacing = 0.0254,
                                       dens_threshold = 0.25,
                                       neigh_sizes = c(0.333, 0.166, 0.5),
                                       eigen_threshold = 0.5,
                                       grid_slice_min = 0.666,
                                       grid_slice_max = 2,
                                       minimum_polygon_area = 0.01,
                                       cylinder_fit_type = "ransac",
                                       output_location = getwd(),
                                       max_dia = 1,
                                       SDvert = 0.25)
}

}
