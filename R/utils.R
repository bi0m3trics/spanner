.onUnload <- function(libpath)
{
  library.dynam.unload("spanner", libpath)
}


#' Calculates eigen decomposition metrics for fixed neighborhood point cloud data
#' @description This function calculates twelve (plus the first and second PCA) for
#' several point geometry-related metrics (listed below) in parallel using C++ for
#' a user-specified radius.
#' @param las LAS Normalized las object.
#' @param radius numeric the radius of the neighborhood
#' @param ncpu integer the number of cpu's to be used in parallelfor the calculation
#'
#' @section List of available point metrics:
#' \loadmathjax
#'  \itemize{
#'  \item \code{eLargest}: first eigenvalue, \mjeqn{\lambda_{1}}{ASCII representation}
#'  \item \code{eMedium}: second eigenvalue, \mjeqn{\lambda_{2}}{ASCII representation}
#'  \item \code{eSmallest}: third eigenvalue, \mjeqn{\lambda_{3}}{ASCII representation}
#'  \item \code{eSum}: sum of eigenvalues, \mjeqn{\sum_{i=1}^{n=3} \lambda_{i}}{ASCII representation}
#'  \item \code{Curvature}: surface variation, \mjeqn{\lambda_{3} / \sum_{i=1}^{n=3} \lambda_{i}}{ASCII representation}
#'  \item \code{Omnivariance}: high values correspond to spherical features and low values to planes or linear features, \mjeqn{(\lambda_{1} * \lambda_{2} * \lambda_{3})^{1/3}}{ASCII representation}
#'  \item \code{Anisotropy}: relationships between the directions of the point distribution, \mjeqn{(\lambda_{1} - \lambda_{3}) / \lambda_{1}}{ASCII representation}
#'  \item \code{Eigentropy}: entropy in the eigenvalues, \mjeqn{- \sum_{i=1}^{n=3} \lambda_{i} * ln(\lambda_{i})}{ASCII representation}
#'  \item \code{Linearity}: linear saliency, \mjeqn{(\lambda_{1} - \lambda_{2}) / \lambda_{1}}{ASCII representation}
#'  \item \code{Verticality}: vertical saliency, \mjeqn{1-abs(\langle (0,0,1),e_3\rangle)}{ASCII representation}
#'  \item \code{Planarity}: planar saliency, \mjeqn{(\lambda_{2} - \lambda_{3}) / \lambda_{1}}{ASCII representation}
#'  \item \code{Sphericity}: spherical saliency, \mjeqn{\lambda_{3} / \lambda_{1}}{ASCII representation}
#'  \item \code{Nx,Ny,Nz}: 3 components of the normal vector (smallest eigenvector)
#'  \item \code{SurfaceVariation}: surface variation (change of curvature), same as Curvature
#'  \item \code{ChangeCurvature}: alternative name for surface variation
#'  \item \code{SurfaceDensity}: 2D point density using circle area, \mjeqn{k / (\pi R^{2})}{ASCII representation}
#'  \item \code{VolumeDensity}: 3D point density using sphere volume, \mjeqn{k / (\frac{4}{3}\pi R^{3})}{ASCII representation}
#'  \item \code{MomentOrder1}: 1st order moment from CloudCompare, projection onto 2nd eigenvector, \mjeqn{m_{1}^{2} / m_{2}}{ASCII representation}
#'  \item \code{NormalChangeRate}: normal change rate, same as Curvature, \mjeqn{\lambda_{3} / \sum_{i=1}^{n=3} \lambda_{i}}{ASCII representation}
#'  \item \code{Roughness}: distance from query point to fitted plane, \mjeqn{|\vec{d} \cdot \vec{n}|}{ASCII representation}
#'  \item \code{MeanCurvature}: mean curvature from quadric surface fitting, \mjeqn{H = \frac{(1+f_{y}^{2})f_{xx} - 2f_{x}f_{y}f_{xy} + (1+f_{x}^{2})f_{yy}}{2(1+f_{x}^{2}+f_{y}^{2})^{3/2}}}{ASCII representation}
#'  \item \code{GaussianCurvature}: Gaussian curvature from quadric surface fitting, \mjeqn{K = \frac{f_{xx}f_{yy} - f_{xy}^{2}}{(1+f_{x}^{2}+f_{y}^{2})^{2}}}{ASCII representation}
#'  \item \code{PCA1}: eigenvector projection variance normalized by eigensum, \mjeqn{\sigma_{PC1}^{2} / \sum_{i=1}^{n=3} \lambda_{i}}{ASCII representation}
#'  \item \code{PCA2}: eigenvector projection variance normalized by eigensum, \mjeqn{\sigma_{PC2}^{2} / \sum_{i=1}^{n=3} \lambda_{i}}{ASCII representation}
#'  \item \code{NumNeighbors}: number of points in the spherical neighborhood, \mjeqn{k}{ASCII representation}
#' }
#' @return A labeled data.table of point metrics for each point in the LAS object
#'
#' @examples
#' \dontrun{
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(LASfile)
#' eigen = eigen_metrics(las, radius=2, ncpu=4)
#' }
#'
#' @export
eigen_metrics = function(las = las, radius=0.1, ncpu = 8){
  # Input validation
  if (!inherits(las, "LAS")) {
    stop("las must be a LAS object")
  }
  if (!is.numeric(radius) || length(radius) != 1) {
    stop("radius must be a single numeric value")
  }
  if (radius <= 0) {
    stop("radius must be positive")
  }
  if (!is.numeric(ncpu) || length(ncpu) != 1 || ncpu < 1) {
    stop("ncpu must be a positive integer")
  }
  
  temp = C_eigen_in_sphere(las, radius = radius, ncpu = ncpu)
  data.table::setDT(temp)
  cols<-c("eLargest","eMedium","eSmallest","eSum","Curvature","Omnivariance",
          "Anisotropy","Eigentropy","Linearity","Verticality","Planarity",
          "Sphericity", "Nx", "Ny", "Nz", "SurfaceVariation", "ChangeCurvature",
          "SurfaceDensity", "VolumeDensity", "MomentOrder1", "NormalChangeRate",
          "Roughness", "MeanCurvature", "GaussianCurvature", "PCA1", "PCA2", "NumNeighbors")
  data.table::setnames(temp, cols)
  return(temp)
}


#' Point cloud cylinder fitting as per de Conto et al. 2017 as implemented here: https://github.com/tiagodc/TreeLS
#' @description Fits a cylinder on a set of 3D points.
#' @param las LAS normalized and segmented las object.
#' @param method method for estimating the cylinder parameters. Currently available: \code{"nm"}, \code{"irls"}, \code{"ransac"} and \code{"bf"}.
#' @param n number of points selected on every RANSAC iteration.
#' @param inliers expected proportion of inliers among stem segments' point cloud chunks.
#' @param conf confidence level.
#' @param max_angle used when \code{method == "bf"}. The maximum tolerated deviation, in degrees, from an absolute vertical line (Z = c(0,0,1)).
#' @param n_best estimate optimal RANSAC parameters as the median of the \code{n_best} estimations with lowest error.
#' @return vector of parameters
#'
#' @examples
#' \dontrun{
#' # Define the cylinder attributes
#' npts = 500
#' cyl_length = 0.5
#' radius = 0.2718
#'
#' # Generate the X,Y,Z values
#' Z=runif(n = npts, min = 0, max = cyl_length)
#' angs = runif(npts, 0, 2*pi)
#' X = sin(angs)*radius
#' Y = cos(angs)*radius
#'
#' # Creation of a LAS object out of external data
#' cloud <- LAS(data.frame(X,Y,Z))
#'
#' # Fit a cylinder and retrun the information
#' cyl_par = spanner::cylinderFit(cloud, method = 'ransac', n=5, inliers=.9,
#'                                conf=.95, max_angle=30, n_best=20)
#' }
#' @export
cylinderFit = function(las, method = 'ransac', n=5, inliers=.9, conf=.95, max_angle=30, n_best=20){
  if(nrow(las@data) < 3) return(NULL)
  if(method == 'ransac' & nrow(las@data) <= n) method = 'nm'
  pars = cppCylinderFit(las2xyz(las), method, n, conf, inliers, max_angle, n_best)
  if(method == 'bf'){
    pars[3] = pars[3]
    names(pars) = c('x','y','radius', 'err', 'ax', 'ay')
  }else{
    pars[5] = pars[5]
    pars = c(pars, as.double(apply(las@data[,c('X','Y','Z')], 2, function(x) sum(range(x))/2)))
    names(pars) = c('rho','theta','phi', 'alpha', 'radius', 'err', 'px', 'py', 'pz')
  }
  pars = as.data.frame(t(pars))
  return(pars)
}

#' Convert LAS object to XYZ matrix
#' @description Extracts the X, Y, and Z coordinates from a LAS object and returns them as a matrix.
#' @param las LAS object to convert
#' @return A numeric matrix with three columns (X, Y, Z) containing the point coordinates
#'
#' @examples
#' \dontrun{
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(LASfile)
#' xyz_matrix <- las2xyz(las)
#' head(xyz_matrix)
#' }
#' @export
las2xyz = function(las){
  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")
  las = as.matrix(las@data[,c('X','Y','Z')])
  return(las)
}

#' Spanner color palette
#' @description Returns a named vector of colors for use in spanner visualizations.
#' The palette includes 10 distinct colors suitable for categorical data visualization.
#' @return A named character vector of hex color codes
#'
#' @examples
#' \dontrun{
#' # Get the palette
#' colors <- spanner_pal()
#'
#' # Use in a plot
#' barplot(1:10, col = spanner_pal(), names.arg = names(spanner_pal()), las = 2)
#' }
#' @export
spanner_pal <- function() {
  c(
    Red_CMYK = "#F52220",
    Tangerine = "#E8801A",
    Flax = "#E1DA8A",
    Kelly_green = "#5EA530",
    Brunswick_green = "#114232",
    Keppel = "#2BB4A2",
    Dark_slate_gray = "#115A5D",
    Magenta_dye = "#B20F66",
    Chocolate_cosmos = "#560F11",
    Cornell_red = "#C0181A"
  )
}
