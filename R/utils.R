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
#'  \item \code{Linearity}: linear saliency, \mjeqn{(\lambda_{1} + \lambda_{2}) / \lambda_{1}}{ASCII representation}
#'  \item \code{Verticality}: vertical saliency, \mjeqn{1-abs(\langle (0,0,1),e_3\rangle)}{ASCII representation}
#'  \item \code{Planarity}: planar saliency, \mjeqn{(\lambda_{2} + \lambda_{3}) / \lambda_{1}}{ASCII representation}
#'  \item \code{Nx,Ny,Nz}: 3 components of the normal vector, {ASCII representation}
#' }
#' @return A labeled data.table of point metrics for each point in the LAS object
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(LASfile)
#' eigen = eigen_metrics(las, radius=2, ncpu=4)
#'
#' @export
eigen_metrics = function(las = las, radius=0.1, ncpu = 8){
  temp = spanner:::C_eigen_in_sphere(las, radius = radius, ncpu = ncpu)
  data.table::setDT(temp)
  cols<-c("eLargest","eMedium","eSmallest","eSum","Curvature","Omnivariance",
          "Anisotropy","Eigentropy","Linearity","Verticality","Planarity",
          "Sphericity", "Nx", "Ny", "Nz")
  data.table::setnames(temp, cols)
  pca_basic = princomp(temp[,c("eLargest","eMedium","eSmallest")])
  temp$PCA1 = pca_basic$scores[,1]
  temp$PCA2 = pca_basic$scores[,2]
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
#' @examples
#' # Define the cylinder attributes
#' npts = 500
#' cyl_length = 0.5
#' radius = 0.2718
#'
#' # Generate the X,Y,Z values
#' Z=runif(n = n, min = 0, max = cyl_length)
#' angs = runif(n, 0, 2*pi)
#' X = sin(angs)*radius
#' Y = cos(angs)*radius
#'
#' # Creation of a LAS object out of external data
#' cloud <- LAS(data.frame(X,Y,Z))
#'
#' # Fit a cylinder and retrun the information
#' cyl_par = spanner::cylinderFit(cloud, method = 'ransac', n=5, inliers=.9,
#'                                conf=.95, max_angle=30, n_best=20)
#' @export
cylinderFit = function(las, method = 'ransac', n=5, inliers=.9, conf=.95, max_angle=30, n_best=20){
  if(nrow(las@data) < 3) return(NULL)
  if(method == 'ransac' & nrow(las@data) <= n) method = 'nm'
  pars = cppCylinderFit(las %>% las2xyz, method, n, conf, inliers, max_angle, n_best)
  if(method == 'bf'){
    pars[3] = pars[3]
    names(pars) = c('x','y','radius', 'err', 'ax', 'ay')
  }else{
    pars[5] = pars[5]
    pars %<>% c(apply(las@data[,.(X,Y,Z)], 2, function(x) sum(range(x))/2) %>% as.double)
    names(pars) = c('rho','theta','phi', 'alpha', 'radius', 'err', 'px', 'py', 'pz')
  }
  pars = pars %>% t %>% as.data.frame
  return(pars)
}

#' @export
las2xyz = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  las = as.matrix(las@data[,c('X','Y','Z')])
  return(las)
}

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
