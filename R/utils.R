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
#'  \item \code{Sphericity}: spherical saliency, \mjeqn{\lambda_{3} / \lambda_{1}}{ASCII representation}
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
          "Sphericity")
  data.table::setnames(temp, cols)
  pca_basic = princomp(temp[,c("eLargest","eMedium","eSmallest")])
  temp$PCA1 = pca_basic$scores[,1]
  temp$PCA2 = pca_basic$scores[,2]
  return(temp)
}
