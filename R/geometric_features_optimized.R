#' Test the simple geometric features function
#'
#' @param las LAS object from lidR package
#' @param radius Search radius for neighbor computation
#' @param max_neighbors Maximum number of neighbors to consider
#' @return List of geometric features
#' @export
test_geometric_features <- function(las, radius, max_neighbors = 50) {
    # Source the C++ function if not available
    if (!exists("C_geometric_features_simple", mode = "function")) {
        message("Loading C++ function...")
        sourceCpp(system.file("src", "geometric_features_simple.cpp", package = "spanner"))
    }

    # Call C++ function
    features <- C_geometric_features_simple(las, radius, max_neighbors)

    # Return as data frame for easier handling
    return(as.data.frame(features))
}
#'   computes all available features. Options include:
#'   - Eigenvalue features: "lambda1", "lambda2", "lambda3", "sum_eigenvalues"
#'   - Shape features: "omnivariance", "eigenentropy", "anisotropy", "planarity",
#'     "linearity", "sphericity", "verticality", "surface_variation"
#'   - Normal/PCA: "nx", "ny", "nz", "pca1", "pca2"
#'   - Curvature: "mean_curvature", "gaussian_curvature", "normal_change_rate"
#'   - Surface: "roughness", "signed_roughness"
#'   - Density: "num_neighbors", "surface_density", "volume_density"
#'   - Other: "first_order_moment", "height_above_ground", "relative_height", "intensity_variance"
#'
#' @return A data.table with all computed geometric features
#'
#' @section Feature Definitions:
#' \loadmathjax
#'
#' **Eigenvalue-based features:**
#' \itemize{
#'   \item \code{lambda1}: First eigenvalue, \mjeqn{\lambda_{1}}{ASCII representation}
#'   \item \code{lambda2}: Second eigenvalue, \mjeqn{\lambda_{2}}{ASCII representation}
#'   \item \code{lambda3}: Third eigenvalue, \mjeqn{\lambda_{3}}{ASCII representation}
#'   \item \code{sum_eigenvalues}: Sum of eigenvalues, \mjeqn{\sum_{i=1}^{3} \lambda_{i}}{ASCII representation}
#'   \item \code{omnivariance}: \mjeqn{(\lambda_{1} \lambda_{2} \lambda_{3})^{1/3}}{ASCII representation}
#'   \item \code{eigenentropy}: \mjeqn{-\sum_{i=1}^{3} \lambda_{i} \ln(\lambda_{i})}{ASCII representation}
#'   \item \code{anisotropy}: \mjeqn{(\lambda_{1} - \lambda_{3}) / \lambda_{1}}{ASCII representation}
#'   \item \code{planarity}: \mjeqn{(\lambda_{2} - \lambda_{3}) / \lambda_{1}}{ASCII representation}
#'   \item \code{linearity}: \mjeqn{(\lambda_{1} - \lambda_{2}) / \lambda_{1}}{ASCII representation}
#'   \item \code{sphericity}: \mjeqn{\lambda_{3} / \lambda_{1}}{ASCII representation}
#'   \item \code{verticality}: \mjeqn{1 - |\langle (0,0,1), e_3 \rangle|}{ASCII representation}
#'   \item \code{surface_variation}: \mjeqn{\lambda_{3} / \sum_{i=1}^{3} \lambda_{i}}{ASCII representation} (curvature)
#' }
#'
#' **Curvature features:**
#' \itemize{
#'   \item \code{mean_curvature}: Mean curvature from quadric surface fitting
#'   \item \code{gaussian_curvature}: Gaussian curvature from quadric surface fitting
#'   \item \code{normal_change_rate}: Same as surface_variation (CloudCompare definition)
#' }
#'
#' **Surface quality:**
#' \itemize{
#'   \item \code{roughness}: Distance from point to best-fit plane
#'   \item \code{signed_roughness}: Signed distance from point to best-fit plane
#' }
#'
#' **Density features:**
#' \itemize{
#'   \item \code{num_neighbors}: Number of neighbors within radius
#'   \item \code{surface_density}: \mjeqn{N / (\pi r^2)}{ASCII representation}
#'   \item \code{volume_density}: \mjeqn{N / (4/3 \pi r^3)}{ASCII representation}
#' }
#'
#' **Other features:**
#' \itemize{
#'   \item \code{first_order_moment}: Average distance to local centroid
#'   \item \code{height_above_ground}: Height above minimum Z in neighborhood
#'   \item \code{relative_height}: Normalized height in neighborhood (0-1)
#'   \item \code{intensity_variance}: Variance of intensity values in neighborhood (if available)
#'   \item \code{nx}, \code{ny}, \code{nz}: Components of surface normal vector
#'   \item \code{pca1}, \code{pca2}: First two principal component scores
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- readLAS(LASfile)
#'
#' # Compute all features (fastest)
#' features <- geometric_features_optimized(las, radius = 0.1, ncpu = 8)
#'
#' # Compute only specific features
#' shape_features <- geometric_features_optimized(las, radius = 0.1, ncpu = 8,
#'   features = c("linearity", "planarity", "sphericity", "roughness"))
#'
#' # Compare performance with original
#' system.time(old <- eigen_metrics(las, radius = 0.1, ncpu = 8))
#' system.time(new <- geometric_features_optimized(las, radius = 0.1, ncpu = 8))
#' }
#'
#' @export
geometric_features_optimized <- function(las, radius = 0.1, ncpu = 8, features = NULL) {

  # Validate inputs
  if (!inherits(las, "LAS")) {
    stop("Input must be a LAS object")
  }

  if (radius <= 0) {
    stop("Radius must be positive")
  }

  if (ncpu < 1) {
    ncpu <- 1
  }

  # Available feature names
  all_features <- c(
    "lambda1", "lambda2", "lambda3", "sum_eigenvalues", "omnivariance", "eigenentropy",
    "anisotropy", "planarity", "linearity", "sphericity", "verticality", "surface_variation",
    "pca1", "pca2", "nx", "ny", "nz", "mean_curvature", "gaussian_curvature",
    "normal_change_rate", "roughness", "signed_roughness", "num_neighbors",
    "surface_density", "volume_density", "first_order_moment", "height_above_ground",
    "relative_height", "intensity_variance"
  )

  # Validate feature selection
  if (!is.null(features)) {
    if (!all(features %in% all_features)) {
      invalid <- features[!features %in% all_features]
      stop(paste("Invalid feature names:", paste(invalid, collapse = ", ")))
    }
  }

  # Call optimized C++ function
  cat("Computing geometric features with nanoflann spatial indexing...\n")
  start_time <- Sys.time()

  # Use the correct C_geometric_features_simple function
  result <- C_geometric_features_simple(las, radius = radius, max_neighbors = 50, ncpu = ncpu)

  end_time <- Sys.time()
  cat(sprintf("Completed in %.2f seconds\n", as.numeric(end_time - start_time, units = "secs")))

  # Convert to data.table and keep original column names from C_geometric_features_simple
  result_dt <- data.table::as.data.table(result)

  # Filter to requested features if specified
  if (!is.null(features)) {
    available_features <- intersect(features, names(result_dt))
    if (length(available_features) > 0) {
      result_dt <- result_dt[, ..available_features]
    } else {
      warning("None of the requested features are available")
    }
  }

  return(result_dt)
}

#' Compute geometric features with CloudCompare compatibility
#'
#' @description Alias for geometric_features_optimized() that provides CloudCompare-style
#' feature computation with the exact same interface as CloudCompare's "Compute geometric features" tool.
#'
#' @inheritParams geometric_features_optimized
#' @return A data.table with CloudCompare-compatible feature names
#'
#' @examples
#' \dontrun{
#' # CloudCompare-style usage
#' features <- cloudcompare_features(las, radius = 0.1)
#'
#' # Specific features like CloudCompare interface
#' curvature_features <- cloudcompare_features(las, radius = 0.1,
#'   features = c("mean_curvature", "gaussian_curvature", "roughness"))
#' }
#'
#' @export
cloudcompare_features <- function(las, radius = 0.1, ncpu = 8, features = NULL) {
  return(geometric_features_optimized(las, radius, ncpu, features))
}
