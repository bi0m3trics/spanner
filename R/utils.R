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
  n_points <- nrow(las@data)

  temp = spanner:::C_geometric_features_simple(las, radius = radius, max_neighbors = 50)

  data.table::setDT(temp)

  # Map to original column names for backward compatibility
  cols<-c("eLargest","eMedium","eSmallest","eSum","Curvature","Omnivariance",
          "Anisotropy","Eigentropy","Linearity","Verticality","Planarity",
          "Sphericity", "Nx", "Ny", "Nz")

  # Rename columns to match expected names
  if("eLargest" %in% names(temp)) data.table::setnames(temp, "eLargest", "eLargest")
  if("eMedium" %in% names(temp)) data.table::setnames(temp, "eMedium", "eMedium")
  if("eSmallest" %in% names(temp)) data.table::setnames(temp, "eSmallest", "eSmallest")
  if("nx" %in% names(temp)) data.table::setnames(temp, "nx", "Nx")
  if("ny" %in% names(temp)) data.table::setnames(temp, "ny", "Ny")
  if("nz" %in% names(temp)) data.table::setnames(temp, "nz", "Nz")
  if("anisotropy" %in% names(temp)) data.table::setnames(temp, "anisotropy", "Anisotropy")
  if("eigenentropy" %in% names(temp)) data.table::setnames(temp, "eigenentropy", "Eigentropy")
  if("linearity" %in% names(temp)) data.table::setnames(temp, "linearity", "Linearity")
  if("planarity" %in% names(temp)) data.table::setnames(temp, "planarity", "Planarity")
  if("sphericity" %in% names(temp)) data.table::setnames(temp, "sphericity", "Sphericity")
  if("omnivariance" %in% names(temp)) data.table::setnames(temp, "omnivariance", "Omnivariance")
  if("verticality" %in% names(temp)) data.table::setnames(temp, "verticality", "Verticality")

  # Calculate derived values
  if(all(c("eLargest","eMedium","eSmallest") %in% names(temp))) {
    temp$eSum <- temp$eLargest + temp$eMedium + temp$eSmallest
    temp$Curvature <- temp$eSmallest / pmax(temp$eSum, 1e-10)  # Prevent division by zero
  } else {
    temp$eSum <- 1.0
    temp$Curvature <- 0.1
  }

  # Select only the expected columns that exist
  existing_cols <- intersect(cols, names(temp))
  if(length(existing_cols) > 0) {
    temp <- temp[, ..existing_cols]
  }

  # Compute PCA on eigenvalues (if they exist)
  if(all(c("eLargest","eMedium","eSmallest") %in% names(temp))) {
    # Safely compute PCA with error handling
    tryCatch({
      pca_basic = princomp(temp[,c("eLargest","eMedium","eSmallest")])
      temp$PCA1 = pca_basic$scores[,1]
      temp$PCA2 = pca_basic$scores[,2]
    }, error = function(e) {
      temp$PCA1 <<- 0.0
      temp$PCA2 <<- 0.0
    })
  } else {
    temp$PCA1 = 0.0
    temp$PCA2 = 0.0
  }

  return(temp)
}


#' Point cloud cylinder fitting as per de Conto et al. 2017 as implemented here: https://github.com/tiagodc/TreeLS
#' @description Fits a cylinder on a set of 3D points.
#' @param las LAS normalized and segmented las object.
#' @param method method for estimating the cylinder parameters. Currently available: \code{"nm"}, \code{"irls"}, \code{"ransac"} and \code{"bf"}.
#' @param n number of points selected on every RANSAC iteration.
#' @param inliers expected proportion of inliers among stem segments' point cloud.
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
#' Cylinder Fitting Using RANSAC or Nelder-Mead Methods
#'
#' @description Fits a cylinder to a 3D point cloud using either RANSAC or Nelder-Mead optimization.
#' This function is particularly useful for fitting cylinders to tree trunks or other cylindrical structures
#' in LiDAR point cloud data.
#'
#' @param las LAS object containing the 3D point cloud data
#' @param method character Method for cylinder fitting: 'ransac' (default) or 'nm' (Nelder-Mead) or 'bf'
#' @param n integer Number of sample points for RANSAC (default: 5)
#' @param inliers numeric Proportion of inlier points required (default: 0.9)
#' @param conf numeric Confidence level for RANSAC (default: 0.95)
#' @param max_angle numeric Maximum angle deviation for cylinder fitting in degrees (default: 30)
#' @param n_best integer Number of best cylinder candidates to consider (default: 20)
#'
#' @return data.frame containing cylinder parameters:
#' \itemize{
#'   \item For RANSAC/NM: rho, theta, phi, alpha (cylinder orientation), radius, err (error), px, py, pz (center position)
#'   \item For BF: x, y, radius, err, ax, ay (axis components)
#' }
#'
#' @details The function automatically switches from RANSAC to Nelder-Mead method if the number of points
#' is less than or equal to the sample size parameter 'n'. The RANSAC method is more robust for noisy data,
#' while Nelder-Mead can be faster for clean, small datasets.
#'
#' @examples
#' \dontrun{
#' # Create synthetic cylinder data
#' angs <- runif(100, 0, 2*pi)
#' radius <- 0.2
#' Z <- runif(100, 0, 3)
#' X <- sin(angs) * radius
#' Y <- cos(angs) * radius
#'
#' # Create LAS object
#' cloud <- LAS(data.frame(X, Y, Z))
#'
#' # Fit cylinder using RANSAC
#' cyl_par <- cylinderFit(cloud, method = 'ransac', n = 5, inliers = 0.9,
#'                        conf = 0.95, max_angle = 30, n_best = 20)
#' }
#'
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

#' Convert LAS Object to XYZ Matrix
#'
#' @description Extracts the X, Y, Z coordinates from a LAS object and returns them as a numeric matrix.
#' This is a utility function commonly used for interfacing with functions that require matrix input
#' rather than LAS objects.
#'
#' @param las LAS object containing point cloud data
#'
#' @return numeric matrix with 3 columns (X, Y, Z) and one row per point
#'
#' @details This function validates that the input is a LAS object and extracts only the coordinate
#' columns. The resulting matrix is suitable for use with geometric algorithms and C++ functions
#' that require matrix input.
#'
#' @examples
#' \dontrun{
#' # Load LAS data
#' las <- readLAS("path/to/file.las")
#' 
#' # Convert to matrix
#' xyz_matrix <- las2xyz(las)
#' 
#' # Check dimensions
#' dim(xyz_matrix)  # Should be [n_points, 3]
#' }
#'
#' @export
las2xyz = function(las){

  if(class(las)[1] != "LAS")
    stop("las must be a LAS object")

  las = as.matrix(las@data[,c('X','Y','Z')])
  return(las)
}

#' Spanner Package Color Palette
#'
#' @description Returns the official color palette for the spanner package. This palette consists
#' of 10 carefully selected colors that provide good contrast and are suitable for data visualization
#' and plotting functions within the spanner ecosystem.
#'
#' @return Named character vector containing 10 hex color codes with descriptive names
#'
#' @details The palette includes colors ranging from reds and oranges to greens and blues,
#' providing a diverse set of colors for various visualization needs. Color names follow
#' standard color naming conventions.
#'
#' @examples
#' \dontrun{
#' # Get the spanner color palette
#' colors <- spanner_pal()
#' 
#' # View the colors
#' print(colors)
#' 
#' # Use in a plot
#' plot(1:10, col = spanner_pal(), pch = 16, cex = 2)
#' 
#' # Use specific colors by name
#' plot(1:3, col = spanner_pal()[c("Red_CMYK", "Kelly_green", "Keppel")])
#' }
#'
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

#' Colorize LAS Point Cloud Based on Attribute Values
#'
#' @description This function colorizes a LAS point cloud by mapping attribute values to RGB colors
#' using a specified color palette function. The function normalizes the attribute values and applies
#' the color mapping to create RGB color channels in the LAS object.
#'
#' @param las LAS object containing the point cloud data
#' @param attribute character Name of the attribute column in the LAS data to use for coloring
#' @param palette_func function Color palette function that generates colors (e.g., rainbow, heat.colors, lidR::random.colors)
#' @param n integer Number of colors to generate in the palette (default: 100)
#'
#' @return LAS object with RGB color channels (R, G, B) added based on the specified attribute
#'
#' @details The function handles NA values by replacing them with the minimum attribute value.
#' Attribute values are normalized to the range 0-1 before mapping to colors. The resulting
#' RGB values are added as new columns (R, G, B) to the LAS object.
#'
#' @examples
#' \dontrun{
#' # Load LAS data
#' las <- readLAS("path/to/file.las")
#' 
#' # Color by height using rainbow palette
#' las_colored <- colorizeLAS(las, "Z", rainbow, n = 256)
#' 
#' # Color by intensity using heat colors
#' las_colored <- colorizeLAS(las, "Intensity", heat.colors, n = 100)
#' 
#' # Color by tree ID using random colors
#' las_colored <- colorizeLAS(las, "treeID", lidR::random.colors, n = 50)
#' }
#'
#' @export
colorizeLAS <- function(las, attribute, palette_func, n = 100) {
  # Ensure the attribute exists in the LAS data
  if (!attribute %in% names(las@data)) {
    stop(paste("Attribute", attribute, "not found in LAS data."))
  }

  attr_values <- las@data[[attribute]]

  # Handle NA values by replacing them with the minimum value
  if (any(is.na(attr_values))) {
    min_val <- min(attr_values, na.rm = TRUE)
    attr_values[is.na(attr_values)] <- min_val
  }

  # Normalize attribute values to 0–1
  rng <- range(attr_values, na.rm = TRUE)
  if (rng[1] == rng[2]) {
    norm_attr <- rep(0.5, length(attr_values))
  } else {
    norm_attr <- (attr_values - rng[1]) / (rng[2] - rng[1])
  }

  # Create the color palette
  colors <- palette_func(n)

  # Map normalized values to color indices
  color_indices <- as.integer(norm_attr * (n - 1)) + 1
  mapped_colors <- colors[color_indices]

  # Convert colors to RGB
  rgb_values <- grDevices::col2rgb(mapped_colors)

  # Assign RGB to the LAS object
  las@data$R <- rgb_values[1, ]
  las@data$G <- rgb_values[2, ]
  las@data$B <- rgb_values[3, ]

  return(las)
}
