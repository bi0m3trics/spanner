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

#' Colorize a LAS object using multiple methods
#' @description Colors a LAS object using one of three methods: attribute-based coloring,
#' raster-based RGB coloring, or PCV (Portion de Ciel Visible) ambient occlusion.
#'
#' @param las LAS object to colorize
#' @param method Character string specifying the coloring method: "attr" for attribute-based,
#'   "rgb" for raster-based, or "pcv" for ambient occlusion. Default is "attr".
#' @param attribute_name Character string specifying the attribute name (required for method="attr").
#'   The attribute must exist in the LAS data.
#' @param palette Character vector of at least two colors for the color ramp (used with method="attr", "pcv", or "ssao").
#'   Colors can be hex codes (e.g., "#FF0000") or named colors (e.g., "red"). Default is grayscale.
#' @param raster_path Character string or vector of paths to raster files (required for method="rgb").
#'   Can be a single RGB raster or three separate rasters for R, G, and B channels.
#' @param radius Numeric radius for neighborhood search in PCV calculation (method="pcv"). Default is 1.0.
#' @param num_directions Integer number of directional rays for PCV calculation (method="pcv"). Default is 60.
#' @param kernel_size Integer kernel size in pixels for SSAO sampling (method="ssao"). Default is 5.
#' @param pixel_size Numeric resolution of the depth map in spatial units (method="ssao"). Default is 0.1.
#' @param num_samples Integer number of samples per point for SSAO (method="ssao"). Default is 16.
#' @param ncpu Integer number of CPUs to use for parallel processing (method="pcv" or "ssao"). Default is 4.
#'
#' @return A LAS object with updated R, G, and B fields based on the selected method.
#'
#' @details
#' The function supports four coloring methods:
#' \describe{
#'   \item{attr}{Attribute-based coloring: normalizes attribute values and maps them to colors using the palette.}
#'   \item{rgb}{Raster-based coloring: extracts RGB values from georeferenced raster(s) that align with the point cloud.
#'     Requires matching CRS between LAS and raster. Can use a single 3-band RGB raster or three separate rasters.}
#'   \item{pcv}{PCV (Portion de Ciel Visible): computes 3D ambient occlusion by calculating sky visibility for each point.
#'     Based on the algorithm from Duguet & Girardeau-Montaut (2004). More accurate but slower than SSAO.}
#'   \item{ssao}{SSAO (Screen Space Ambient Occlusion): fast ambient occlusion using 2D depth map techniques.
#'     Projects points to a depth buffer and calculates occlusion based on depth differences. Much faster than PCV.}
#' }
#'
#' @examples
#' \donttest{
#'
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package="spanner")
#' las <- readLAS(LASfile, select = "xyz")
#'
#' # Attribute-based coloring
#' las_colored <- colorize_las(las, method="attr", attribute_name="Z",
#'                             palette=c("blue", "green", "yellow", "red"))
#'
#' # Raster-based coloring with RGB file
#' rgb_file <- system.file("extdata", "UAS_Clip_RGB.tif", package="spanner")
#' las_colored <- colorize_las(las, method="rgb", raster_path=rgb_file)
#'
#' # PCV ambient occlusion (slow, high quality)
#' las_colored <- colorize_las(las, method="pcv", radius=1.0,
#'                             num_directions=30, palette=c("black", "white"))
#'
#' # SSAO ambient occlusion (faster alternative to PCV)
#' las_colored <- colorize_las(las, method="ssao", pixel_size=0.1,
#'                             kernel_size=5, num_samples=16, palette=c("black", "white"), ncpu=8)
#' }
#'
#' @export
#'
#' @export
colorize_las <- function(las, method = "attr", attribute_name = NULL, palette = c("black", "white"),
                        raster_path = NULL, radius = 1.0, num_directions = 60,
                        kernel_size = 5, pixel_size = 0.1, num_samples = 16, ncpu = 4) {

  # Check if the input is a valid LAS object
  if (!inherits(las, "LAS")) {
    stop("The input is not a valid LAS object.")
  }

  # Validate method
  method <- tolower(method)
  if (!method %in% c("attr", "rgb", "pcv", "ssao")) {
    stop("Method must be one of 'attr', 'rgb', 'pcv', or 'ssao'.")
  }

  # Method-specific processing
  if (method == "attr") {
    # Attribute-based coloring (original implementation)
    if (is.null(attribute_name)) {
      stop("attribute_name must be specified when method='attr'.")
    }

    if (!attribute_name %in% names(las@data)) {
      stop(paste("Attribute", attribute_name, "not found in the LAS object."))
    }

    if (!is.character(palette) || length(palette) < 2) {
      stop("The palette must be a character vector with at least two colors.")
    }

    # Extract and normalize attribute values
    attribute_values <- las@data[[attribute_name]]
    normalized_values <- (attribute_values - min(attribute_values, na.rm = TRUE)) /
      (max(attribute_values, na.rm = TRUE) - min(attribute_values, na.rm = TRUE))

    # Map to colors
    colors <- grDevices::colorRampPalette(palette)(256)
    rgb_colors <- colors[as.numeric(cut(normalized_values, breaks = 256))]
    rgb_matrix <- grDevices::col2rgb(rgb_colors)

    las@data$R <- rgb_matrix[1, ]
    las@data$G <- rgb_matrix[2, ]
    las@data$B <- rgb_matrix[3, ]

  } else if (method == "rgb") {
    # Raster-based coloring
    if (is.null(raster_path)) {
      stop("raster_path must be specified when method='rgb'.")
    }

    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required for raster-based coloring. Please install it.")
    }

    # Load raster(s)
    if (length(raster_path) == 1) {
      # Single RGB raster
      raster <- terra::rast(raster_path)
      if (terra::nlyr(raster) < 3) {
        stop("Raster must have at least 3 bands for RGB coloring.")
      }
      r_band <- raster[[1]]
      g_band <- raster[[2]]
      b_band <- raster[[3]]
    } else if (length(raster_path) == 3) {
      # Separate R, G, B rasters
      r_band <- terra::rast(raster_path[1])
      g_band <- terra::rast(raster_path[2])
      b_band <- terra::rast(raster_path[3])
    } else {
      stop("raster_path must be either a single file or a vector of 3 files (R, G, B).")
    }

    # Check CRS compatibility
    las_crs <- sf::st_crs(las)
    raster_crs <- terra::crs(r_band, describe = TRUE)$code

    if (!is.na(las_crs$epsg) && !is.null(raster_crs)) {
      if (las_crs$epsg != as.numeric(raster_crs)) {
        warning("CRS mismatch between LAS and raster. Results may be incorrect.")
      }
    }

    # Extract RGB values for each point
    coords <- las@data[, c("X", "Y")]

    r_values <- terra::extract(r_band, coords, method = "simple")[, 2]
    g_values <- terra::extract(g_band, coords, method = "simple")[, 2]
    b_values <- terra::extract(b_band, coords, method = "simple")[, 2]

    # Handle NA values (points outside raster extent)
    r_values[is.na(r_values)] <- 0
    g_values[is.na(g_values)] <- 0
    b_values[is.na(b_values)] <- 0

    # Normalize to 0-255 if needed
    if (max(r_values, na.rm = TRUE) <= 1) {
      r_values <- r_values * 255
      g_values <- g_values * 255
      b_values <- b_values * 255
    }

    las@data$R <- as.integer(r_values)
    las@data$G <- as.integer(g_values)
    las@data$B <- as.integer(b_values)

  } else if (method == "pcv") {
    # PCV (Portion de Ciel Visible) ambient occlusion
    if (!is.character(palette) || length(palette) < 2) {
      stop("The palette must be a character vector with at least two colors.")
    }

    # Compute PCV values
    message("Computing PCV ambient occlusion...")
    pcv_values <- compute_pcv(las, radius = radius, num_directions = num_directions, ncpu = ncpu)

    # Normalize PCV values
    normalized_values <- (pcv_values - min(pcv_values, na.rm = TRUE)) /
      (max(pcv_values, na.rm = TRUE) - min(pcv_values, na.rm = TRUE))

    # Map to colors
    colors <- grDevices::colorRampPalette(palette)(256)
    rgb_colors <- colors[as.numeric(cut(normalized_values, breaks = 256))]
    rgb_matrix <- grDevices::col2rgb(rgb_colors)

    las@data$R <- rgb_matrix[1, ]
    las@data$G <- rgb_matrix[2, ]
    las@data$B <- rgb_matrix[3, ]

  } else if (method == "ssao") {
    # SSAO (Screen Space Ambient Occlusion) - fast ambient occlusion
    if (!is.character(palette) || length(palette) < 2) {
      stop("The palette must be a character vector with at least two colors.")
    }

    # Compute SSAO values
    message("Computing SSAO (Screen Space Ambient Occlusion)...")
    ssao_values <- compute_ssao(las, kernel_size = kernel_size, pixel_size = pixel_size,
                                num_samples = num_samples, ncpu = ncpu)

    # Normalize SSAO values
    normalized_values <- (ssao_values - min(ssao_values, na.rm = TRUE)) /
      (max(ssao_values, na.rm = TRUE) - min(ssao_values, na.rm = TRUE))

    # Map to colors
    colors <- grDevices::colorRampPalette(palette)(256)
    rgb_colors <- colors[as.numeric(cut(normalized_values, breaks = 256))]
    rgb_matrix <- grDevices::col2rgb(rgb_colors)

    las@data$R <- rgb_matrix[1, ]
    las@data$G <- rgb_matrix[2, ]
    las@data$B <- rgb_matrix[3, ]
  }

  return(las)
}

#' Compute PCV (Portion de Ciel Visible) for point cloud
#' @description Calculates ambient occlusion using sky visibility algorithm.
#' Based on Duguet & Girardeau-Montaut (2004).
#'
#' @param las LAS object
#' @param radius Numeric radius for neighborhood search
#' @param num_directions Integer number of directional rays to cast
#' @param ncpu Integer number of CPUs to use for parallel processing
#'
#' @return Numeric vector of PCV values (sky visibility) for each point
#'
#' @details
#' The PCV algorithm computes the visible portion of the sky from each point by:
#' \itemize{
#'   \item Casting rays in multiple directions around each point
#'   \item Computing the maximum elevation angle (horizon) in each direction
#'   \item Calculating the average sky visibility across all directions
#' }
#'
#' This function uses an optimized C++ implementation with OpenMP parallelization
#' for improved performance on large point clouds.
#'
#' @keywords internal
compute_pcv <- function(las, radius = 1.0, num_directions = 60, ncpu = 4) {
  # Call C++ implementation with S4 LAS object
  # Uses lidR::GridPartition spatial index for fast neighbor lookup
  pcv_values <- cppComputePCV(las, radius = radius, num_directions = num_directions, ncpu = ncpu)

  return(pcv_values)
}

#' Compute SSAO (Screen Space Ambient Occlusion) for point cloud
#' @description Fast ambient occlusion using 2D depth buffer technique.
#' Much faster than PCV as it works on projected depth maps.
#'
#' @param las LAS object
#' @param kernel_size Integer kernel size in pixels for sampling
#' @param pixel_size Numeric resolution of the depth map in spatial units
#' @param num_samples Integer number of samples per point
#' @param ncpu Integer number of CPUs to use for parallel processing
#'
#' @return Numeric vector of SSAO values (ambient occlusion) for each point
#'
#' @details
#' The SSAO algorithm computes ambient occlusion by:
#' \itemize{
#'   \item Projecting the point cloud to a 2D depth map (grid)
#'   \item For each point, sampling the depth buffer around it
#'   \item Calculating occlusion based on depth differences
#'   \item Applying distance and angle-based falloff
#' }
#'
#' This is significantly faster than full 3D PCV because it only requires
#' 2D image processing operations rather than 3D neighbor searches and ray tracing.
#'
#' @keywords internal
compute_ssao <- function(las, kernel_size = 5, pixel_size = 0.1, num_samples = 16, ncpu = 4) {
  # Call C++ implementation with S4 LAS object
  # Uses optimized depth map projection with pre-computed bounds
  ssao_values <- cppComputeSSAO(las, kernel_size = kernel_size, pixel_size = pixel_size,
                                num_samples = num_samples, ncpu = ncpu)

  return(ssao_values)
}

#' Download NAIP Imagery for LiDAR Extent
#' @description Downloads NAIP (National Agriculture Imagery Program) imagery from
#' Microsoft Planetary Computer STAC API for the extent of a LAS/LAZ point cloud.
#'
#' @param las A LAS object or path to a LAS/LAZ file
#' @param output_path Character string specifying output file path for the downloaded imagery.
#'   If NULL, creates a file named based on the input LAS file.
#' @param year_range Character vector of length 2 specifying date range for NAIP imagery
#'   in format c("YYYY-MM-DD", "YYYY-MM-DD"). Default is c("2018-01-01", "2023-12-31").
#' @param buffer Numeric value to buffer the extent in meters. Default is 0.
#' @param overwrite Logical, whether to overwrite existing output file. Default is FALSE.
#'
#' @return Character string of the output file path, or NULL if download failed.
#'
#' @details
#' This function queries the Microsoft Planetary Computer STAC API to find and download
#' NAIP imagery that overlaps with the extent of the input LAS file. The imagery is
#' automatically cropped to match the LiDAR extent and saved as a GeoTIFF.
#'
#' NAIP imagery is typically 4-band (RGB + NIR) with 0.6m or 1m resolution, collected
#' annually or biannually across the continental United States.
#'
#' Requires the \code{rstac} package for STAC API access.
#'
#' @examples
#' \dontrun{
#' # Load example LAS file
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package="spanner")
#' las <- readLAS(LASfile)
#' 
#' # Download NAIP for a LAS file
#' naip_path <- download_naip_for_las(las, output_path = "naip_imagery.tif")
#'
#' # Download with buffer and specific year range
#' naip_path <- download_naip_for_las(las, buffer = 10,
#'                                    year_range = c("2020-01-01", "2023-12-31"))
#'
#' # Then use with colorize_las
#' las_colored <- colorize_las(las, method = "rgb", raster_path = naip_path)
#' }
#'
#' @export
download_naip_for_las <- function(las, output_path = NULL,
                                   year_range = c("2018-01-01", "2023-12-31"),
                                   buffer = 0, overwrite = FALSE) {

  # Check for required packages
  if (!requireNamespace("rstac", quietly = TRUE)) {
    stop("Package 'rstac' is required for NAIP download. Install with: install.packages('rstac')")
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required. Install with: install.packages('terra')")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Load LAS if it's a file path
  if (is.character(las)) {
    if (!requireNamespace("lidR", quietly = TRUE)) {
      stop("Package 'lidR' is required to read LAS files. Install with: install.packages('lidR')")
    }
    message("Reading LAS file: ", las)
    las_obj <- lidR::readLAS(las)
  } else if (inherits(las, "LAS")) {
    las_obj <- las
  } else {
    stop("las must be a LAS object or path to a LAS/LAZ file")
  }

  # Get extent and convert to lat/lon
  las_extent <- sf::st_bbox(las_obj)
  las_crs <- sf::st_crs(las_obj)

  # Check if CRS is defined
  if (is.na(las_crs)) {
    message("LAS object has no CRS defined. Cannot download NAIP imagery.")
    message("Use lidR::st_crs(las) <- <EPSG_CODE> to assign a CRS first.")
    return(NULL)
  }

  # Apply buffer if specified
  if (buffer > 0) {
    las_extent[c("xmin", "ymin")] <- las_extent[c("xmin", "ymin")] - buffer
    las_extent[c("xmax", "ymax")] <- las_extent[c("xmax", "ymax")] + buffer
  }

  # Convert to WGS84 (EPSG:4326) for STAC query
  bbox_sf <- sf::st_as_sfc(las_extent, crs = las_crs)
  bbox_wgs84 <- sf::st_transform(bbox_sf, 4326)
  bbox <- sf::st_bbox(bbox_wgs84)

  message("Querying Microsoft Planetary Computer for NAIP imagery...")
  message("Extent (WGS84): ", paste(round(bbox, 5), collapse = ", "))

  # Connect to Microsoft Planetary Computer STAC
  stac_obj <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")

  # Search for NAIP imagery
  items <- tryCatch({
    search_obj <- rstac::stac_search(
      stac_obj,
      collections = "naip",
      bbox = as.numeric(bbox),
      datetime = paste(year_range, collapse = "/")
    )
    post_obj <- rstac::post_request(search_obj)
    rstac::items_fetch(post_obj)
  }, error = function(e) {
    stop("STAC query failed: ", e$message)
  })

  if (length(items$features) == 0) {
    message("No NAIP imagery found for this extent and date range.")
    return(NULL)
  }

  message("Found ", length(items$features), " NAIP scene(s)")

  # Get the most recent scene
  item <- items$features[[1]]
  message("Using scene from: ", item$properties$datetime)
  message("Scene ID: ", item$id)

  # Sign the URL for access
  signed_items <- rstac::items_sign(items, rstac::sign_planetary_computer())
  rgb_url <- paste0("/vsicurl/", signed_items$features[[1]]$assets$image$href)

  message("Downloading and cropping imagery...")

  # Read raster
  r <- terra::rast(rgb_url)

  # Convert LAS extent to raster CRS
  bbox_raster_crs <- sf::st_transform(bbox_sf, sf::st_crs(r))

  # Crop to extent
  r_crop <- terra::crop(r, bbox_raster_crs)

  # Generate output path if not provided
  if (is.null(output_path)) {
    if (is.character(las)) {
      base_name <- tools::file_path_sans_ext(basename(las))
      output_path <- file.path(dirname(las), paste0(base_name, "_NAIP.tif"))
    } else {
      output_path <- file.path(getwd(), "naip_imagery.tif")
    }
  }

  # Check if file exists
  if (file.exists(output_path) && !overwrite) {
    stop("Output file already exists. Set overwrite=TRUE to replace it.")
  }

  # Save
  terra::writeRaster(r_crop, output_path, overwrite = overwrite)

  message("NAIP imagery saved to: ", output_path)
  message("Bands: ", terra::nlyr(r_crop))
  message("Resolution: ", paste(terra::res(r_crop), collapse = " x "), " m")
  message("Dimensions: ", paste(dim(r_crop)[2:1], collapse = " x "), " pixels")

  return(invisible(output_path))
}

#' Merge RGB colors from two colorized LAS objects
#' @description Blends the RGB values from two LAS objects to create a new composite coloring.
#' Useful for combining different coloring methods (e.g., ambient occlusion with raster RGB).
#'
#' @param las1 First LAS object with R, G, B fields
#' @param las2 Second LAS object with R, G, B fields (must have same number of points as las1)
#' @param alpha Numeric value between 0 and 1 controlling the blend ratio.
#'   0 = all las1 colors, 1 = all las2 colors, 0.5 = equal blend. Default is 0.5.
#' @param method Character string specifying blend method: "alpha" for alpha blending,
#'   "multiply" for multiplicative blending, "screen" for screen blending,
#'   "overlay" for overlay blending. Default is "alpha".
#'
#' @return A LAS object (copy of las1) with merged R, G, and B fields
#'
#' @details
#' Blending methods:
#' \describe{
#'   \item{alpha}{Simple linear interpolation: (1-alpha)*las1 + alpha*las2}
#'   \item{multiply}{Multiplicative blend (darkens): (las1 * las2) / 255}
#'   \item{screen}{Screen blend (lightens): 255 - ((255-las1) * (255-las2)) / 255}
#'   \item{overlay}{Overlay blend: combines multiply and screen based on base color}
#' }
#'
#' Common use cases:
#' \itemize{
#'   \item Combine ambient occlusion (PCV/SSAO) with aerial RGB for realistic shading
#'   \item Blend attribute coloring with terrain colors
#'   \item Overlay multiple visualization layers
#' }
#'
#' @examples
#' \dontrun{
#' # Load example LAS file
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package="spanner")
#' las <- readLAS(LASfile)
#' 
#' # Combine SSAO ambient occlusion with aerial RGB
#' las_ao <- colorize_las(las, method="ssao", palette=c("black", "white"))
#' rgb_file <- system.file("extdata", "UAS_Clip_RGB.tif", package="spanner")
#' las_rgb <- colorize_las(las, method="rgb", raster_path=rgb_file)
#' las_merged <- merge_las_colors(las_ao, las_rgb, alpha=0.3, method="multiply")
#'
#' # Blend attribute coloring with RGB at 50/50
#' las_height <- colorize_las(las, method="attr", attribute_name="Z",
#'                            palette=c("blue", "red"))
#' las_merged <- merge_las_colors(las_height, las_rgb, alpha=0.5)
#' }
#'
#' @export
merge_las_colors <- function(las1, las2, alpha = 0.5, method = "alpha") {

  # Validate inputs
  if (!inherits(las1, "LAS") || !inherits(las2, "LAS")) {
    stop("Both las1 and las2 must be LAS objects")
  }

  if (nrow(las1@data) != nrow(las2@data)) {
    stop("las1 and las2 must have the same number of points")
  }

  if (!all(c("R", "G", "B") %in% names(las1@data))) {
    stop("las1 must have R, G, B fields. Use colorize_las() first.")
  }

  if (!all(c("R", "G", "B") %in% names(las2@data))) {
    stop("las2 must have R, G, B fields. Use colorize_las() first.")
  }

  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be a numeric value between 0 and 1")
  }

  method <- tolower(method)
  if (!method %in% c("alpha", "multiply", "screen", "overlay")) {
    stop("method must be one of 'alpha', 'multiply', 'screen', or 'overlay'")
  }

  # Create output LAS (copy of las1)
  las_out <- las1

  # Extract RGB values
  r1 <- as.numeric(las1@data$R)
  g1 <- as.numeric(las1@data$G)
  b1 <- as.numeric(las1@data$B)

  r2 <- as.numeric(las2@data$R)
  g2 <- as.numeric(las2@data$G)
  b2 <- as.numeric(las2@data$B)

  # Apply blending method
  if (method == "alpha") {
    # Simple alpha blending
    r_out <- (1 - alpha) * r1 + alpha * r2
    g_out <- (1 - alpha) * g1 + alpha * g2
    b_out <- (1 - alpha) * b1 + alpha * b2

  } else if (method == "multiply") {
    # Multiplicative blending (darkens)
    r_out <- (r1 * r2) / 255
    g_out <- (g1 * g2) / 255
    b_out <- (b1 * b2) / 255

  } else if (method == "screen") {
    # Screen blending (lightens)
    r_out <- 255 - ((255 - r1) * (255 - r2)) / 255
    g_out <- 255 - ((255 - g1) * (255 - g2)) / 255
    b_out <- 255 - ((255 - b1) * (255 - b2)) / 255

  } else if (method == "overlay") {
    # Overlay blending (combines multiply and screen)
    overlay_channel <- function(base, blend) {
      ifelse(base < 128,
             2 * base * blend / 255,
             255 - 2 * (255 - base) * (255 - blend) / 255)
    }
    r_out <- overlay_channel(r1, r2)
    g_out <- overlay_channel(g1, g2)
    b_out <- overlay_channel(b1, b2)
  }

  # Clamp values to 0-255 range
  las_out@data$R <- as.integer(pmin(pmax(r_out, 0), 255))
  las_out@data$G <- as.integer(pmin(pmax(g_out, 0), 255))
  las_out@data$B <- as.integer(pmin(pmax(b_out, 0), 255))

  return(las_out)
}

#' Create animated GIF of rotating 3D point cloud
#' @description Generates a 360-degree rotating animation of a LAS point cloud using rgl
#' and saves it as an animated GIF.
#'
#' @param las LAS object to visualize. Should have R, G, B fields for color.
#' @param output_path Character string specifying output GIF file path.
#'   Default is "pointcloud_rotation.gif".
#' @param duration Numeric duration of the animation in seconds. Default is 12.
#' @param rpm Numeric rotations per minute for the spin. Default is 5.
#' @param background Character string specifying background color. Default is "white".
#' @param axis Character specifying rotation axis: "z" for vertical rotation (default),
#'   "x" for horizontal rotation, "y" for front-to-back rotation.
#' @param show_axis Logical whether to show axes. Default is TRUE.
#' @param show_legend Logical whether to show legend. Default is TRUE.
#' @param screen_size Numeric vector of length 2 specifying window dimensions as c(width, height).
#'   Default is c(800, 600).
#' @param overwrite Logical whether to overwrite existing output file. Default is FALSE.
#'
#' @return Character string of the output file path (invisible)
#'
#' @details
#' This function creates a smooth 360-degree rotation animation by:
#' \itemize{
#'   \item Plotting the point cloud using lidR's plot function with RGB colors
#'   \item Using rgl's movie3d and spin3d to create smooth rotation
#'   \item Saving the result as an animated GIF
#' }
#'
#' The rotation speed is controlled by the rpm (rotations per minute) parameter.
#' The total duration determines how long the animation will be.
#'
#' Requires the \code{rgl} package and lidR for plotting.
#'
#' @examples
#' \dontrun{
#' # Load example LAS file
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package="spanner")
#' las <- readLAS(LASfile)
#' 
#' # Create basic rotation GIF with attribute coloring
#' las_colored <- colorize_las(las, method="attr", attribute_name="Z")
#' create_rotation_gif(las_colored, output_path="rotation.gif")
#'
#' # High quality with specific settings
#' create_rotation_gif(las_colored,
#'                     output_path="highres_rotation.gif",
#'                     duration=15,
#'                     rpm=10,
#'                     background="black",
#'                     show_axis=FALSE,
#'                     show_legend=FALSE)
#'
#' # Rotate around X axis for side-to-side view
#' create_rotation_gif(las_colored, axis="x")
#' }
#'
#' @export
create_rotation_gif <- function(las,
                                output_path = "pointcloud_rotation.gif",
                                duration = 12,
                                rpm = 5,
                                background = "white",
                                axis = "z",
                                show_axis = TRUE,
                                show_legend = TRUE,
                                screen_size = c(800, 600),
                                overwrite = FALSE) {

  # Check for required packages
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required. Install with: install.packages('rgl')")
  }
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop("Package 'lidR' is required. Install with: install.packages('lidR')")
  }
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("Package 'magick' is required. Install with: install.packages('magick')")
  }

  # Validate inputs
  if (!inherits(las, "LAS")) {
    stop("las must be a LAS object")
  }

  # Validate screen_size
  if (!is.numeric(screen_size) || length(screen_size) != 2) {
    stop("screen_size must be a numeric vector of length 2 (e.g., c(800, 600))")
  }

  if (file.exists(output_path) && !overwrite) {
    stop("Output file already exists. Set overwrite=TRUE to replace it.")
  }

  axis <- tolower(axis)
  if (!axis %in% c("x", "y", "z")) {
    stop("axis must be one of 'x', 'y', or 'z'")
  }

  # Check for RGB colors
  has_rgb <- all(c("R", "G", "B") %in% names(las@data))

  if (!has_rgb) {
    warning("No RGB fields found in LAS object. Using default coloring. Use colorize_las() first for custom colors.")
  }

  message("Creating rotation animation...")
  message(sprintf("  Duration: %d seconds", duration))
  message(sprintf("  Rotation speed: %d RPM", rpm))
  message(sprintf("  Rotation axis: %s", toupper(axis)))

  # Plot the LAS object using lidR's default plot settings
  # Only specify color="RGB" if RGB fields exist
  if (has_rgb) {
    lidR::plot(las, color = "RGB")
  } else {
    lidR::plot(las)
  }

  # Set window size after plot is created
  rgl::par3d(windowRect = c(0, 0, screen_size[1], screen_size[2]))

  # Set up rotation axis vector
  axis_vector <- switch(axis,
                       "x" = c(1, 0, 0),
                       "y" = c(0, 1, 0),
                       "z" = c(0, 0, 1))

  # Extract just the filename without extension for PNG files
  output_name <- tools::file_path_sans_ext(basename(output_path))
  output_dir_path <- dirname(output_path)

  # If output_path is just a filename, use current directory
  if (output_dir_path == ".") {
    output_dir_path <- getwd()
  }

  # Create temporary directory for frames
  temp_dir <- tempfile()
  dir.create(temp_dir)

  message("Rendering animation...")

  # Calculate number of frames and angles
  n_frames <- ceiling(duration * 10)  # 10 fps for reasonable file size

  message(sprintf("Capturing %d frames...", n_frames))

  # Use play3d to capture frames automatically
  # This is simpler and more reliable than manual capture
  spin_function <- rgl::spin3d(axis = axis_vector, rpm = rpm)

  # Capture frames using play3d
  # Use webshot=FALSE to avoid Chrome/webshot2 errors and use direct rgl snapshot
  rgl::movie3d(
    spin_function,
    duration = duration,
    fps = 10,
    dir = temp_dir,
    movie = "temp",
    clean = FALSE,  # Don't clean up - we need the PNGs
    verbose = FALSE,
    webshot = FALSE  # Use rgl.snapshot() directly, avoid webshot2/Chrome
  )

  # Close rgl window
  rgl::close3d()

  message("Combining frames into GIF...")

  # Read frames with magick
  frame_files <- list.files(temp_dir, pattern = "*.png", full.names = TRUE)
  frame_files <- sort(frame_files)

  frames <- magick::image_read(frame_files)

  # Create animation
  fps <- n_frames / duration
  frames_animated <- magick::image_animate(frames, fps = fps, dispose = "previous")

  # Construct output path
  actual_output <- file.path(output_dir_path, paste0(output_name, ".gif"))

  # Write GIF
  magick::image_write(frames_animated, actual_output, format = "gif")

  # Cleanup temporary directory
  unlink(temp_dir, recursive = TRUE)

  message("Animation saved to: ", actual_output)
  message("Duration: ", duration, " seconds")
  if (file.exists(actual_output)) {
    message("File size: ", round(file.size(actual_output) / 1024 / 1024, 2), " MB")
  }

  return(invisible(actual_output))
}

