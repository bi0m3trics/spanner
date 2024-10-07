# Define a function to create a circular kernel for focal operations
getCircularKernel <- function(radius)
{
  # Calculate the side length of the kernel
  kernel_side <- 2 * as.integer(radius) + 1
  # Create a matrix representing the x- and y-coordinates of the kernel
  kernel_y <- matrix(rep(radius:-radius, kernel_side), ncol=kernel_side)
  kernel_x <- -t(kernel_y)
  # Calculate the Euclidean distance from the center for each cell in the kernel
  kernel   <- matrix(as.matrix(dist(cbind(as.vector(kernel_x), as.vector(kernel_y))))[as.integer((kernel_side^2) / 2) + 1,], ncol=kernel_side)
  # Set cells within the radius to 0 and cells outside the radius to 1
  kernel[kernel <= radius] <- 0
  kernel[kernel > 0]  <- 1
  # Invert the kernel so that cells within the radius are 1 and cells outside are 0
  kernel <- 1 - kernel
  # Return the kernel
  return(kernel)
}

#' Process rasters based on suitability, gap, and spur parameters
#'
#' `process_rasters_patchmorph` processes an input raster by reclassifying it based on suitability levels,
#' and then applying gap and spur distance transformations to generate a list of processed rasters.
#'
#' @param input_raster RasterLayer The input raster to be processed.
#' @param suitList numeric A vector of suitability levels for reclassification.
#' @param gapList numeric A vector of gap distances for processing.
#' @param spurList numeric A vector of spur distances for processing.
#' @return list A list of processed rasters with names indicating the suitability, gap, and spur parameters used.
#'
#' @examples
#' \dontrun{
#' # Define input parameters
#' las <- lidR::readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))
#' input_raster <- lidR::rasterize_canopy(las, res = 1, lidR::pitfree(c(0,2,5,10,15), c(0, 2)))
#' suitList <- c(0, 2, 32)
#' gapList <- seq(1, 8, by = 1)
#' spurList <- seq(1, 8, by = 1)
#'
#' # Process the rasters
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList)
#'
#' # Plot the first processed raster
#' plot(processed_rasters[[1]])
#' }
#'
#' @export
process_rasters_patchmorph <- function(input_raster, suitList, gapList, spurList) {

  # Initialize a list to store the processed rasters
  rasters <- list()

  # Loop over the suitability levels
  for (suit in suitList) {
    # Reclassify the raster based on the current suitability level
    rSuit <- classify(input_raster, matrix(c(-Inf, suit, NA, suit, Inf, 1), ncol=3, byrow=TRUE), right = FALSE)

    # Loop over the gap distances
    for (cntGap in gapList) {
      gapDist <- cntGap

      # Loop over the spur distances
      for (cntSpur in spurList) {
        spurDist <- cntSpur

        # Calculate the Euclidean distance for each cell in the suitability raster
        rDist <- terra::distance(rSuit)
        # plot(rDist)
        # Create a circular kernel based on the spur distance
        spurKernal <- getCircularKernel(spurDist / res(input_raster)[1])

        # Apply the maximum focal operation to the distance raster using the spur kernel
        rFocal <- focal(rDist, w = spurKernal, fun = max, na.rm = TRUE)
        # plot(rFocal)
        # Reclassify the raster based on the spur distance
        rSpur <- classify(rFocal, matrix(c(-Inf, spurDist/2, NA, spurDist/2, Inf, 1), ncol=3, byrow=TRUE), right = FALSE)
        # plot(rSpur)
        # Calculate the Euclidean distance for each cell in the spur raster
        rDist2 <- terra::distance(rSpur)
        # plot(rDist2)
        # Create a circular kernel based on the gap distance
        gapKernal <- getCircularKernel(gapDist / res(input_raster)[1])

        # Apply the maximum focal operation to the distance raster using the gap kernel
        rFocal2 <- focal(rDist2, w = gapKernal, fun = max, na.rm = TRUE)
        # plot(rFocal2)
        # Reclassify the raster based on the gap distance
        rGap <- classify(rFocal2, matrix(c(-Inf, gapDist/2, NA, gapDist/2, Inf, 1), ncol=3, byrow=TRUE), right = FALSE)
        # plot(rGap)
        # Print a message indicating the current parameters
        message(paste("Processing Suitability:", suit, "Gap:", cntGap, "Spur:", cntSpur))

        # Assign a name to the raster based on the current parameters
        names(rGap) <- paste("suit", suit, "gap", cntGap, "spur", cntSpur, sep = "_")
        # plot(rGap)
        # Add the raster to the list
        rasters[[length(rasters) + 1]] <- rGap
      }
    }
  }

  # Return the list of processed rasters
  return(rasters)
}

#' Plot a raster by its name
#'
#' `plot_raster_by_name` plots a raster from a list of rasters based on the provided raster name.
#'
#' @param rasters list A list of rasters.
#' @param raster_name character The name of the raster to be plotted.
#' @return NULL This function does not return a value. It plots the raster if found.
#'
#' @examples
#' \dontrun{
#' # Define input parameters
#' las <- lidR::readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))
#' input_raster <- lidR::rasterize_canopy(las, res = 1, lidR::pitfree(c(0,2,5,10,15), c(0, 2)))
#' suitList <- c(0, 2, 32)
#' gapList <- seq(1, 8, by = 1)
#' spurList <- seq(1, 8, by = 1)
#'
#' # Process the rasters
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList)
#'
#' # Plot a raster by its name
#' plot_raster_by_name(processed_rasters, "suit_2_gap_2_spur_6")
#' }
#'
#' @export
plot_raster_by_name <- function(rasters, raster_name) {
  # Find the index of the raster in the list by its name
  index <- which(sapply(rasters, function(r) raster_name %in% names(r)))

  # Check if the raster with the given name exists in the list
  if (length(index) > 0) {
    # If it exists, select the raster from the list using its index
    r <- rasters[[index]]

    # Plot the selected raster
    plot(r)
  } else {
    # If it doesn't exist, print a message indicating that no raster was found with the given name
    print(paste("No raster found with name", raster_name))
  }
}

#' Sum rasters by suitability level
#'
#' `sum_rasters_by_suitability` sums rasters from a list based on their suitability levels.
#'
#' @param rasters list A list of rasters.
#' @param suitList numeric A vector of suitability levels.
#' @return list A list of summed rasters for each suitability level.
#'
#' @examples
#' \dontrun{
#' # Define input parameters
#' las <- lidR::readLAS(system.file("extdata", "MixedConifer.laz", package="lidR"))
#' input_raster <- lidR::rasterize_canopy(las, res = 1, lidR::pitfree(c(0,2,5,10,15), c(0, 2)))
#' suitList <- c(0, 2, 32)
#' gapList <- seq(1, 8, by = 1)
#' spurList <- seq(1, 8, by = 1)
#'
#' # Process the rasters
#' processed_rasters <- process_rasters_patchmorph(input_raster, suitList, gapList, spurList)
#'
#'
#' Sum rasters by suitability level
#' summed_rasters <- sum_rasters_by_suitability(processed_rasters, suitList)
# '
#' # Call the plot_raster_by_name function to plot the raster named "suit_2_sum"
#' plot_raster_by_name(rSumList, "suit_2_sum")
#' }
#'
#' @export
sum_rasters_by_suitability <- function(rasters, suitList) {
  # Initialize an empty list to store the summed rasters
  rSumList <- list()

  # Loop over the unique suitability levels
  for (suit in unique(suitList)) {
    # Filter the rasters in the list by the current suitability level
    suitRasters <- rasters[sapply(rasters, function(x) any(grepl(paste0("suit_", suit), names(x))))]

    # Check if there are any rasters for the current suitability level
    if (length(suitRasters) > 0) {
      # If there are, initialize a SpatRaster with the same dimensions as the input rasters
      rSum <- rast(suitRasters[[1]])

      # Set all values in the SpatRaster to 0
      values(rSum) <- 0

      # Loop over the rasters for the current suitability level
      for (i in seq_along(suitRasters)) {
        # Add the values of the current raster to the sum, replacing NaN values with 0
        rSum <- rSum + ifel(is.na(suitRasters[[i]]), 0, suitRasters[[i]])
      }

      # Assign a name to the summed raster based on the current suitability level
      names(rSum) <- paste("suit", suit, "sum", sep = "_")

      # Add the summed raster to the list
      rSumList[[length(rSumList) + 1]] <- rSum
    }
  }

  # Return the list of summed rasters
  return(rSumList)
}
