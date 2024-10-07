#' A simple function to download some example MLS/TLS data in a ponderosa pine forest north of Flagstaff AZ
#'
#' Note: none of these datasets have been downsampled, classified for ground points, or normalized, all of which
#' must be done before using them with 'Raster_Eigen_TreeLocations' or 'Segment_Graph' functions in the 'spanner' package.
#'
#' @param character The basename of the example dataset to be downloaded. Possibilities include:
#' * "PineExampleA" - A 1ha exmaple dataset pre-downsampled using a 0.01 cm voxel grid
#' * "DensePatchA" - A small, dense ponderosa pine patch with many, small trees
#' * "DensePatchB"  - A small, dense ponderosa pine patch with few, medium-sized trees
#' * "SparsePatchA"  - A small, open ponderosa pine patch with few, medium-sized trees
#' * "SparsePatchB" - A small, open ponderosa pine patch with few, largely spaced trees
#' * "ZebcamExample" - A big, RGB colored dataset
#'
#' @return Nothing Downloads and saves a specified .laz file to the "extdata" folder of the 'spanner' package
#'
#' @examples
#'
#'
#' \dontrun{
#' # Specify one of six possible datasets to download
#' getExampleData("SparsePatchA")
#' }
#' @export
getExampleData <- function(lazName=c("PineExampleA", "DensePatchA","DensePatchB","SparsePatchA","SparsePatchB","ZebcamExample")) {
  # Define the spannerPath
  spannerPath <- find.package("spanner", lib.loc = NULL, quiet = FALSE, verbose = getOption("verbose"))

  # Define the extdata path
  extdataPath <- file.path(spannerPath, "extdata")

  # Check if extdata directory exists, if not create it
  if (!dir.exists(extdataPath)) {
    dir.create(extdataPath, recursive = TRUE)
  }

  # Construct the URL
  url <- paste0("http://quantitativeecology.org/files/spanner/", lazName, ".laz")

  # Define the destination file path
  destfile <- file.path(extdataPath, paste0(lazName, ".laz"))

  # Try to download the file using a valid method
  tryCatch({
    download.file(url, destfile, method = "libcurl", mode = "wb")
    # Return the path to the downloaded file
    return(destfile)
  }, error = function(e) {
    message("Failed to download the file: ", e$message)
    return(NULL)
  })
}
