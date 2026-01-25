## usethis namespace: start
#' @useDynLib spanner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases density dist na.omit princomp quantile
#' @importFrom sf st_as_sf st_coordinates
#' @importFrom terra classify focal ifel rast res values<-
#' @importFrom lidR filter_poi pixel_metrics random_per_voxel LAS
#' @importFrom mathjaxr preview_rd
## usethis namespace: end

# Declare global variables to avoid R CMD check NOTEs for non-standard evaluation
utils::globalVariables(c("X", "Y", "Z", "TreeID", "treeID", "Classification", 
                         "verticality", "null", "f", "p2r"))

# Configure lidR to use terra package for raster operations
.onLoad <- function(libname, pkgname) {
  options(lidR.raster.default = "terra")
}

