## usethis namespace: start
#' @useDynLib spanner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases density dist na.omit princomp quantile
#' @importFrom sf st_as_sf st_coordinates
#' @importFrom terra classify focal ifel rast res values<-
#' @importFrom lidR filter_poi grid_metrics random_per_voxel LAS
## usethis namespace: end

# Declare global variables to avoid R CMD check NOTEs for non-standard evaluation
utils::globalVariables(c("X", "Y", "Z", "TreeID", "treeID", "Classification", 
                         "verticality", "null", "f", "p2r"))
