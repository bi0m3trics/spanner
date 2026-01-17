## usethis namespace: start
#' @useDynLib spanner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats complete.cases density dist na.omit princomp quantile
#' @importFrom sf st_as_sf st_coordinates
#' @importFrom terra classify focal ifel rast res values<-
#' @importFrom lidR filter_poi grid_metrics random_per_voxel LAS
#' @importFrom mathjaxr preview_rd
## usethis namespace: end

# Declare global variables to avoid R CMD check NOTEs for non-standard evaluation
utils::globalVariables(c("X", "Y", "Z", "TreeID", "treeID", "Classification", 
                         "verticality", "f", "p2r"))

# Note: The src/ directory contains C++ header files (.hpp) alongside source files (.cpp).
# These header files (algorithms.hpp, classes.hpp, methods.hpp, utils.hpp, and optim.hpp)
# are required for compilation and contain template definitions, inline functions, and
# class/function declarations used by the corresponding .cpp implementation files.
# The optim.hpp file is part of the OptimLib header-only optimization library.
