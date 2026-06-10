# Utility helpers -----------------------------------------------------------

.assert_scalar_numeric <- function(x, name, lower = -Inf, lower_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be a numeric scalar.", name), call. = FALSE)
  }
  if (lower_open) {
    if (!(x > lower)) {
      stop(sprintf("'%s' must be > %s.", name, lower), call. = FALSE)
    }
  } else {
    if (!(x >= lower)) {
      stop(sprintf("'%s' must be >= %s.", name, lower), call. = FALSE)
    }
  }
}

.assert_scalar_integerish <- function(x, name, lower = 1L) {
  .assert_scalar_numeric(x, name, lower = lower)
  if (abs(x - round(x)) > 1e-8) {
    stop(sprintf("'%s' must be an integer-like scalar.", name), call. = FALSE)
  }
}

.extract_xyz <- function(las) {
  data <- las@data
  if (!all(c("X", "Y", "Z") %in% names(data))) {
    stop("LAS must contain X, Y, and Z coordinates.", call. = FALSE)
  }
  list(
    X = as.numeric(data$X),
    Y = as.numeric(data$Y),
    Z = as.numeric(data$Z)
  )
}

.match_seed_z <- function(seed_x, seed_y, xyz) {
  if (length(seed_x) == 0L) {
    return(numeric())
  }

  point_xy <- cbind(xyz$X, xyz$Y)
  out <- numeric(length(seed_x))
  for (i in seq_along(seed_x)) {
    dx <- point_xy[, 1] - seed_x[i]
    dy <- point_xy[, 2] - seed_y[i]
    j <- which.min(dx * dx + dy * dy)
    out[i] <- xyz$Z[j]
  }
  out
}

.derive_seeds_lmf <- function(las, hmin) {
  ws_fun <- function(z) {
    pmin(8, pmax(2, 0.08 * z + 2))
  }

  ttops <- tryCatch(
    lidR::locate_trees(las, lidR::lmf(ws = ws_fun, hmin = hmin)),
    error = function(e) NULL
  )

  if (is.null(ttops)) {
    return(data.frame(X = numeric(), Y = numeric(), Z = numeric()))
  }
  if (nrow(ttops) == 0L) {
    return(data.frame(X = numeric(), Y = numeric(), Z = numeric()))
  }

  coords <- sf::st_coordinates(ttops)
  out <- data.frame(
    X = as.numeric(coords[, 1]),
    Y = as.numeric(coords[, 2]),
    Z = rep(NA_real_, nrow(coords))
  )

  z_name <- intersect(c("Z", "z", "Height", "height", "tree_height"), names(ttops))
  if (length(z_name) > 0L) {
    out$Z <- as.numeric(ttops[[z_name[1L]]])
  }

  out
}

.coerce_seed_xyz <- function(seeds, las, hmin) {
  xyz <- .extract_xyz(las)

  if (is.null(seeds)) {
    seeds_df <- .derive_seeds_lmf(las, hmin = hmin)
  } else if (inherits(seeds, "sf")) {
    coords <- sf::st_coordinates(seeds)
    seeds_df <- data.frame(
      X = as.numeric(coords[, 1]),
      Y = as.numeric(coords[, 2]),
      Z = NA_real_
    )
    z_name <- intersect(c("Z", "z", "Height", "height", "tree_height"), names(seeds))
    if (length(z_name) > 0L) {
      seeds_df$Z <- as.numeric(seeds[[z_name[1L]]])
    }
  } else if (is.data.frame(seeds)) {
    if (!all(c("X", "Y") %in% names(seeds))) {
      stop("When 'seeds' is a data.frame, columns X and Y are required.", call. = FALSE)
    }
    seeds_df <- data.frame(
      X = as.numeric(seeds$X),
      Y = as.numeric(seeds$Y),
      Z = if ("Z" %in% names(seeds)) as.numeric(seeds$Z) else NA_real_
    )
  } else {
    stop("'seeds' must be NULL, an sf object, or a data.frame with X/Y columns.", call. = FALSE)
  }

  if (nrow(seeds_df) == 0L) {
    return(seeds_df)
  }

  ok <- is.finite(seeds_df$X) & is.finite(seeds_df$Y)
  seeds_df <- seeds_df[ok, , drop = FALSE]
  if (nrow(seeds_df) == 0L) {
    return(seeds_df)
  }

  miss_z <- !is.finite(seeds_df$Z)
  if (any(miss_z)) {
    seeds_df$Z[miss_z] <- .match_seed_z(seeds_df$X[miss_z], seeds_df$Y[miss_z], xyz)
  }

  seeds_df <- seeds_df[is.finite(seeds_df$Z) & seeds_df$Z >= hmin, , drop = FALSE]
  if (nrow(seeds_df) == 0L) {
    return(seeds_df)
  }

  ord <- order(-seeds_df$Z, seeds_df$X, seeds_df$Y)
  seeds_df <- seeds_df[ord, , drop = FALSE]
  rownames(seeds_df) <- NULL
  seeds_df
}

.check_algorithm_output <- function(ids, n) {
  if (!is.integer(ids)) {
    storage.mode(ids) <- "integer"
  }
  if (length(ids) != n) {
    stop("Internal error: backend returned an ID vector with invalid length.", call. = FALSE)
  }
  bad <- !(is.na(ids) | ids > 0L)
  if (any(bad)) {
    stop("Internal error: backend returned non-positive IDs.", call. = FALSE)
  }
  ids
}

.resolve_profile_id <- function(profile) {
  if (is.function(profile)) {
    # Custom functions are currently accepted at the API level but
    # mapped to the default profile in C++ for determinism/performance.
    return(1L)
  }

  profile <- match.arg(profile, c("beta", "parabolic", "cone"))
  switch(
    profile,
    beta = 1L,
    parabolic = 2L,
    cone = 3L
  )
}

#' Allometric Li-Style Geodesic Segmentation
#'
#' Builds a point-cloud-based individual tree segmentation algorithm object
#' compatible with [lidR::segment_trees()]. The method is inspired by Li et al.
#' (2012), with allometric crown constraints, moving crown centers, and
#' geodesic growing over a local kNN graph.
#'
#' @param seeds Optional seed treetops as `sf` points or `data.frame` with `X`,
#'   `Y`, and optional `Z`.
#' @param hmin Minimum point height to consider for assignment.
#' @param k Number of graph neighbors per point.
#' @param max_jump Maximum local jump distance used in graph construction.
#' @param bin_height Vertical bin size guiding top-down growth.
#' @param crown_a Allometric coefficient in crown radius envelope.
#' @param crown_b Allometric exponent in crown radius envelope.
#' @param max_crown_radius Optional hard cap for crown radius (`NULL` for none).
#' @param center_drift Maximum per-step drift for moving crown center.
#' @param density_radius Radius used to estimate local density gaps.
#' @param gap_penalty Weight for low-density/gap transitions.
#' @param geodesic_threshold Maximum accumulated geodesic cost for assignment.
#' @param use_connected_components If `TRUE`, penalize transitions across sparse
#'   connected components.
#' @param min_component_size Minimum in-crown connected component size to keep.
#' @param omp Logical or `NULL`. If `NULL` (default), OpenMP is enabled on
#'   Windows/Linux and disabled on macOS. If `TRUE`/`FALSE`, this overrides the
#'   OS default.
#' @param crown_profile Named allometric profile (`"beta"`, `"parabolic"`,
#'   `"cone"`) or a custom function.
#' @param verbose If `TRUE`, emit high-level progress messages in R.
#'
#' @return A `PointCloudBased` `IndividualTreeSegmentation` algorithm object.
#' @export
#' @examples
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package = "spanner")
#' las <- lidR::readLAS(LASfile, select = "xyz")
#' las <- lidR::decimate_points(las, lidR::random_per_voxel(0.33, 1L))
#' las <- lidR::classify_ground(las, lidR::ptd())
#' las <- lidR::normalize_height(las, lidR::tin())
#' las_geodesic <- lidR::segment_trees(
#'   las,
#'   allometric_li_geodesic(
#'     seeds = NULL,
#'     hmin = 2,
#'     k = 12,
#'     max_jump = 1.5,
#'     bin_height = 1,
#'     crown_a = 0.35,
#'     crown_b = 0.65,
#'     max_crown_radius = NULL,
#'     center_drift = 1.5,
#'     density_radius = 1,
#'     gap_penalty = 3,
#'     geodesic_threshold = 5,
#'     use_connected_components = TRUE,
#'     min_component_size = 10,
#'     crown_profile = c("beta", "parabolic", "cone"),
#'     verbose = FALSE
#'   )
#' )
#' lidR::plot(las_geodesic, color = "treeID")
allometric_li_geodesic <- function(
  seeds = NULL,
  hmin = 2,
  k = 12,
  max_jump = 1.5,
  bin_height = 1.0,
  crown_a = 0.35,
  crown_b = 0.65,
  max_crown_radius = NULL,
  center_drift = 1.5,
  density_radius = 1.0,
  gap_penalty = 3,
  geodesic_threshold = 5,
  use_connected_components = TRUE,
  min_component_size = 3,
  omp = NULL,
  crown_profile = c("beta", "parabolic", "cone"),
  verbose = FALSE
) {
  .assert_scalar_numeric(hmin, "hmin", lower = 0)
  .assert_scalar_integerish(k, "k", lower = 2)
  .assert_scalar_numeric(max_jump, "max_jump", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(bin_height, "bin_height", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(crown_a, "crown_a", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(crown_b, "crown_b", lower = 0, lower_open = TRUE)
  if (!is.null(max_crown_radius)) {
    .assert_scalar_numeric(max_crown_radius, "max_crown_radius", lower = 0, lower_open = TRUE)
  }
  .assert_scalar_numeric(center_drift, "center_drift", lower = 0)
  .assert_scalar_numeric(density_radius, "density_radius", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(gap_penalty, "gap_penalty", lower = 0)
  .assert_scalar_numeric(geodesic_threshold, "geodesic_threshold", lower = 0, lower_open = TRUE)
  .assert_scalar_integerish(min_component_size, "min_component_size", lower = 1)
  if (!is.logical(use_connected_components) || length(use_connected_components) != 1L) {
    stop("'use_connected_components' must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.null(omp) && (!is.logical(omp) || length(omp) != 1L || is.na(omp))) {
    stop("'omp' must be NULL or TRUE/FALSE.", call. = FALSE)
  }

  sysname <- tolower(Sys.info()[["sysname"]])
  omp_default <- !identical(sysname, "darwin")
  use_omp <- if (is.null(omp)) omp_default else isTRUE(omp)

  profile_id <- .resolve_profile_id(crown_profile)
  max_cr <- if (is.null(max_crown_radius)) NA_real_ else as.numeric(max_crown_radius)

  f <- function(las) {
    if (!methods::is(las, "LAS")) {
      stop("This algorithm expects a lidR LAS object.", call. = FALSE)
    }

    n <- lidR::npoints(las)
    if (n == 0L) {
      return(integer())
    }

    xyz <- .extract_xyz(las)
    if (!all(is.finite(xyz$X)) || !all(is.finite(xyz$Y)) || !all(is.finite(xyz$Z))) {
      stop("Point cloud contains non-finite coordinates; clean input before allometric_li_geodesic().", call. = FALSE)
    }

    seeds_df <- .coerce_seed_xyz(seeds, las, hmin = hmin)

    if (nrow(seeds_df) == 0L) {
      return(rep(NA_integer_, n))
    }

    if (!all(is.finite(seeds_df$X)) || !all(is.finite(seeds_df$Y)) || !all(is.finite(seeds_df$Z))) {
      stop("Seed coordinates contain non-finite values.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message(sprintf("allometric_li_geodesic(): %d points, %d seeds", n, nrow(seeds_df)))
    }

    ids <- .Call(
      "cpp_allometric_li_geodesic",
      xyz$X,
      xyz$Y,
      xyz$Z,
      seeds_df$X,
      seeds_df$Y,
      seeds_df$Z,
      as.numeric(hmin),
      as.integer(k),
      as.numeric(max_jump),
      as.numeric(bin_height),
      as.numeric(crown_a),
      as.numeric(crown_b),
      as.numeric(max_cr),
      as.numeric(center_drift),
      as.numeric(density_radius),
      as.numeric(gap_penalty),
      as.numeric(geodesic_threshold),
      as.integer(isTRUE(use_connected_components)),
      as.integer(min_component_size),
      as.integer(profile_id),
      PACKAGE = "spanner"
    )

    .check_algorithm_output(ids, n)
  }

  lidR::plugin_its(f, omp = use_omp, raster_based = FALSE)
}

#' Allometric Random-Walker Segmentation
#'
#' Builds a seeded graph-diffusion style algorithm object for
#' [lidR::segment_trees()]. Treetops are fixed labels, while unlabeled points
#' iteratively propagate labels using weighted neighborhood support and
#' allometric penalties. The implementation uses a sparse local-candidate loop
#' to keep runtime practical on larger point clouds.
#'
#' @inheritParams allometric_li_geodesic
#' @param alpha Distance weight.
#' @param beta Vertical-separation weight.
#' @param gamma Eigenfeature-difference weight.
#' @param delta Allometry penalty weight.
#' @param eta Density-gap weight.
#' @param eigen_radius Radius for local eigenfeature estimation.
#' @param probability_threshold Minimum winning label confidence required to
#'   assign a point.
#' @param max_iterations Maximum number of label-propagation iterations.
#' @param tolerance Relative convergence tolerance.
#'
#' For faster runs, reduce `max_iterations` and/or increase `tolerance`.
#' For denser data, smaller `k` values often trade a little quality for a
#' meaningful speedup.
#'
#' @return A `PointCloudBased` `IndividualTreeSegmentation` algorithm object.
#' @export
#' @examples
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package = "spanner")
#' las <- lidR::readLAS(LASfile, select = "xyz")
#' las <- lidR::decimate_points(las, lidR::random_per_voxel(0.33, 1L))
#' las <- lidR::classify_ground(las, lidR::ptd())
#' las <- lidR::normalize_height(las, lidR::tin())
#' las_rw <- lidR::segment_trees(
#'   las,
#'   allometric_random_walker(
#'     seeds = NULL,
#'     hmin = 2,
#'     k = 12,
#'     alpha = 1,
#'     beta = 0.5,
#'     gamma = 0.5,
#'     delta = 1.5,
#'     eta = 1,
#'     crown_a = 0.35,
#'     crown_b = 0.65,
#'     eigen_radius = 1,
#'     density_radius = 1,
#'     probability_threshold = 0.35,
#'     max_iterations = 150,
#'     tolerance = 1e-05,
#'     crown_profile = c("beta", "parabolic", "cone"),
#'     verbose = FALSE
#'   )
#' )
#' lidR::plot(las_rw, color = "treeID")
allometric_random_walker <- function(
  seeds = NULL,
  hmin = 2,
  k = 12,
  alpha = 1,
  beta = 0.5,
  gamma = 0.5,
  delta = 2,
  eta = 1,
  crown_a = 0.35,
  crown_b = 0.65,
  eigen_radius = 1.0,
  density_radius = 1.0,
  probability_threshold = 0.5,
  max_iterations = 500,
  tolerance = 1e-5,
  crown_profile = c("beta", "parabolic", "cone"),
  verbose = FALSE
) {
  .assert_scalar_numeric(hmin, "hmin", lower = 0)
  .assert_scalar_integerish(k, "k", lower = 2)
  .assert_scalar_numeric(alpha, "alpha", lower = 0)
  .assert_scalar_numeric(beta, "beta", lower = 0)
  .assert_scalar_numeric(gamma, "gamma", lower = 0)
  .assert_scalar_numeric(delta, "delta", lower = 0)
  .assert_scalar_numeric(eta, "eta", lower = 0)
  .assert_scalar_numeric(crown_a, "crown_a", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(crown_b, "crown_b", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(eigen_radius, "eigen_radius", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(density_radius, "density_radius", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(probability_threshold, "probability_threshold", lower = 0)
  if (probability_threshold > 1) {
    stop("'probability_threshold' must be in [0, 1].", call. = FALSE)
  }
  .assert_scalar_integerish(max_iterations, "max_iterations", lower = 1)
  .assert_scalar_numeric(tolerance, "tolerance", lower = 0)
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE/FALSE.", call. = FALSE)
  }

  profile_id <- .resolve_profile_id(crown_profile)

  f <- function(las) {
    if (!methods::is(las, "LAS")) {
      stop("This algorithm expects a lidR LAS object.", call. = FALSE)
    }

    n <- lidR::npoints(las)
    if (n == 0L) {
      return(integer())
    }

    xyz <- .extract_xyz(las)
    seeds_df <- .coerce_seed_xyz(seeds, las, hmin = hmin)
    if (nrow(seeds_df) == 0L) {
      return(rep(NA_integer_, n))
    }

    if (isTRUE(verbose)) {
      message(sprintf("allometric_random_walker(): %d points, %d seeds", n, nrow(seeds_df)))
    }

    ids <- .Call(
      "cpp_allometric_random_walker",
      xyz$X,
      xyz$Y,
      xyz$Z,
      seeds_df$X,
      seeds_df$Y,
      seeds_df$Z,
      as.numeric(hmin),
      as.integer(k),
      as.numeric(alpha),
      as.numeric(beta),
      as.numeric(gamma),
      as.numeric(delta),
      as.numeric(eta),
      as.numeric(crown_a),
      as.numeric(crown_b),
      as.numeric(eigen_radius),
      as.numeric(density_radius),
      as.numeric(probability_threshold),
      as.integer(max_iterations),
      as.numeric(tolerance),
      as.integer(profile_id),
      PACKAGE = "spanner"
    )

    .check_algorithm_output(ids, n)
  }

  lidR::plugin_its(f, omp = TRUE, raster_based = FALSE)
}

#' Allometric Supervoxel Segmentation
#'
#' Builds a supervoxel-graph segmentation algorithm for
#' [lidR::segment_trees()]. Occupied voxels are segmented first and resulting
#' IDs are propagated back to points. The implementation includes small-island
#' cleanup and backfilling so isolated unlabeled speckles are less likely to be
#' left behind inside crowns.
#'
#' @inheritParams allometric_li_geodesic
#' @param voxel_size Voxel edge size in map units.
#' @param k_voxel Number of voxel neighbors in graph construction.
#' @param min_voxel_points Minimum points per voxel to keep as active.
#' @param component_gap Maximum centroid-gap for component continuity.
#' @param merge_threshold Maximum local voxel-edge coherence cost.
#' @param anisotropy_weight Weight for anisotropy differences.
#' @param verticality_weight Weight for verticality differences.
#' @param density_weight Weight for voxel density-gap penalties.
#' @param allometry_weight Weight for allometric penalties.
#'
#' @return A `PointCloudBased` `IndividualTreeSegmentation` algorithm object.
#' @export
#' @examples
#' LASfile <- system.file("extdata", "ALS_Clip.laz", package = "spanner")
#' las <- lidR::readLAS(LASfile, select = "xyz")
#' las <- lidR::decimate_points(las, lidR::random_per_voxel(0.33, 1L))
#' las <- lidR::classify_ground(las, lidR::ptd())
#' las <- lidR::normalize_height(las, lidR::tin())
#' las_sv <- lidR::segment_trees(
#'   las,
#'   allometric_supervoxel_segment(
#'     seeds = NULL,
#'     hmin = 2,
#'     voxel_size = 1.2,
#'     k_voxel = 18,
#'     min_voxel_points = 3,
#'     crown_a = 0.35,
#'     crown_b = 0.65,
#'     component_gap = 1.5,
#'     merge_threshold = 6,
#'     anisotropy_weight = 0.5,
#'     verticality_weight = 0.75,
#'     density_weight = 1,
#'     allometry_weight = 1.9,
#'     crown_profile = c("beta", "parabolic", "cone"),
#'     verbose = FALSE
#'   )
#' )
#' lidR::plot(las_sv, color = "treeID")
allometric_supervoxel_segment <- function(
  seeds = NULL,
  hmin = 2,
  voxel_size = 1.2,
  k_voxel = 18,
  min_voxel_points = 3,
  crown_a = 0.35,
  crown_b = 0.65,
  component_gap = 1.5,
  merge_threshold = 6,
  anisotropy_weight = 0.5,
  verticality_weight = 0.75,
  density_weight = 1,
  allometry_weight = 1.9,
  crown_profile = c("beta", "parabolic", "cone"),
  verbose = FALSE
) {
  .assert_scalar_numeric(hmin, "hmin", lower = 0)
  .assert_scalar_numeric(voxel_size, "voxel_size", lower = 0, lower_open = TRUE)
  .assert_scalar_integerish(k_voxel, "k_voxel", lower = 1)
  .assert_scalar_integerish(min_voxel_points, "min_voxel_points", lower = 1)
  .assert_scalar_numeric(crown_a, "crown_a", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(crown_b, "crown_b", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(component_gap, "component_gap", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(merge_threshold, "merge_threshold", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(anisotropy_weight, "anisotropy_weight", lower = 0)
  .assert_scalar_numeric(verticality_weight, "verticality_weight", lower = 0)
  .assert_scalar_numeric(density_weight, "density_weight", lower = 0)
  .assert_scalar_numeric(allometry_weight, "allometry_weight", lower = 0)
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE/FALSE.", call. = FALSE)
  }

  profile_id <- .resolve_profile_id(crown_profile)

  f <- function(las) {
    if (!methods::is(las, "LAS")) {
      stop("This algorithm expects a lidR LAS object.", call. = FALSE)
    }

    n <- lidR::npoints(las)
    if (n == 0L) {
      return(integer())
    }

    xyz <- .extract_xyz(las)
    seeds_df <- .coerce_seed_xyz(seeds, las, hmin = hmin)
    if (nrow(seeds_df) == 0L) {
      return(rep(NA_integer_, n))
    }

    if (isTRUE(verbose)) {
      message(sprintf("allometric_supervoxel_segment(): %d points, %d seeds", n, nrow(seeds_df)))
    }

    ids <- .Call(
      "cpp_allometric_supervoxel_segment",
      xyz$X,
      xyz$Y,
      xyz$Z,
      seeds_df$X,
      seeds_df$Y,
      seeds_df$Z,
      as.numeric(hmin),
      as.numeric(voxel_size),
      as.integer(k_voxel),
      as.integer(min_voxel_points),
      as.numeric(crown_a),
      as.numeric(crown_b),
      as.numeric(component_gap),
      as.numeric(merge_threshold),
      as.numeric(anisotropy_weight),
      as.numeric(verticality_weight),
      as.numeric(density_weight),
      as.numeric(allometry_weight),
      as.integer(profile_id),
      PACKAGE = "spanner"
    )

    .check_algorithm_output(ids, n)
  }

  lidR::plugin_its(f, omp = TRUE, raster_based = FALSE)
}
