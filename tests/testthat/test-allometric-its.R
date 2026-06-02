read_mixed_conifer <- function() {
  LASfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
  skip_if(!nzchar(LASfile), "lidR extdata MixedConifer.laz is unavailable")

  poi <- "-drop_z_below 0 -inside 481280 3812940 481320 3812980"
  las <- lidR::readLAS(LASfile, select = "xyz", filter = poi)
  skip_if(lidR::is.empty(las), "MixedConifer subset is empty")
  las
}

test_that("mixed conifer fixture is readable", {
  skip_if_not_installed("lidR")

  las <- read_mixed_conifer()
  expect_s4_class(las, "LAS")
  expect_gt(lidR::npoints(las), 0)
})

test_that("algorithms return valid integer IDs", {
  skip_if_not_installed("lidR")

  las <- read_mixed_conifer()
  n <- lidR::npoints(las)

  algos <- list(
    spanner::allometric_li_geodesic(hmin = 2, k = 8, geodesic_threshold = 8),
    spanner::allometric_random_walker(hmin = 2, k = 8, max_iterations = 80, probability_threshold = 0.4),
    spanner::allometric_supervoxel_segment(hmin = 2, voxel_size = 0.75, merge_threshold = 6)
  )

  for (algo in algos) {
    ids <- algo(las)
    expect_type(ids, "integer")
    expect_length(ids, n)
    expect_true(all(is.na(ids) | ids > 0L))
  }
})

test_that("algorithms are accepted by segment_trees", {
  skip_if_not_installed("lidR")

  las <- read_mixed_conifer()

  algos <- list(
    spanner::allometric_li_geodesic(hmin = 2, k = 8, geodesic_threshold = 8),
    spanner::allometric_random_walker(hmin = 2, k = 8, max_iterations = 60, probability_threshold = 0.35),
    spanner::allometric_supervoxel_segment(hmin = 2, voxel_size = 0.8, merge_threshold = 6)
  )

  for (algo in algos) {
    expect_true(methods::is(algo, "PointCloudBased"))
    expect_true(methods::is(algo, "IndividualTreeSegmentation"))

    seg <- lidR::segment_trees(las, algo)
    expect_true("treeID" %in% names(seg@data))
    expect_length(seg@data$treeID, lidR::npoints(seg))
  }
})
