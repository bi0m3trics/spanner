test_that("eigen_metrics basic functionality", {
  skip_if_not_installed("lidR")
  skip_on_cran()
  
  # Create a simple test LAS object
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Subset to small area for speed
  las <- lidR::clip_rectangle(las, 
                               min(las$X), min(las$Y),
                               min(las$X) + 2, min(las$Y) + 2)
  
  expected_cols <- c("Linearity", "Planarity", "Sphericity", "Verticality")

  # --- sphere-only mode ---
  result_r <- eigen_metrics(las, r = 0.5, ncpu = 1)
  expect_s3_class(result_r, "data.table")
  expect_true(any(expected_cols %in% names(result_r)))
  expect_equal(nrow(result_r), nrow(las@data))
  for (col in intersect(names(result_r), expected_cols))
    expect_type(result_r[[col]], "double")

  # --- knn-only mode ---
  result_k <- eigen_metrics(las, k = 10L, ncpu = 1)
  expect_s3_class(result_k, "data.table")
  expect_equal(nrow(result_k), nrow(las@data))

  # --- knn + radius cap mode ---
  result_rk <- eigen_metrics(las, r = 0.5, k = 10L, ncpu = 1)
  expect_s3_class(result_rk, "data.table")
  expect_equal(nrow(result_rk), nrow(las@data))
})

test_that("eigen_metrics validates input", {
  skip_if_not_installed("lidR")
  skip_on_cran()
  
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  las <- lidR::clip_rectangle(las, 
                               min(las$X), min(las$Y),
                               min(las$X) + 1, min(las$Y) + 1)
  
  # Must supply at least one of r or k
  expect_error(eigen_metrics(las, ncpu = 1), "at least one of")

  # Negative radius should error
  expect_error(eigen_metrics(las, r = -1, ncpu = 1), "r must be positive")
  
  # Zero radius should error
  expect_error(eigen_metrics(las, r = 0, ncpu = 1), "r must be positive")

  # Non-positive k should error
  expect_error(eigen_metrics(las, k = 0L, ncpu = 1), "k must be a positive integer")
  
  # Invalid ncpu should error
  expect_error(eigen_metrics(las, r = 1, ncpu = -1), "ncpu must be a positive integer")
})


