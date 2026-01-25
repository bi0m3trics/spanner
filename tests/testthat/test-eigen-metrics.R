test_that("eigen_metrics basic functionality", {
  skip_if_not_installed("lidR")
  skip_on_cran()
  
  # Create a simple test LAS object
  LASfile <- system.file("extdata", "TLSSparseCloud_xyzOnly.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Subset to small area for speed
  las <- lidR::clip_rectangle(las, 
                               min(las$X), min(las$Y),
                               min(las$X) + 2, min(las$Y) + 2)
  
  # Calculate eigen metrics
  result <- eigen_metrics(las, radius = 0.5, ncpu = 1)
  
  # Test that result is a data.table
  expect_s3_class(result, "data.table")
  
  # Test that it has expected columns (at least some key ones)
  expected_cols <- c("Linearity", "Planarity", "Sphericity", "Verticality")
  expect_true(any(expected_cols %in% names(result)))
  
  # Test that number of rows matches input
  expect_equal(nrow(result), nrow(las@data))
  
  # Test that values are numeric
  for (col in names(result)) {
    if (col %in% expected_cols) {
      expect_type(result[[col]], "double")
    }
  }
})

test_that("eigen_metrics validates input", {
  skip_if_not_installed("lidR")
  skip_on_cran()
  
  # Test with invalid radius
  LASfile <- system.file("extdata", "TLSSparseCloud_xyzOnly.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  las <- lidR::clip_rectangle(las, 
                               min(las$X), min(las$Y),
                               min(las$X) + 1, min(las$Y) + 1)
  
  # Negative radius should error
  expect_error(eigen_metrics(las, radius = -1, ncpu = 1), "radius must be positive")
  
  # Zero radius should error
  expect_error(eigen_metrics(las, radius = 0, ncpu = 1), "radius must be positive")
  
  # Invalid ncpu should error
  expect_error(eigen_metrics(las, radius = 1, ncpu = -1), "ncpu must be a positive integer")
})
