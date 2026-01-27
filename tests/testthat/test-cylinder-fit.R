test_that("cylinderFit RANSAC returns expected structure", {
  skip_if_not_installed("lidR")
  skip_on_cran()
  
  # Create test point cloud data representing a cylinder
  # Generate points on a cylinder surface
  set.seed(123)
  n_points <- 100
  radius <- 0.2
  height <- 2
  
  theta <- runif(n_points, 0, 2 * pi)
  z <- runif(n_points, 0, height)
  x <- radius * cos(theta) + rnorm(n_points, 0, 0.01)
  y <- radius * sin(theta) + rnorm(n_points, 0, 0.01)
  
  xyz_df <- data.frame(X = x, Y = y, Z = z)
  las_obj <- lidR::LAS(xyz_df)
  
  # Test RANSAC fitting
  result <- cylinderFit(las_obj, method = "ransac", n_best = 10, n = 10, 
                        inliers = 0.8, conf = 0.95)
  
  # Test that result has expected components
  expect_type(result, "list")
  expect_s3_class(result, "data.frame")
  expect_true("radius" %in% names(result))
  
  # Test that radius is approximately correct (within reasonable tolerance)
  expect_true(abs(result$radius - radius) < 0.15)
  
  # Test that center coordinates exist
  expect_true(all(c("px", "py", "pz") %in% names(result)))
})

test_that("cylinderFit IRLS returns expected structure", {
  skip_if_not_installed("lidR")
  skip_on_cran()
  
  # Create test point cloud data
  set.seed(456)
  n_points <- 50
  radius <- 0.15
  height <- 1.5
  
  theta <- runif(n_points, 0, 2 * pi)
  z <- runif(n_points, 0, height)
  x <- radius * cos(theta) + rnorm(n_points, 0, 0.005)
  y <- radius * sin(theta) + rnorm(n_points, 0, 0.005)
  
  xyz_df <- data.frame(X = x, Y = y, Z = z)
  las_obj <- lidR::LAS(xyz_df)
  
  # Test IRLS fitting
  result <- cylinderFit(las_obj, method = "irls")
  
  # Test that result has expected components (actual column names from implementation)
  expect_type(result, "list")
  expect_true(all(c("rho", "theta", "phi", "alpha", "radius", "err") %in% names(result)))
  
  # Test that radius is positive
  expect_true(result$radius > 0)
  
  # Test that error exists
  expect_true("err" %in% names(result))
})

test_that("cylinderFit handles invalid input", {
  # Test with too few points
  xyz_df <- data.frame(X = c(1, 2), Y = c(1, 2), Z = c(1, 2))
  las_obj <- lidR::LAS(xyz_df)
  
  # With too few points, should return NULL or handle gracefully
  result <- cylinderFit(las_obj, method = "ransac")
  expect_true(is.null(result) || is.data.frame(result))
})
