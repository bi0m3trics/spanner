# Helper functions for tests

# Create a simple cylindrical point cloud for testing
create_test_cylinder <- function(radius = 0.2, height = 2, n_points = 100, noise = 0.01) {
  set.seed(123)
  theta <- runif(n_points, 0, 2 * pi)
  z <- runif(n_points, 0, height)
  x <- radius * cos(theta) + rnorm(n_points, 0, noise)
  y <- radius * sin(theta) + rnorm(n_points, 0, noise)
  
  data.frame(X = x, Y = y, Z = z)
}

# Create a simple test LAS object
create_test_las <- function(n_points = 100, extent = 10) {
  set.seed(456)
  las_data <- data.frame(
    X = runif(n_points, 0, extent),
    Y = runif(n_points, 0, extent),
    Z = runif(n_points, 0, 5)
  )
  
  lidR::LAS(las_data)
}

# Check if a value is within tolerance
expect_near <- function(actual, expected, tolerance = 0.1) {
  expect_true(abs(actual - expected) < tolerance,
              info = sprintf("Expected %f to be near %f (tolerance: %f)", 
                           actual, expected, tolerance))
}
