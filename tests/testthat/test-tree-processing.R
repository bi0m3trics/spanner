test_that("process_tree_data handles basic input", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("sf")
  skip_on_cran()
  
  # Create minimal test data
  # Simple sf object for tree locations
  treeData <- sf::st_sf(
    TreeID = c(1, 2),
    X = c(10, 15),
    Y = c(10, 15),
    Z = c(0, 0),
    Radius = c(0.2, 0.25),
    Error = c(0.01, 0.01),
    geometry = sf::st_sfc(
      sf::st_point(c(10, 10)),
      sf::st_point(c(15, 15)),
      crs = 26912
    )
  )
  
  # Create a simple LAS object with treeID
  n_points <- 100
  las_data <- data.frame(
    X = c(rnorm(50, 10, 0.5), rnorm(50, 15, 0.5)),
    Y = c(rnorm(50, 10, 0.5), rnorm(50, 15, 0.5)),
    Z = c(runif(50, 0, 5), runif(50, 0, 6)),
    treeID = c(rep(1, 50), rep(2, 50))
  )
  
  # Create LAS object
  segmentedLAS <- lidR::LAS(las_data)
  sf::st_crs(segmentedLAS) <- 26912
  
  # Suppress messages during testing
  result <- suppressMessages(
    process_tree_data(treeData, segmentedLAS, return_sf = FALSE)
  )
  
  # Test that result is an sf object
  expect_s3_class(result, "sf")
  
  # Test that it has expected columns
  expect_true(all(c("TreeID", "height", "crown_area", "diameter", 
                    "crown_base_height", "crown_volume") %in% names(result)))
  
  # Test that diameter is calculated correctly
  expect_equal(result$diameter[1], treeData$Radius[1] * 2)
  
  # Test that height values are not NA
  expect_false(any(is.na(result$height)))
  
  # Test return_sf = TRUE
  result_sf <- suppressMessages(
    process_tree_data(treeData, segmentedLAS, return_sf = TRUE)
  )
  expect_s3_class(result_sf, "sf")
})

test_that("estimate_crown_base_height calculates reasonable values", {
  # Test with simple z values
  z <- seq(0, 10, length.out = 100)
  
  cbh <- spanner:::estimate_crown_base_height(z, threshold = 0.05)
  
  # Test that result is numeric
  expect_type(cbh, "double")
  
  # Test that result is within range of z values
  expect_true(cbh >= min(z) && cbh <= max(z))
})
