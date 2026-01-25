test_that("spanner_pal returns correct palette structure", {
  # Get the palette
  pal <- spanner_pal()
  
  # Test that it returns a character vector
  expect_type(pal, "character")
  
  # Test that it has 10 colors
  expect_equal(length(pal), 10)
  
  # Test that all colors are valid hex codes
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
  
  # Test that it has names
  expect_true(!is.null(names(pal)))
  expect_equal(length(names(pal)), 10)
})

test_that("las2xyz converts LAS to data frame correctly", {
  skip_if_not_installed("lidR")
  
  # Create a simple test LAS object
  LASfile <- system.file("extdata", "TLSSparseCloud_xyzOnly.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Convert to xyz
  xyz <- las2xyz(las)
  
  # Test that it returns a matrix
  expect_true(is.matrix(xyz))
  
  # Test that it has X, Y, Z columns
  expect_true(all(c("X", "Y", "Z") %in% colnames(xyz)))
  
  # Test that the number of rows matches
  expect_equal(nrow(xyz), length(las$X))
  
  # Test that values are numeric
  expect_type(xyz[,1], "double")
  expect_type(xyz[,2], "double")
  expect_type(xyz[,3], "double")
})
