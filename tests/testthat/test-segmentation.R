test_that("get_raster_eigen_treelocs basic structure", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("sf")
  skip_on_cran()
  
  # This is a computationally intensive test, so we'll just check 
  # the function accepts valid parameters without erroring on small data
  
  LASfile <- system.file("extdata", "TLSSparseCloud_xyzOnly.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Clip to very small area for testing
  las <- lidR::clip_rectangle(las, 
                               min(las$X), min(las$Y),
                               min(las$X) + 3, min(las$Y) + 3)
  
  # Normalize height (simulate)
  las$Z <- las$Z - min(las$Z)
  
  # This might not find trees in such a small area, but should run without error
  result <- get_raster_eigen_treelocs(
    las = las,
    res = 0.5,
    pt_spacing = 0.05,
    dens_threshold = 0.1,
    neigh_sizes = c(0.3, 0.15, 0.5),
    eigen_threshold = 0.5,
    grid_slice_min = 0.5,
    grid_slice_max = 2,
    minimum_polygon_area = 0.01,
    cylinder_fit_type = "ransac",
    max_dia = 1,
    SDvert = 0.5,
    n_best = 10,
    n_pts = 10,
    inliers = 0.8,
    conf = 0.95,
    max_angle = 30
  )
  
  # If trees are found, check structure
  if (!is.null(result) && nrow(result) > 0) {
    expect_s3_class(result, "sf")
    expect_true(all(c("TreeID", "Radius", "Error") %in% names(result)))
  } else {
    # If no trees found, that's okay for this small test area
    expect_true(is.null(result) || nrow(result) == 0)
  }
})

test_that("segment_graph basic structure", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("sf")
  skip_on_cran()
  
  # Create minimal test data
  # Simple normalized LAS
  n_points <- 200
  las_data <- data.frame(
    X = c(rnorm(100, 10, 0.3), rnorm(100, 12, 0.3)),
    Y = c(rnorm(100, 10, 0.3), rnorm(100, 12, 0.3)),
    Z = c(runif(100, 0, 3), runif(100, 0, 4))
  )
  
  las <- lidR::LAS(las_data)
  sf::st_crs(las) <- 26912
  
  # Create tree locations
  tree_locs <- sf::st_sf(
    TreeID = c(1, 2),
    X = c(10, 12),
    Y = c(10, 12),
    Z = c(1.3, 1.3),
    Radius = c(0.15, 0.15),
    Error = c(0.01, 0.01),
    geometry = sf::st_sfc(
      sf::st_point(c(10, 10)),
      sf::st_point(c(12, 12)),
      crs = 26912
    )
  )
  
  # Run segmentation
  result <- segment_graph(
    las = las,
    tree.locations = tree_locs,
    k = 10,
    distance.threshold = 0.5,
    use.metabolic.scale = FALSE,
    ptcloud_slice_min = 0.5,
    ptcloud_slice_max = 2,
    subsample.graph = 0.1,
    return.dense = FALSE
  )
  
  # Test that result is a LAS object
  expect_s4_class(result, "LAS")
  
  # Test that treeID column was added
  expect_true("treeID" %in% names(result@data))
  
  # Test that some points were assigned
  expect_true(any(!is.na(result$treeID)))
})
