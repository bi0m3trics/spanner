test_that("package loads without errors", {
  expect_true("spanner" %in% loadedNamespaces())
})

test_that("package has DESCRIPTION file", {
  desc_file <- system.file("DESCRIPTION", package = "spanner")
  expect_true(file.exists(desc_file))
})

test_that("exported functions are available", {
  # Check main exported functions
  expect_true(exists("get_raster_eigen_treelocs"))
  expect_true(exists("segment_graph"))
  expect_true(exists("process_tree_data"))
  expect_true(exists("process_rasters_patchmorph"))
  expect_true(exists("sum_rasters_by_suitability"))
  expect_true(exists("plot_raster_by_name"))
  expect_true(exists("eigen_metrics"))
  expect_true(exists("cylinderFit"))
  expect_true(exists("las2xyz"))
  expect_true(exists("spanner_pal"))
})

test_that("package data files exist", {
  # Check for example data files
  tls_file <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
  mls_file <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  
  expect_true(file.exists(tls_file))
  expect_true(file.exists(mls_file))
})

test_that("package has proper version", {
  desc <- packageDescription("spanner")
  expect_true(!is.null(desc$Version))
  expect_match(desc$Version, "^\\d+\\.\\d+\\.\\d+")
})


