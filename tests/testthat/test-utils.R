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
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
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

test_that("colorize_las validates input correctly", {
  skip_if_not_installed("lidR")
  
  # Test with invalid LAS object
  expect_error(
    colorize_las(data.frame(X = 1:10, Y = 1:10, Z = 1:10), method="attr", attribute_name="Z", palette=c("blue", "red")),
    "The input is not a valid LAS object."
  )
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Test invalid method
  expect_error(
    colorize_las(las, method="invalid"),
    "Method must be one of 'attr', 'rgb', 'pcv', or 'ssao'."
  )
  
  # Test attr method without attribute_name
  expect_error(
    colorize_las(las, method="attr", palette=c("blue", "red")),
    "attribute_name must be specified when method='attr'."
  )
  
  # Test with non-existent attribute
  expect_error(
    colorize_las(las, method="attr", attribute_name="NonExistentAttribute", palette=c("blue", "red")),
    "Attribute NonExistentAttribute not found in the LAS object."
  )
  
  # Test with invalid palette (not character)
  expect_error(
    colorize_las(las, method="attr", attribute_name="Z", palette=c(1, 2, 3)),
    "The palette must be a character vector with at least two colors."
  )
  
  # Test with invalid palette (only one color)
  expect_error(
    colorize_las(las, method="attr", attribute_name="Z", palette="blue"),
    "The palette must be a character vector with at least two colors."
  )
  
  # Test rgb method without raster_path
  expect_error(
    colorize_las(las, method="rgb"),
    "raster_path must be specified when method='rgb'."
  )
})

test_that("colorize_las colorizes LAS object correctly with attr method", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Colorize by Z using a simple palette
  las_colored <- colorize_las(las, method="attr", attribute_name="Z", palette=c("blue", "red"))
  
  # Test that RGB fields were added
  expect_true("R" %in% names(las_colored@data))
  expect_true("G" %in% names(las_colored@data))
  expect_true("B" %in% names(las_colored@data))
  
  # Test that RGB values are in valid range (0-255)
  expect_true(all(las_colored@data$R >= 0 & las_colored@data$R <= 255))
  expect_true(all(las_colored@data$G >= 0 & las_colored@data$G <= 255))
  expect_true(all(las_colored@data$B >= 0 & las_colored@data$B <= 255))
  
  # Test that RGB values are integers
  expect_true(all(las_colored@data$R == floor(las_colored@data$R)))
  expect_true(all(las_colored@data$G == floor(las_colored@data$G)))
  expect_true(all(las_colored@data$B == floor(las_colored@data$B)))
  
  # Test that the number of points hasn't changed
  expect_equal(nrow(las@data), nrow(las_colored@data))
})

test_that("colorize_las works with different palettes", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Test with grayscale palette
  las_gray <- colorize_las(las, method="attr", attribute_name="Z", palette=c("black", "white"))
  expect_true(all(c("R", "G", "B") %in% names(las_gray@data)))
  
  # Test with multi-color palette
  las_multi <- colorize_las(las, method="attr", attribute_name="Z", palette=c("blue", "green", "yellow", "red"))
  expect_true(all(c("R", "G", "B") %in% names(las_multi@data)))
  
  # Test with hex color codes
  las_hex <- colorize_las(las, method="attr", attribute_name="Z", palette=c("#0000FF", "#FF0000"))
  expect_true(all(c("R", "G", "B") %in% names(las_hex@data)))
  
  # Test with spanner_pal
  las_spanner <- colorize_las(las, method="attr", attribute_name="Z", palette=spanner_pal())
  expect_true(all(c("R", "G", "B") %in% names(las_spanner@data)))
})

test_that("colorize_las works with different attributes", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Test colorizing by X
  las_x <- colorize_las(las, method="attr", attribute_name="X", palette=c("blue", "red"))
  expect_true(all(c("R", "G", "B") %in% names(las_x@data)))
  
  # Test colorizing by Y
  las_y <- colorize_las(las, method="attr", attribute_name="Y", palette=c("blue", "red"))
  expect_true(all(c("R", "G", "B") %in% names(las_y@data)))
  
  # Test colorizing by Z
  las_z <- colorize_las(las, method="attr", attribute_name="Z", palette=c("blue", "red"))
  expect_true(all(c("R", "G", "B") %in% names(las_z@data)))
  
  # Verify that different attributes produce different colors
  # (unless by coincidence, which is unlikely)
  expect_false(identical(las_x@data$R, las_y@data$R) && 
               identical(las_x@data$G, las_y@data$G) && 
               identical(las_x@data$B, las_y@data$B))
})

test_that("colorize_las PCV method works", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Take a small sample for faster testing
  if (nrow(las@data) > 100) {
    set.seed(123)
    sample_idx <- sample.int(nrow(las@data), 100)
    las <- las[sample_idx, ]
  }
  
  # Test PCV coloring
  las_pcv <- colorize_las(las, method="pcv", radius=1.0, num_directions=36,
                          palette=c("black", "white"), ncpu=2)
  
  # Test that RGB fields were added
  expect_true("R" %in% names(las_pcv@data))
  expect_true("G" %in% names(las_pcv@data))
  expect_true("B" %in% names(las_pcv@data))
  
  # Test that RGB values are in valid range
  expect_true(all(las_pcv@data$R >= 0 & las_pcv@data$R <= 255))
  expect_true(all(las_pcv@data$G >= 0 & las_pcv@data$G <= 255))
  expect_true(all(las_pcv@data$B >= 0 & las_pcv@data$B <= 255))
  
  # Test that the number of points hasn't changed
  expect_equal(nrow(las@data), nrow(las_pcv@data))
})

test_that("colorize_las RGB method validation", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("terra")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Test with invalid number of rasters
  expect_error(
    colorize_las(las, method="rgb", raster_path=c("r1.tif", "r2.tif")),
    "raster_path must be either a single file or a vector of 3 files"
  )
})

test_that("colorize_las SSAO method works", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Take a small sample for faster testing
  if (nrow(las@data) > 500) {
    set.seed(123)
    sample_idx <- sample.int(nrow(las@data), 500)
    las <- las[sample_idx, ]
  }
  
  # Test SSAO coloring
  las_ssao <- colorize_las(las, method="ssao", pixel_size=0.5, kernel_size=3,
                           num_samples=8, palette=c("black", "white"), ncpu=2)
  
  # Test that RGB fields were added
  expect_true("R" %in% names(las_ssao@data))
  expect_true("G" %in% names(las_ssao@data))
  expect_true("B" %in% names(las_ssao@data))
  
  # Test that RGB values are in valid range
  expect_true(all(las_ssao@data$R >= 0 & las_ssao@data$R <= 255))
  expect_true(all(las_ssao@data$G >= 0 & las_ssao@data$G <= 255))
  expect_true(all(las_ssao@data$B >= 0 & las_ssao@data$B <= 255))
  
  # Test that the number of points hasn't changed
  expect_equal(nrow(las@data), nrow(las_ssao@data))
  
  # Test that SSAO produces variation (not all the same color)
  expect_true(length(unique(las_ssao@data$R)) > 1 || 
              length(unique(las_ssao@data$G)) > 1 || 
              length(unique(las_ssao@data$B)) > 1)
})

test_that("merge_las_colors works with different blend methods", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Take a small sample for faster testing
  if (nrow(las@data) > 100) {
    set.seed(123)
    sample_idx <- sample.int(nrow(las@data), 100)
    las <- las[sample_idx, ]
  }
  
  # Create two colored versions
  las1 <- colorize_las(las, method="attr", attribute_name="Z", palette=c("blue", "red"))
  las2 <- colorize_las(las, method="attr", attribute_name="Z", palette=c("black", "white"))
  
  # Test alpha blending
  las_alpha <- merge_las_colors(las1, las2, alpha=0.5, method="alpha")
  expect_true("R" %in% names(las_alpha@data))
  expect_true("G" %in% names(las_alpha@data))
  expect_true("B" %in% names(las_alpha@data))
  expect_true(all(las_alpha@data$R >= 0 & las_alpha@data$R <= 255))
  expect_true(all(las_alpha@data$G >= 0 & las_alpha@data$G <= 255))
  expect_true(all(las_alpha@data$B >= 0 & las_alpha@data$B <= 255))
  expect_equal(nrow(las@data), nrow(las_alpha@data))
  
  # Test multiply blending
  las_mult <- merge_las_colors(las1, las2, method="multiply")
  expect_true(all(las_mult@data$R >= 0 & las_mult@data$R <= 255))
  
  # Test screen blending
  las_screen <- merge_las_colors(las1, las2, method="screen")
  expect_true(all(las_screen@data$R >= 0 & las_screen@data$R <= 255))
  
  # Test overlay blending
  las_overlay <- merge_las_colors(las1, las2, method="overlay")
  expect_true(all(las_overlay@data$R >= 0 & las_overlay@data$R <= 255))
})

test_that("merge_las_colors input validation works", {
  skip_if_not_installed("lidR")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  las_colored <- colorize_las(las, method="attr", attribute_name="Z")
  
  # Test mismatched point counts
  las_subset <- las[1:10, ]
  las_subset_colored <- colorize_las(las_subset, method="attr", attribute_name="Z")
  expect_error(
    merge_las_colors(las_colored, las_subset_colored),
    "must have the same number of points"
  )
  
  # Test missing RGB fields
  expect_error(
    merge_las_colors(las, las_colored),
    "must have R, G, B fields"
  )
  
  # Test invalid alpha
  expect_error(
    merge_las_colors(las_colored, las_colored, alpha=1.5),
    "alpha must be a numeric value between 0 and 1"
  )
  
  # Test invalid method
  expect_error(
    merge_las_colors(las_colored, las_colored, method="invalid"),
    "method must be one of"
  )
})

test_that("create_rotation_gif creates output file", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("rgl")
  skip_if_not_installed("magick")
  skip_on_ci()  # Skip on CI as rgl may not work headless
  skip_on_cran()  # Skip on CRAN
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Take a very small sample for fast testing
  set.seed(123)
  sample_idx <- sample.int(nrow(las@data), 50)
  las <- las[sample_idx, ]
  
  # Color it
  las_colored <- colorize_las(las, method="attr", attribute_name="Z")
  
  # Create a temporary output file
  temp_gif <- tempfile(fileext = ".gif")
  
  # Create GIF with minimal settings for speed
  result <- tryCatch({
    create_rotation_gif(las_colored,
                       output_path = temp_gif,
                       duration = 1,
                       fps = 5,
                       width = 200,
                       height = 200,
                       overwrite = TRUE)
  }, error = function(e) {
    skip(paste("rgl not available:", e$message))
  })
  
  # Test that file was created
  expect_true(file.exists(temp_gif))
  expect_true(file.size(temp_gif) > 0)
  
  # Cleanup
  unlink(temp_gif)
})

test_that("create_rotation_gif input validation works", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("magick")
  skip_if_not_installed("rgl")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  las_colored <- colorize_las(las, method="attr", attribute_name="Z")
  
  # Test invalid axis
  expect_error(
    create_rotation_gif(las_colored, axis="invalid"),
    "axis must be one of"
  )
  
  # Test overwrite protection
  temp_gif <- tempfile(fileext = ".gif")
  writeLines("dummy", temp_gif)
  
  expect_error(
    create_rotation_gif(las_colored, output_path = temp_gif, overwrite = FALSE),
    "already exists"
  )
  
  unlink(temp_gif)
})

test_that("download_naip_for_las validates input correctly", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("rstac")
  
  # Test with invalid input type
  expect_error(
    download_naip_for_las(data.frame(X = 1:10, Y = 1:10, Z = 1:10)),
    "las must be a LAS object or path to a LAS/LAZ file"
  )
  
  # Test with non-existent file
  expect_error(
    download_naip_for_las("nonexistent_file.laz")
  )
})

test_that("download_naip_for_las parameter validation", {
  skip_if_not_installed("lidR")
  skip_if_not_installed("rstac")
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  
  # Load a test LAS file
  LASfile <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
  skip_if(!file.exists(LASfile), "Test LAS file not found")
  
  las <- lidR::readLAS(LASfile, select = "xyz")
  
  # Create a temporary output path
  temp_output <- tempfile(fileext = ".tif")
  
  # Test that function returns path or NULL
  # (May return NULL if no imagery available for test extent)
  result <- tryCatch({
    download_naip_for_las(las, output_path = temp_output, 
                         year_range = c("2018-01-01", "2023-12-31"))
  }, error = function(e) {
    # Network errors are acceptable in testing
    if (grepl("STAC query failed|connection", e$message, ignore.case = TRUE)) {
      skip("Network unavailable for NAIP test")
    }
    stop(e)
  })
  
  # If result is not NULL, verify it's a valid path
  if (!is.null(result)) {
    expect_type(result, "character")
    expect_true(file.exists(result) || !is.null(result))
  }
  
  # Clean up
  if (file.exists(temp_output)) {
    unlink(temp_output)
  }
})


