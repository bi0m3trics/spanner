# Tests for the bole-segmentation pipeline:
#   classify_stem_points() → segment_bole() → compute_bole_volume()
#   → assess_tree_quality() → refit_trees()

# ============================================================================
# Shared test fixtures (built once for the entire file)
# ============================================================================

# Build a synthetic multi-tree LAS with known cylinder geometry
.make_bole_las <- function(n_trees = 2, n_pts_per_tree = 300,
                            radius = 0.15, tree_ht = 4.0, noise = 0.01) {
  set.seed(42)
  rows <- vector("list", n_trees)
  for (ti in seq_len(n_trees)) {
    cx <- ti * 3.0; cy <- 0.0
    theta <- runif(n_pts_per_tree, 0, 2 * pi)
    z     <- runif(n_pts_per_tree, 0.05, tree_ht)
    x     <- cx + radius * cos(theta) + rnorm(n_pts_per_tree, 0, noise)
    y     <- cy + radius * sin(theta) + rnorm(n_pts_per_tree, 0, noise)
    rows[[ti]] <- data.table::data.table(
      X = x, Y = y, Z = z,
      treeID = ti
    )
  }
  all_pts <- data.table::rbindlist(rows)
  # Add some ground points (treeID = 0, Z near 0)
  ground <- data.table::data.table(
    X = runif(50, 0, (n_trees + 1) * 3),
    Y = runif(50, -1, 1),
    Z = runif(50, 0, 0.04),
    treeID = 0L
  )
  all_pts <- data.table::rbindlist(list(all_pts, ground))
  las <- lidR::LAS(all_pts)
  las
}

# Build a minimal sf tree_locs from the synthetic LAS
.make_tree_locs <- function(n_trees = 2) {
  df <- data.frame(
    TreeID = seq_len(n_trees),
    X      = seq_len(n_trees) * 3.0,
    Y      = rep(0.0, n_trees),
    Radius = rep(0.15, n_trees),
    Height = rep(4.0, n_trees)
  )
  sf::st_as_sf(df, coords = c("X", "Y"), crs = 26912)
}


# ============================================================================
# classify_stem_points()
# ============================================================================

test_that("classify_stem_points returns a LAS with Stem and WoodProb columns", {
  las      <- .make_bole_las()
  tl       <- .make_tree_locs()

  result   <- classify_stem_points(las, tl, method = "eigen",
                                   z_min = 0.05, z_max = 4.0,
                                   search_radius = 0.5,
                                   neigh_k = 30L, ncpu = 1L)

  expect_true(inherits(result, "LAS"))
  expect_true("Stem" %in% names(result@data))
  expect_true("WoodProb" %in% names(result@data))
  expect_type(result@data$Stem, "logical")
  expect_true(is.numeric(result@data$WoodProb))
  expect_true(all(result@data$WoodProb >= 0 & result@data$WoodProb <= 1,
                  na.rm = TRUE))
  # At least some points should be classified as stem
  expect_true(sum(result@data$Stem, na.rm = TRUE) > 0)
})

test_that("classify_stem_points fails gracefully if treeID is missing", {
  las_no_id <- lidR::LAS(data.frame(X = 0, Y = 0, Z = 0))
  tl <- .make_tree_locs()
  expect_error(classify_stem_points(las_no_id, tl), "treeID")
})

test_that("classify_stem_points lewos method runs without error", {
  las <- .make_bole_las()
  tl  <- .make_tree_locs()
  result <- classify_stem_points(las, tl, method = "lewos",
                                 z_min = 0.05, z_max = 4.0,
                                 neigh_k = 30L, ncpu = 1L,
                                 n_propagation = 1L, k_propagation = 5L)
  expect_true("Stem" %in% names(result@data))
})


# ============================================================================
# segment_bole()
# ============================================================================

test_that("segment_bole returns a data.frame with expected columns", {
  las <- .make_bole_las()
  tl  <- .make_tree_locs()

  seg <- segment_bole(las, tl, algorithm = "ransac",
                      z_min = 0.05, z_max = 4.0,
                      dz = 0.5, overlap = 0.0,
                      use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 10L, n_best = 5L)

  expect_true(is.data.frame(seg))
  required <- c("TreeID", "Segment", "Z_low", "Z_high", "Z_mid",
                "Radius", "Diameter", "Valid", "Fit_method")
  for (col in required) {
    expect_true(col %in% names(seg),
                info = paste("Missing column:", col))
  }
  # Should have at least one valid segment per tree
  valid_per_tree <- tapply(seg$Valid, seg$TreeID, any, na.rm = TRUE)
  expect_true(all(valid_per_tree))
})

test_that("segment_bole radius estimates are close to ground truth", {
  las <- .make_bole_las(radius = 0.15, n_pts_per_tree = 500, noise = 0.005)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, algorithm = "ransac",
                      z_min = 0.05, z_max = 4.0,
                      dz = 1.0, overlap = 0.0,
                      use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 10L,
                      n_ransac = 20L, n_best = 10L)
  valid <- seg[!is.na(seg$Valid) & seg$Valid, ]
  expect_true(nrow(valid) > 0)
  median_r <- median(valid$Radius, na.rm = TRUE)
  expect_true(abs(median_r - 0.15) < 0.03,
              info = sprintf("Median radius %.4f, expected ~0.15", median_r))
})

test_that("segment_bole respects z_min and z_max", {
  las <- .make_bole_las()
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 1.0, z_max = 2.0, dz = 0.5,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 3L,
                      n_ransac = 5L, n_best = 3L)
  expect_true(all(seg$Z_mid >= 1.0 & seg$Z_mid <= 2.5))
})

test_that("segment_bole uses Stem column when present", {
  las <- .make_bole_las()
  tl  <- .make_tree_locs()
  # Add a Stem column that marks ALL points as non-stem → expect no valid segs
  las@data$Stem <- FALSE
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 0.5,
                      overlap = 0.0, use_stem_column = TRUE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 5L, n_best = 3L)
  expect_true(all(!seg$Valid))
})


# ============================================================================
# compute_bole_volume()
# ============================================================================

test_that("compute_bole_volume returns one row per tree", {
  las <- .make_bole_las(n_pts_per_tree = 400, noise = 0.005)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 0.5,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 10L, n_best = 5L)
  vol <- compute_bole_volume(seg, method = "cylinder")

  expect_equal(nrow(vol), 2L)
  expect_true("Volume_m3" %in% names(vol))
  expect_true(all(vol$Volume_m3 > 0, na.rm = TRUE))
})

test_that("compute_bole_volume smalian formula gives positive volumes", {
  las <- .make_bole_las(n_pts_per_tree = 400, noise = 0.005)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 1.0,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 10L, n_best = 5L)
  vol <- compute_bole_volume(seg, method = "smalian")
  valid_vols <- vol$Volume_m3[!is.na(vol$Volume_m3)]
  expect_true(all(valid_vols > 0))
})

test_that("compute_bole_volume cylinder volume is near analytical value", {
  # 2 trees × radius=0.15 × height=4 = pi*0.15^2*4 ≈ 0.2827 m3 each
  las <- .make_bole_las(n_pts_per_tree = 1000, noise = 0.003)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 0.25,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 30L, n_best = 10L)
  vol <- compute_bole_volume(seg, method = "cylinder")
  expect_true(all(!is.na(vol$Volume_m3)))
  # Allow ±30% tolerance given RANSAC noise and edge effects
  expect_true(all(abs(vol$Volume_m3 - 0.2827) < 0.09),
              info = paste("Volumes:", paste(round(vol$Volume_m3, 4),
                                             collapse = ", ")))
})


# ============================================================================
# assess_tree_quality()
# ============================================================================

test_that("assess_tree_quality returns an sf with quality_class column", {
  las <- .make_bole_las(n_pts_per_tree = 400, noise = 0.005)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 0.5,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 10L, n_best = 5L)
  qual <- assess_tree_quality(tl, seg)

  expect_true(inherits(qual, "sf"))
  expect_true("quality_class" %in% names(qual))
  expect_true("fit_score" %in% names(qual))
  expect_true(is.factor(qual$quality_class))
  expect_true(all(levels(qual$quality_class) %in% c("good", "marginal", "bad")))
})

test_that("assess_tree_quality fit_score is between 0 and 1", {
  las <- .make_bole_las(n_pts_per_tree = 400, noise = 0.005)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 0.5,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 10L, n_best = 5L)
  qual <- assess_tree_quality(tl, seg)
  scores <- qual$fit_score[!is.na(qual$fit_score)]
  expect_true(all(scores >= 0 & scores <= 1))
})

test_that("assess_tree_quality with vol_table NULL auto-computes volume", {
  las <- .make_bole_las(n_pts_per_tree = 200, noise = 0.005)
  tl  <- .make_tree_locs()
  seg <- segment_bole(las, tl, z_min = 0.05, z_max = 4.0, dz = 0.5,
                      overlap = 0.0, use_stem_column = FALSE,
                      search_radius = 0.5, min_pts = 5L,
                      n_ransac = 10L, n_best = 5L)
  expect_no_error(assess_tree_quality(tl, seg, vol_table = NULL))
})
