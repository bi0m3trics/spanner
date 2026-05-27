# =============================================================================
# spanner – LAI / PAD / eigen_metrics / transmittance test script
#
# Tests:
#   1. Load LAS and inspect
#   2. compute_lai()                  – MacArthur-Horn LAI/LAD raster (LAS)
#   3. compute_lai()                  – MacArthur-Horn LAI/LAD raster (LAScatalog)
#   4. compute_pad_voxels()           – Beer-Lambert 3-D PAD voxel grid
#   5. eigen_metrics()                – per-point eigen decomposition metrics
#   6. compute_transmittance_raster() – column transmittance (τ) + vostokR correction
#
# Input: Strawberry_20260428_452000_3809000_fixed_preprocessed.laz
#        (assumed height-normalized; Z >= 0 above ground)
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
devtools::load_all("C:/GitHubRepos/spanner")  # load source directly (no rebuild needed)
library(lidR)
library(data.table)
library(terra)

LASfile <- "C:/GitHubRepos/spanner/ignoreme/Strawberry_20260428_452000_3809000_fixed_preprocessed.laz"
stopifnot(file.exists(LASfile))

lidR::set_lidr_threads(8L)

# =============================================================================
# 1. Load and inspect
# =============================================================================
cat("\n=== 1. Loading LAS ===\n")
las <- lidR::readLAS(LASfile, select = "xyzrc")  # r = ReturnNumber (required by compute_lai)
cat("Points     :", nrow(las@data), "\n")
cat("Z range    :", round(min(las@data$Z), 2), "–", round(max(las@data$Z), 2), "m\n")
cat("CRS        :", sf::st_crs(las)$input, "\n\n")

# Confirm the cloud is already height-normalized (Z ≥ 0).
# If not, uncomment:
# las <- lidR::normalize_height(las, lidR::tin())
# las <- lidR::filter_poi(las, Z >= 0)

# =============================================================================
# 2. compute_lai() – LAS input
# =============================================================================
cat("=== 2. compute_lai() — LAS ===\n")

lai_las <- compute_lai(las, res_xy = 10, dz = 1, k = 0.5, min_height = 0)

cat("Class      :", class(lai_las), "\n")
cat("Layers     :", names(lai_las), "\n")
cat("Dimensions :", paste(dim(lai_las), collapse = " × "), "\n")
cat("LAI  – min/mean/max:",
    round(terra::global(lai_las[["LAI"]],     "min",  na.rm = TRUE)[[1]], 3), "/",
    round(terra::global(lai_las[["LAI"]],     "mean", na.rm = TRUE)[[1]], 3), "/",
    round(terra::global(lai_las[["LAI"]],     "max",  na.rm = TRUE)[[1]], 3), "\n")
cat("LAD_mean – min/mean/max:",
    round(terra::global(lai_las[["LAD_mean"]], "min",  na.rm = TRUE)[[1]], 3), "/",
    round(terra::global(lai_las[["LAD_mean"]], "mean", na.rm = TRUE)[[1]], 3), "/",
    round(terra::global(lai_las[["LAD_mean"]], "max",  na.rm = TRUE)[[1]], 3), "\n\n")

terra::plot(lai_las[["LAI"]],     main = "LAI (MacArthur-Horn) – LAS")
terra::plot(lai_las[["LAD_mean"]], main = "Mean LAD – LAS")

# =============================================================================
# 3. compute_lai() – LAScatalog input
# =============================================================================
cat("=== 3. compute_lai() — LAScatalog ===\n")

ctg <- lidR::readLAScatalog(LASfile)
lidR::opt_chunk_size(ctg)   <- 200          # 200 m tiles
lidR::opt_chunk_buffer(ctg) <- 10           # 10 m buffer (>= res_xy)
lidR::opt_progress(ctg)     <- TRUE

# Suppress writing output to disk (in-memory merge via catalog_apply)
lidR::opt_output_files(ctg) <- ""

# future::plan(future::multisession, workers = 4)  # uncomment for parallel chunks

lai_ctg <- compute_lai(ctg, res_xy = 10, dz = 1, k = 0.5, min_height = 0)

cat("Class      :", class(lai_ctg), "\n")
cat("Layers     :", names(lai_ctg), "\n")
cat("Dimensions :", paste(dim(lai_ctg), collapse = " × "), "\n")
cat("LAI  – min/mean/max:",
    round(terra::global(lai_ctg[["LAI"]],     "min",  na.rm = TRUE)[[1]], 3), "/",
    round(terra::global(lai_ctg[["LAI"]],     "mean", na.rm = TRUE)[[1]], 3), "/",
    round(terra::global(lai_ctg[["LAI"]],     "max",  na.rm = TRUE)[[1]], 3), "\n\n")

terra::plot(lai_ctg[["LAI"]], main = "LAI (MacArthur-Horn) – LAScatalog")

# Sanity check: LAS and catalog results should be essentially identical.
# Small differences are expected at tile edges due to buffer cropping.
diff_r <- terra::resample(lai_ctg[["LAI"]], lai_las[["LAI"]]) - lai_las[["LAI"]]
cat("Max absolute LAI diff (LAS vs catalog):",
    round(terra::global(abs(diff_r), "max", na.rm = TRUE)[[1]], 4), "\n\n")

# =============================================================================
# 4. compute_pad_voxels()
# =============================================================================
cat("=== 4. compute_pad_voxels() ===\n")

pad <- compute_pad_voxels(las,
                          vox_size         = 0.5,
                          G                = 0.5,
                          mean_path_factor = 0.843,
                          min_height       = 0,
                          max_occlusion    = 0.8,
                          min_pad          = 0.01,
                          max_pad          = 6.0,
                          ncpu             = 4L)

cat("Class      :", class(pad), "\n")
cat("Columns    :", paste(names(pad), collapse = ", "), "\n")
cat("Voxel rows :", nrow(pad), "\n")
cat("\nVoxel class breakdown:\n")
print(pad[, .N, by = VoxelClass][order(-N)])

cat("\nPAD (foliage voxels only) – min/mean/max:\n")
fol <- pad[VoxelClass == "foliage", PAD]
cat(" ", round(min(fol),  3), "/",
        round(mean(fol), 3), "/",
        round(max(fol),  3), "\n")

# Column LAI from PAD voxels
lai_pad <- pad[VoxelClass == "foliage",
               .(LAI_pad = sum(PAD * 0.5, na.rm = TRUE)),
               by = .(X, Y)]
cat("\nColumn LAI (PAD-based) – min/mean/max:\n")
cat(" ", round(min(lai_pad$LAI_pad),  3), "/",
        round(mean(lai_pad$LAI_pad), 3), "/",
        round(max(lai_pad$LAI_pad),  3), "\n")

# Vertical PAD profile
profile <- pad[, .(PAD_mean     = mean(PAD[VoxelClass == "foliage"], na.rm = TRUE),
                   foliage_frac = mean(VoxelClass == "foliage"),
                   empty_frac   = mean(VoxelClass == "empty")),
               by = .(HAG_bin = floor(Z / 0.5) * 0.5)]
cat("\nVertical profile (top 10 rows by height):\n")
print(profile[order(-HAG_bin)][1:min(10, .N)])

# =============================================================================
# 5. eigen_metrics()
# =============================================================================
cat("\n=== 5. eigen_metrics() ===\n")

# Thin the cloud first if it is very dense (optional but speeds things up).
las_thin <- lidR::decimate_points(las, lidR::random_per_voxel(res = 0.05, n = 1L))
cat("Thinned point count:", nrow(las_thin@data), "\n")

# kNN neighbourhood
cat("\n-- kNN (k = 20) --\n")
em_k <- eigen_metrics(las_thin, k = 20L, ncpu = 8L)
cat("Columns :", paste(names(em_k), collapse = ", "), "\n")
cat("Rows    :", nrow(em_k), "\n")

# Check E1x/E1y/E1z are present
stopifnot(all(c("E1x", "E1y", "E1z") %in% names(em_k)))
cat("E1x/E1y/E1z present: TRUE\n")
cat("E1z summary (should be near ±1 for vertical structures):\n")
print(summary(em_k$E1z))

# Spherical neighbourhood
cat("\n-- Spherical (r = 0.1 m) --\n")
em_r <- eigen_metrics(las_thin, r = 0.1, ncpu = 8L)
cat("Linearity  – mean:", round(mean(em_r$Linearity,   na.rm = TRUE), 4), "\n")
cat("Verticality– mean:", round(mean(em_r$Verticality, na.rm = TRUE), 4), "\n")
cat("Curvature  – mean:", round(mean(em_r$Curvature,   na.rm = TRUE), 4), "\n")

# kNN + radius cap
cat("\n-- kNN with radius cap (r = 0.2, k = 15) --\n")
em_rk <- eigen_metrics(las_thin, r = 0.2, k = 15L, ncpu = 8L)
cat("Rows         :", nrow(em_rk), "\n")
cat("NumNeighbors – mean:", round(mean(em_rk$NumNeighbors, na.rm = TRUE), 1), "\n")

# =============================================================================
# 6. compute_transmittance_raster() – canopy τ and vostokR correction
# =============================================================================
cat("\n=== 6. compute_transmittance_raster() ===\n")

tau_r <- compute_transmittance_raster(pad,
                                      G        = 0.5,
                                      vox_size = 0.5,
                                      max_pad  = 6.0,   # same as compute_pad_voxels()
                                      crs      = sf::st_crs(las)$wkt)

cat("Class      :", class(tau_r), "\n")
cat("Layer name :", names(tau_r), "\n")
cat("Resolution :", terra::res(tau_r), "m\n")
cat("tau – min/mean/max:",
    round(terra::global(tau_r, "min",  na.rm = TRUE)[[1]], 4), "/",
    round(terra::global(tau_r, "mean", na.rm = TRUE)[[1]], 4), "/",
    round(terra::global(tau_r, "max",  na.rm = TRUE)[[1]], 4), "\n\n")

terra::plot(tau_r, main = "Column transmittance (tau) – ALS Beer-Lambert",
            range = c(0, 1))

# ---- Correct a vostokR solar_potential raster (if available) ---------------
# vostokR::calculate_solar_potential_cpp() returns Wh m⁻² day⁻¹ per point.
# After converting to a SpatRaster (solar_potential), multiply by tau:
#
#   tau_rs      <- terra::resample(tau_r, solar_potential, method = "bilinear")
#   solar_floor <- solar_potential * tau_rs
#   terra::plot(solar_floor, main = "Solar potential at forest floor (Wh/m²/day)")
#
# If solar_potential and tau_r share the same grid (same origin + resolution),
# no resample is needed — multiply directly.

cat("=== Section 6 complete ===\n")

cat("\n=== All tests complete ===\n")
