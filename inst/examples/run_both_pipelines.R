# =============================================================================
# spanner – Full Pipeline Example Script
#
# Demonstrates and tests both pipelines using the bundled TLS example data.
#
# Pipeline 1 (original):
#   get_raster_eigen_treelocs() → segment_graph() → process_tree_data()
#
# Pipeline 2 (new bole segmentation):
#   get_raster_eigen_treelocs() → segment_graph()
#     → classify_lw_points()    (bole + branch + other in one eigen pass)
#     → segment_bole()
#     → compute_bole_volume()
#     → assess_tree_quality()
#     → refit_trees()          [optional – only runs if bad/marginal trees exist]
#
# The two pipelines share the tree-detection and segmentation steps, so you
# can run one or both from the same segmented LAS.
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------
# Install/load spanner (build from source if working from the repo):
# devtools::install("C:/GitHubRepos/spanner")

library(spanner)
library(lidR)
library(sf)

# Optional: parallel processing – adjust to your machine
lidR::set_lidr_threads(15L)

# lowell.laz: 25.8 M points, 43×43 m, Z –0.43→27.8 m, 13 771 pts/m²
# No embedded CRS; ReturnNumber = 0 (single-return scanner).
# Use a 20-mm voxel filter to reduce to ~2 500 pts/m² before processing.
LASfile <- "c:/GitHubRepos/spanner/ignoreme/lowell.laz"
stopifnot(file.exists(LASfile))
las     <- readTLSLAS(LASfile, select = "xyzcr", "-filter_with_voxel 0.0254")
sf::st_crs(las) <- 32612   # UTM Zone 12N

# Normalize height above ground and remove noise
las <- classify_ground(las, lidR::ptd())
las <- normalize_height(las, lidR::tin())
las <- classify_noise(las, ivf(0.33,3))
las <- remove_noise(las)
las <- lidR::filter_poi(las, Z >= 0)   # drop below-ground artefacts

las_label <- "lowell.laz"

cat("Dataset:", las_label, "\n")
cat("LAS summary:\n")
print(las)
cat("Point count:", nrow(las@data), "\n")
cat("Z range:", round(min(las@data$Z), 2), "–", round(max(las@data$Z), 2), "m\n\n")

# =============================================================================
# 2. Tree detection  (shared by both pipelines)
# =============================================================================
cat("--- Step 1: Tree detection (get_raster_eigen_treelocs) ---\n")

tree_locs <- get_raster_eigen_treelocs(
  las                  = las,
  res                  = 0.025,
  pt_spacing           = 0.0254,
  dens_threshold       = 0.10,
  neigh_sizes          = c(0.25, 0.15, 0.66),
  eigen_threshold      = 0.65,
  grid_slice_min       = 0.5,
  grid_slice_max       = 2,
  minimum_polygon_area = 0.003,
  cylinder_fit_type    = "ransac",
  max_dia              = 1,
  SDvert               = 0.50,
  n_pts                = 20,
  n_best               = 25,
  inliers              = 0.9,
  conf                 = 0.99,
  max_angle            = 20
)

cat("Trees detected:", nrow(tree_locs), "\n")
cat("Tree attributes:\n")
print(tree_locs)

# -- CHM + tree locations (ggplot2) ------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  chm    <- lidR::rasterize_canopy(las, res = 0.1, lidR::p2r(subcircle = 0.01))
  chm_df <- as.data.frame(chm, xy = TRUE)
  z_col  <- setdiff(names(chm_df), c("x", "y"))[1L]
  names(chm_df)[names(chm_df) == z_col] <- "Height"
  chm_df <- chm_df[!is.na(chm_df$Height), ]

  tl_coords   <- sf::st_coordinates(tree_locs)
  tl_df       <- data.frame(TreeID = tree_locs$TreeID,
                             x = tl_coords[, 1L], y = tl_coords[, 2L],
                             DBH_cm = tree_locs$Radius * 200)
  # Buffer each stem location by its fitted radius → circles at true map scale
  tl_circles  <- sf::st_buffer(tree_locs, dist = tree_locs$Radius)

  print(
    ggplot2::ggplot() +
      ggplot2::geom_raster(data = chm_df, ggplot2::aes(x, y, fill = Height)) +
      ggplot2::scale_fill_viridis_c("Height (m)", option = "D", na.value = NA) +
      ggplot2::geom_sf(data = tl_circles, fill = NA, colour = "black",
                       linewidth = 1.0, inherit.aes = FALSE) +
      ggplot2::geom_text(data = tl_df,
                         ggplot2::aes(x, y, label = paste0("T", TreeID, "\n",
                                                            round(DBH_cm, 1), " cm")),
                         vjust = -1.1, colour = "white", size = 2.8, fontface = "bold") +
      ggplot2::coord_sf() +
      ggplot2::labs(title    = "Tree locations and DBH over canopy height model",
                    subtitle = paste0(nrow(tl_df), " trees detected  |  CHM res = 10 cm  |  ", las_label),
                    x = "Easting (m)", y = "Northing (m)") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                     plot.subtitle = ggplot2::element_text(colour = "grey40"))
  )
  cat("  [plot] Tree locations and DBH rendered over CHM\n")
} else {
  cat("  [plot skipped] Install ggplot2 to see the CHM + tree-location map\n")
}

# =============================================================================
# 3. Individual tree segmentation  (shared by both pipelines)
# =============================================================================
cat("\n--- Step 2: Tree segmentation (segment_graph) ---\n")

las_seg <- segment_graph(
  las,
  tree_locs,
  k                    = 50,
  distance.threshold   = 0.5,
  use.metabolic.scale  = FALSE,
  ptcloud_slice_min    = 0.5,
  ptcloud_slice_max    = 2.5,
  subsample.graph      = 0.1,
  return.dense         = TRUE
)

n_trees_seg <- length(unique(las_seg@data$treeID[las_seg@data$treeID > 0]))
cat("Trees with at least one point:", n_trees_seg, "\n")

# las_seg <- colorize_las(las_seg, method = "attr", attribute_name="treeID",
#                     palette=random.colors(100))
# lidRviewer::view(filter_poi(las_seg, Classification != 2))

# las_seg contains every point from las plus the treeID column; las is no
# longer needed and can be freed before the expensive pipeline steps.
rm(las); gc()

# =============================================================================
# PIPELINE 1 – Original:  process_tree_data()
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("PIPELINE 1 – Original tree processing\n")
cat(strrep("=", 60), "\n")

cat("\n--- Step 3a: Post-process trees (process_tree_data) ---\n")

tree_data <- process_tree_data(
  tree_locs,
  las_seg,
  return_sf = FALSE
)

cat("Tree inventory (Pipeline 1):\n")
print(tree_data)

# Quick summary
cat("\nSummary statistics:\n")
cat("  DBH (cm)  – mean:", round(mean(tree_data$Radius * 200, na.rm = TRUE), 1),
    " SD:", round(sd(tree_data$Radius * 200, na.rm = TRUE), 1), "\n")
cat("  Height (m) – mean:", round(mean(tree_data$height, na.rm = TRUE), 1),
    " SD:", round(sd(tree_data$height, na.rm = TRUE), 1), "\n")

# =============================================================================
# PIPELINE 2 – New leaf–wood classification pipeline
# =============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("PIPELINE 2 – Bole segmentation pipeline\n")
cat(strrep("=", 60), "\n")

# ---- Step 3b: Classify bole + branch points --------------------------------
cat("\n--- Step 3b: Leaf-wood classification (classify_lw_points) ---\n")

las_bole <- classify_lw_points(
  las_seg,
  tree_locs,
  method = "lewos",
  z_min = 0.1,
  z_max = 30,
  search_radius = 0.6,
  neigh_k = 20L,
  # Bole detection uses the COARSE scale only (0.25 m radius gives stable
  # eigenvalues on the cylinder surface).  Fine-scale (0.05 m) is only used
  # for branch linearity — too sparse on a voxelized bole to give reliable
  # Verticality.  Spatial discrimination is handled by search_radius (0.6 m).
  verticality_threshold = 0.75,
  n_propagation         = 5L,
  k_propagation         = 15L,
  ncpu       = 15L,
  voxel_thin = 0.05,   # thin to 5-cm voxels; reused for both eigen passes
  branch_r   = 0.05,   # fine scale: thin branches are very linear at 5 cm
  coarse_r   = 0.25,   # coarse scale: boles stay vertical, branches diverge
  branch_lin_threshold = 0.35,
  branch_cc_min_size   = 10L
)

n_bole   <- sum(las_bole@data$Bole,   na.rm = TRUE)
n_branch <- sum(las_bole@data$Branch, na.rm = TRUE)
n_other  <- sum(las_bole@data$Other,  na.rm = TRUE)
cat(sprintf("  Bole   : %d  (%.1f%%)\n", n_bole,
            100 * n_bole   / nrow(las_bole@data)))
cat(sprintf("  Branch : %d  (%.1f%%)\n", n_branch,
            100 * n_branch / nrow(las_bole@data)))
cat(sprintf("  Other  : %d  (%.1f%%)\n", n_other,
            100 * n_other  / nrow(las_bole@data)))

# colorize inline – avoids copying a full RGB-column LAS into a second object
if (requireNamespace("lidRviewer", quietly = TRUE))
  lidRviewer::view(
    colorize_las(las_bole, method = "attr", attribute_name = "BoleProb",
                 palette = height.colors(100))
  )

# ---- Step 3c: Branch metrics (eigen_metrics + branch_metrics) --------------
cat("\n--- Step 3c: Branch metrics (eigen_metrics → branch_metrics) ---\n")

# branch_metrics() requires eigen_metrics() output that includes E1x/E1y/E1z.
# We run it on non-bole, non-ground points at two neighbourhood scales:
#   Fine  (r = 0.05 m) – captures thin branch linearity
#   Coarse (r = 0.25 m) – same scale used for bole detection; validates that
#                          bole-scale metrics are NOT selected as branches.

las_nonbole <- lidR::filter_poi(las_bole, !Bole & Classification != 2L)
cat(sprintf("  Non-bole/non-ground points : %d\n", nrow(las_nonbole@data)))

# ---- 3c-i  Fine scale (r = 0.05 m) ----------------------------------------
cat("\n  [3c-i] Fine scale (r = 0.05 m)\n")
em_fine <- eigen_metrics(las_nonbole, r = 0.05, ncpu = 15L)

# Verify E1x/E1y/E1z are present before calling branch_metrics
stopifnot(all(c("E1x", "E1y", "E1z") %in% names(em_fine)))
cat("  E1x/E1y/E1z present: TRUE\n")

# Attach height so branch_metrics applies the >1.37 m height filter
em_fine[, Zref := las_nonbole@data$Z]

branch_metrics(em_fine,
               min_angle      = 45,    # ≥ 45° from vertical
               max_angle      = 90,    # up to horizontal
               score_quantile = 0.75)  # top-quartile threshold

cat(sprintf("  Points scored          : %d\n",  nrow(em_fine)))
cat(sprintf("  IsBranchCandidate TRUE : %d  (%.1f%%)\n",
            sum(em_fine$IsBranchCandidate, na.rm = TRUE),
            100 * mean(em_fine$IsBranchCandidate, na.rm = TRUE)))
cat("  AxisAngle quantiles (°)  :", paste(
    round(quantile(em_fine$AxisAngle,   c(0, .25, .5, .75, .9, 1), na.rm = TRUE), 1),
    collapse = "  "), "\n")
cat("  BranchScore quantiles    :", paste(
    round(quantile(em_fine$BranchScore, c(0, .25, .5, .75, .9, 1), na.rm = TRUE), 3),
    collapse = "  "), "\n")

# ---- 3c-ii  Coarse scale (r = 0.25 m) – sensitivity check -----------------
cat("\n  [3c-ii] Coarse scale (r = 0.25 m)\n")
em_coarse <- eigen_metrics(las_nonbole, r = 0.25, ncpu = 15L)
em_coarse[, Zref := las_nonbole@data$Z]

branch_metrics(em_coarse,
               min_angle      = 45,
               max_angle      = 90,
               score_quantile = 0.75)

cat(sprintf("  IsBranchCandidate TRUE : %d  (%.1f%%)\n",
            sum(em_coarse$IsBranchCandidate, na.rm = TRUE),
            100 * mean(em_coarse$IsBranchCandidate, na.rm = TRUE)))
# At the coarse scale bole-like vertical structure dominates, so the branch
# candidate fraction should be lower than at the fine scale.
cat(sprintf("  [check] Fine vs coarse candidate %%: %.1f vs %.1f — fine should be higher\n",
            100 * mean(em_fine$IsBranchCandidate,   na.rm = TRUE),
            100 * mean(em_coarse$IsBranchCandidate, na.rm = TRUE)))

# ---- 3c-iii  score_quantile sensitivity ------------------------------------
cat("\n  [3c-iii] score_quantile sensitivity (fine scale)\n")
for (q in c(0.50, 0.75, 0.90)) {
  em_q <- data.table::copy(em_fine)
  branch_metrics(em_q, min_angle = 45, max_angle = 90, score_quantile = q)
  cat(sprintf("    score_quantile = %.2f  →  %d candidates  (%.1f%%)\n",
              q,
              sum(em_q$IsBranchCandidate, na.rm = TRUE),
              100 * mean(em_q$IsBranchCandidate, na.rm = TRUE)))
}

# ---- 3c-iv  Agreement with classify_lw_points() Branch flag ----------------
cat("\n  [3c-iv] Agreement with classify_lw_points() Branch flag\n")
# em_fine rows correspond 1-to-1 with las_nonbole points.
lw_branch   <- las_nonbole@data$Branch
bm_branch   <- em_fine$IsBranchCandidate

agree       <- sum(lw_branch == bm_branch, na.rm = TRUE)
both_true   <- sum(lw_branch & bm_branch,  na.rm = TRUE)
lw_only     <- sum(lw_branch & !bm_branch, na.rm = TRUE)
bm_only     <- sum(!lw_branch & bm_branch, na.rm = TRUE)

cat(sprintf("  Overall agreement          : %d / %d  (%.1f%%)\n",
            agree, nrow(em_fine), 100 * agree / nrow(em_fine)))
cat(sprintf("  Both flagged branch        : %d\n", both_true))
cat(sprintf("  classify_lw only           : %d\n", lw_only))
cat(sprintf("  branch_metrics only        : %d\n", bm_only))
# Note: perfect agreement is not expected — classify_lw_points uses connected
# components and potentially two-scale eigen; branch_metrics is a standalone
# point-level score without spatial connectivity filtering.

# ---- 3c-v  Per-tree candidate breakdown ------------------------------------
cat("\n  [3c-v] Per-tree IsBranchCandidate breakdown (fine scale, top 10 trees)\n")
em_fine[, treeID := las_nonbole@data$treeID]
per_tree <- em_fine[, .(
  N_pts     = .N,
  N_branch  = sum(IsBranchCandidate, na.rm = TRUE),
  pct       = round(100 * mean(IsBranchCandidate, na.rm = TRUE), 1),
  mean_score = round(mean(BranchScore, na.rm = TRUE), 3),
  mean_angle = round(mean(AxisAngle,   na.rm = TRUE), 1)
), by = treeID][order(-N_branch)][1:min(10L, .N)]
print(per_tree)

# ---- 3c-vi  Visualisation --------------------------------------------------
if (requireNamespace("lidRviewer", quietly = TRUE)) {
  # Map IsBranchCandidate back onto the full las_nonbole cloud for viewing.
  las_nonbole@data[, BranchCandidate := em_fine$IsBranchCandidate]
  lidRviewer::view(
    colorize_las(las_nonbole, method = "attr",
                 attribute_name = "BranchCandidate",
                 palette = c("grey70", "firebrick"))
  )
}

rm(las_nonbole, em_fine, em_coarse, em_q, per_tree, lw_branch, bm_branch)
gc()


cat("\n--- Step 4: Bole segmentation (segment_bole) ---\n")
cat("    Algorithm: ransac  (C++ implementation)\n")

seg_tbl <- segment_bole(
  las_bole,
  tree_locs,
  algorithm       = "ransac",
  z_min           = 0.1,
  z_max           = NULL,     # per-tree maximum Z
  dz              = 1.0,      # 1-m slices (fewer slices = faster)
  overlap         = 0.10,     # 10% overlap between slices
  use_stem_column = TRUE,     # uses Bole column (classify_lw_points); falls back to Stem
  search_radius   = NULL,     # auto
  inlier_tol      = 0.03,
  # n_samples must be > 3 so the system is overdetermined; with exactly 3 points
  # a circle is solved exactly (machine-epsilon error) giving meaningless RANSAC_err
  # and no reliable inlier signal.  10 gives a good overdetermined QR fit.
  # kIterations = (n_best+5) × ceil(log(1-conf)/log(1-inlier_frac^n_ransac))
  # With n_ransac=10, n_best=5, conf=0.95: ~140 iters/slice (vs 440 with n_best=15, conf=0.99)
  n_ransac        = 10L,
  conf            = 0.95,
  inlier_frac     = 0.85,
  n_best          = 5L,
  max_radius      = 1.0,
  min_pts         = 10L,
  # voxel_thin skipped: cloud was already thinned to 5-cm voxels upstream in
  # classify_lw_points(), so re-thinning here only costs time.
  voxel_thin      = NULL,
  ncpu            = 15L    # trees are independent → parallelise; nested OMP disabled internally
)

cat("  Total slices fitted:", nrow(seg_tbl), "\n")
valid_segs <- seg_tbl[!is.na(seg_tbl$Valid) & seg_tbl$Valid, ]
cat("  Valid slices:", nrow(valid_segs), "\n")
cat("  Slices per tree (valid):\n")
print(table(valid_segs$TreeID))

cat("\n  Segment table (first 10 rows):\n")
print(head(seg_tbl, 10), digits = 4)

# seg_tbl has the fitted geometry; the classified point cloud is no longer needed.
rm(las_bole); gc()

# Optionally try with the Hough algorithm (requires houghPrimitives):
# seg_tbl <- segment_bole(las_bole, tree_locs, algorithm = "hough")
# seg_tbl <- segment_bole(las_bole, tree_locs, algorithm = "hough_seed_ransac")

# ---- Step 5: Volume estimation ---------------------------------------------
cat("\n--- Step 5: Volume estimation (compute_bole_volume) ---\n")

vol_cyl     <- compute_bole_volume(seg_tbl, method = "cylinder")
vol_smalian <- compute_bole_volume(seg_tbl, method = "smalian")

cat("  Cylinder formula:\n")
print(vol_cyl)
cat("\n  Smalian formula:\n")
print(vol_smalian)

# Difference between methods
diff_pct <- abs(vol_cyl$Volume_m3 - vol_smalian$Volume_m3) /
  pmax(vol_cyl$Volume_m3, 1e-9) * 100
cat("\n  Volume difference between methods (%):",
    paste(round(diff_pct, 1), collapse = ", "), "\n")

# ---- Step 6: Quality assessment --------------------------------------------
cat("\n--- Step 6: Quality assessment (assess_tree_quality) ---\n")

qual_tbl <- assess_tree_quality(
  tree_locs,
  seg_tbl,
  vol_table          = vol_cyl,
  max_mean_error     = 0.03,   # realistic for 10-pt overdetermined RANSAC on MLS
  max_max_error      = 0.08,
  min_inlier_frac    = 0.40,   # MLS slices typically see 40–60% of points on-circle
  min_bole_coverage  = 0.50,
  dbh_ht_ratio_range = c(0.004, 0.05),
  max_radius_cv      = 0.40,
  min_segments       = 3L
)

cat("  Quality summary:\n")
print(table(qual_tbl$quality_class))
cat("\n  Per-tree quality metrics:\n")
metric_cols <- c("TreeID", "N_segments", "DBH_m", "fit_score",
                 "Mean_RANSAC_err", "Mean_inlier_frac",
                 "Bole_coverage", "quality_class")
metric_cols <- metric_cols[metric_cols %in% names(qual_tbl)]
print(sf::st_drop_geometry(qual_tbl)[, metric_cols], digits = 3)

# ---- Step 7: Iterative refit (only when bad/marginal trees exist) ----------
n_bad <- sum(qual_tbl$quality_class %in% c("bad", "marginal"), na.rm = TRUE)

cat("\n--- Step 7: Iterative refit (refit_trees) ---\n")

if (n_bad == 0L) {
  cat("  All trees are 'good' — no refit needed.\n")
  final_seg  <- seg_tbl
  final_qual <- qual_tbl
} else {
  cat("  Attempting refit for", n_bad,
      "bad/marginal tree(s)...\n")

  refit_out <- refit_trees(
    las_seg,                          # original (not bole-only) LAS
    tree_locs,
    seg_tbl,
    qual_tbl,
    strategies     = c("relax_ransac", "tighten_ransac",
                       "meanshift_bole"),  # dbscan_bole needs dbscan package
    refit_marginal = TRUE,
    max_iterations = 2L,
    z_min          = 0.1,
    dz             = 0.5,
    overlap        = 0.25,
    inlier_tol     = 0.03,
    n_ransac       = 20L,
    inlier_frac    = 0.85,
    n_best         = 25L
  )

  cat("\n  Refit log:\n")
  print(refit_out$refit_log)

  cat("\n  Updated quality summary:\n")
  print(table(refit_out$qual_table$quality_class))

  final_seg  <- refit_out$seg_table
  final_qual <- refit_out$qual_table
}

# =============================================================================
# PIPELINE 2 UNIFIED OUTPUT – process_tree_data with all pipeline-2 tables
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("PIPELINE 2 UNIFIED OUTPUT – process_tree_data with all tables\n")
cat(strrep("=", 60), "\n")

cat("\n--- Step 8: Unified tree table (process_tree_data + pipeline-2 args) ---\n")

tree_data_full <- process_tree_data(
  tree_locs,
  las_seg,
  return_sf  = FALSE,
  seg_table  = final_seg,     # from segment_bole / refit_trees
  vol_table  = vol_cyl,       # from compute_bole_volume
  qual_table = final_qual     # from assess_tree_quality / refit_trees
)

cat("\n  Columns in unified tree table:\n")
cat("  ", paste(names(tree_data_full), collapse = ", "), "\n")

cat("\n  Unified tree table (all rows):\n")
print(sf::st_drop_geometry(tree_data_full), digits = 3)

# Compare Pipeline 1 height vs Pipeline 2 height
if ("height" %in% names(tree_data) && "height" %in% names(tree_data_full)) {
  ht_compare <- merge(
    data.frame(TreeID = tree_data$TreeID, HtP1 = tree_data$height),
    data.frame(TreeID = tree_data_full$TreeID, HtP2 = tree_data_full$height),
    by = "TreeID", all = TRUE
  )
  ht_compare$Delta_m <- ht_compare$HtP2 - ht_compare$HtP1
  cat("\n  Height comparison P1 vs P2 (should match – same segmented LAS):\n")
  print(ht_compare, digits = 3)
}

# Show volume and quality columns that are only in the full table
if ("Volume_m3" %in% names(tree_data_full)) {
  cat("\n  Bole volume and quality (pipeline-2 only columns):\n")
  p2_cols <- intersect(
    c("TreeID", "Volume_m3", "N_segments", "Mean_radius", "Basal_area_m2",
      "Height_modeled", "Mean_RANSAC_err", "DBH_m", "fit_score", "quality_class"),
    names(tree_data_full)
  )
  print(sf::st_drop_geometry(tree_data_full)[, p2_cols, drop = FALSE], digits = 3)
}

# =============================================================================
# 4. Compare Pipeline 1 vs Pipeline 2 DBH estimates
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("DBH comparison – Pipeline 1 vs Pipeline 2\n")
cat(strrep("=", 60), "\n")

# Pipeline 1 uses tree_locs$Radius (from RANSAC cylinder at breast height)
# Pipeline 2 fits a circle per 0.5-m slice; we grab the lowest valid slice
# as the "breast height" estimate
dbh_p1 <- data.frame(
  TreeID  = tree_locs$TreeID,
  DBH_cm_P1 = tree_locs$Radius * 200    # radius → DBH in cm
)

dbh_p2 <- do.call(rbind, lapply(split(final_seg, final_seg$TreeID), function(tr) {
  v <- tr[!is.na(tr$Valid) & tr$Valid, ]
  if (nrow(v) == 0L) return(data.frame(TreeID = tr$TreeID[1], DBH_cm_P2 = NA_real_))
  lowest <- v[which.min(v$Z_mid), ]
  data.frame(TreeID = lowest$TreeID, DBH_cm_P2 = lowest$Radius * 200)
}))

dbh_compare <- merge(dbh_p1, dbh_p2, by = "TreeID")
dbh_compare$Delta_cm <- dbh_compare$DBH_cm_P2 - dbh_compare$DBH_cm_P1
cat("\n")
print(dbh_compare, digits = 3)
cat("\n  Mean |ΔDBH| between pipelines:",
    round(mean(abs(dbh_compare$Delta_cm), na.rm = TRUE), 2), "cm\n")

# =============================================================================
# 5. Final summary
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("FINAL SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat("Pipeline 1 (original):  process_tree_data, crown geometry only\n")
cat("  Trees:          ", nrow(tree_data), "\n")
cat("  Mean DBH (cm):  ",
    round(mean(tree_data$Radius * 200, na.rm = TRUE), 1), "\n")
cat("  Mean Height (m):",
    round(mean(tree_data$height, na.rm = TRUE), 1), "\n")
cat("  Columns:        ", paste(setdiff(names(tree_data), "geometry"), collapse = ", "), "\n")

cat("\nPipeline 2 (bole segmentation) – separate tables:\n")
cat("  Trees detected: ", nrow(tree_locs), "\n")
cat("  Valid segments: ", nrow(valid_segs), "\n")
n_good <- sum(final_qual$quality_class == "good",     na.rm = TRUE)
n_marg <- sum(final_qual$quality_class == "marginal", na.rm = TRUE)
n_bad2 <- sum(final_qual$quality_class == "bad",      na.rm = TRUE)
cat("  Quality: good =", n_good, "| marginal =", n_marg, "| bad =", n_bad2, "\n")
cat("  Total volume (cylinder, m\u00b3):",
    paste(round(vol_cyl$Volume_m3, 4), collapse = " | "), "\n")

cat("\nPipeline 2 unified  – process_tree_data with seg/vol/qual tables:\n")
cat("  Trees:          ", nrow(tree_data_full), "\n")
cat("  Columns:        ", paste(setdiff(names(tree_data_full), "geometry"), collapse = ", "), "\n")
if ("Volume_m3" %in% names(tree_data_full)) {
  cat("  Mean volume (m\u00b3):",
      round(mean(tree_data_full$Volume_m3, na.rm = TRUE), 4), "\n")
}
if ("quality_class" %in% names(tree_data_full)) {
  cat("  Quality summary:\n")
  print(table(tree_data_full$quality_class))
}

cat("\nDone.\n")
