library(spanner)
library(lidR)

# Helper: load, classify ground, normalize, drop ground points
prep_las <- function(path, crs = NULL) {
  las <- readLAS(path)
  if (!is.null(crs)) sf::st_crs(las) <- crs
  las <- classify_ground(las, ptd())
  las <- normalize_height(las, tin(classified = FALSE))
  las <- filter_poi(las, Z > 0)
  las
}

run_pipeline <- function(las, label, map_method, assign_method,
                         stem_method, seg_method) {
  cat("\n", strrep("=", 66), "\n", sep = "")
  cat(" ", label, "\n")
  cat(strrep("=", 66), "\n", sep = "")
  cat("Points after normalization:", npoints(las), "\n")

  # Stage 1 -- find_trees
  cat("\n[Stage 1] find trees\n")
  map <- find_trees(las, method = map_method)
  cat("  Trees detected:", nrow(map), "\n")
  if (nrow(map) == 0L) { cat("  No trees found -- skipping.\n"); return(invisible(NULL)) }

  # Stage 2 -- tree_points
  cat("\n[Stage 2] classify tree points\n")
  las_t <- tree_points(las, map, method = assign_method)
  n_assigned <- sum(!is.na(las_t$treeID) & las_t$treeID > 0L)
  cat("  Points assigned to a tree:", n_assigned, "\n")

  # Stage 3 -- stem_points
  cat("\n[Stage 3] classify stem points\n")
  las_s <- stem_points(las_t, map, method = stem_method)
  cat("  Stem points:", sum(las_s$Stem, na.rm = TRUE), "\n")

  # Stage 4 -- stem_segments
  cat("\n[Stage 4] delinate stem segments\n")
  segs <- stem_segments(las_s, map, method = seg_method)
  cat("  Total segments:", nrow(segs), "  Valid:", sum(segs$Valid), "\n")
  if (sum(segs$Valid) > 0L)
    cat("  Radius range (valid):",
        round(range(segs$Radius[segs$Valid], na.rm = TRUE), 3), "m\n")

  # Stage 5 -- tree_inventory (auto-computes volume from seg_table)
  cat("\n[Stage 5] run tree inventory\n")
  inv <- tree_inventory(map, las_t, seg_table = segs)
  cols_show <- intersect(c("TreeID", "height", "crown_area",
                           "Volume_m3", "Mean_radius"), names(inv))
  print(sf::st_drop_geometry(inv)[, cols_show, drop = FALSE])

  # compute_stem_volume standalone
  cat("\n  compute stem volume:\n")
  vol <- compute_stem_volume(segs)
  print(vol[, c("TreeID", "Volume_m3", "N_segments", "Mean_radius")])

  # assess_tree_quality
  cat("\n  assess tree quality:\n")
  qual <- tryCatch(
    assess_tree_quality(map, segs,
                        max_mean_error = 0.05,
                        max_max_error = 0.2,
                        min_inlier_frac = 0.7,
                        min_stem_coverage = 0.5,
                        dbh_ht_ratio_range = c(0.004, 0.05),
                        max_radius_cv = 0.4,
                        min_segments = 3L),
    error = function(e) { cat("  (skipped:", conditionMessage(e), ")\n"); NULL }
  )
  if (!is.null(qual)) print(table(qual$quality_class))

  invisible(list(map = map, las_t = las_t, las_s = las_s, segs = segs, inv = inv, qual = qual))
}

# ==============================================================================
# MLS -- use_hough + assign_graph + stem_eigen + seg_ransac_cylinder
# ==============================================================================
mls_path <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
las_mls  <- prep_las(mls_path)

res_mls <- run_pipeline(
  las           = las_mls,
  label         = "MLS | use_hough + assign_graph + stem_eigen + seg_ransac_cylinder",
  map_method    = use_hough(min_h = 1, max_h = 3, h_step = 0.5,
                            pixel_size = 0.025, max_d = 0.5),
  assign_method = assign_graph(),
  stem_method   = stem_eigen(neigh_k = 30L, verticality_threshold = 0.65),
  seg_method    = seg_ransac_cylinder()
)

# MLS -- stem_hough as alternative stage-3 method
cat("\n--- MLS: stem_hough variant ---\n")
if (!is.null(res_mls)) {
  las_sh <- stem_points(res_mls$las_t, res_mls$map,
                        method = stem_hough(h_base_lo = 1, h_base_hi = 2.5,
                                           h_step = 0.5, max_d = 0.5,
                                           pixel_size = 0.05, min_density = 0.1,
                                           min_votes = 3L))
  cat("  stem_hough Stem points:", sum(las_sh$Stem, na.rm = TRUE), "\n")
}

# MLS -- seg_irls_cylinder comparison
cat("\n--- MLS: seg_irls_cylinder comparison ---\n")
if (!is.null(res_mls)) {
  map_mls  <- res_mls$map
  las_mls_t <- tree_points(las_mls, map_mls, assign_graph())
  las_mls_s <- stem_points(las_mls_t, map_mls,
                           stem_eigen(neigh_k = 30L, verticality_threshold = 0.65))
  segs_ir <- stem_segments(las_mls_s, map_mls, seg_irls_cylinder())
  cat("  seg_irls_cylinder valid:", sum(segs_ir$Valid), "\n")
}

# ==============================================================================
# TLS -- use_raster_eigen + assign_voronoi + stem_eigen + seg_irls_cylinder
# ==============================================================================
tls_path <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
las_tls  <- prep_las(tls_path, crs = 26912)

res_tls <- run_pipeline(
  las           = las_tls,
  label         = "TLS | use_raster_eigen + assign_voronoi + stem_eigen + seg_irls_cylinder",
  map_method    = use_raster_eigen(
    res = 0.025, pt_spacing = 0.0254, dens_threshold = 0.25,
    neigh_sizes = c(0.25, 0.15, 0.66), eigen_threshold = 0.75,
    grid_slice_min = 0.5, grid_slice_max = 2.5,
    minimum_polygon_area = 0.005, cylinder_fit_type = "ransac",
    max_dia = 1, SDvert = 0.33, n_pts = 20, n_best = 25,
    inliers = 0.9, conf = 0.99, max_angle = 20
  ),
  assign_method = assign_voronoi(),
  stem_method   = stem_eigen(),
  seg_method    = seg_irls_cylinder()
)

# TLS -- use_hough as fast stage-1 alternative
cat("\n--- TLS: use_hough() as stage-1 alternative ---\n")
map_tls_h <- tryCatch(
  find_trees(las_tls, use_hough(min_h = 1, max_h = 3, h_step = 0.5,
                                pixel_size = 0.025, max_d = 0.6)),
  error = function(e) { cat("  (failed:", conditionMessage(e), ")\n"); NULL }
)
if (!is.null(map_tls_h))
  cat("  use_hough on TLS:", nrow(map_tls_h), "trees\n")


# ==============================================================================
# Supplementary tests -- TLS/MLS-relevant functions not covered above
# ==============================================================================

# -- assign_crop() (MLS, alternative stage-2 method) --------------------------
cat("\n--- assign_crop() ---\n")
if (!is.null(res_mls)) {
  las_crop <- tree_points(las_mls, res_mls$map,
                          method = assign_crop(l = 3))
  n_crop <- sum(!is.na(las_crop$treeID) & las_crop$treeID > 0L)
  cat("  assign_crop points assigned:", n_crop, "\n")
}

# -- process_tree_data() (v1 API; auto-computes volume from seg_table) ---------
cat("\n--- process_tree_data() ---\n")
if (!is.null(res_tls)) {
  ptd_res <- tryCatch(
    process_tree_data(res_tls$map, res_tls$las_t,
                      seg_table  = res_tls$segs,
                      qual_table = res_tls$qual,
                      return_sf  = FALSE),
    error = function(e) { cat("  (failed:", conditionMessage(e), ")\n"); NULL }
  )
  if (!is.null(ptd_res))
    cat("  process_tree_data rows:", nrow(ptd_res),
        "  cols:", paste(names(ptd_res), collapse = ", "), "\n")
}

# -- eigen_metrics() + branch_metrics() (TLS) ----------------------------------
cat("\n--- eigen_metrics() + branch_metrics() ---\n")
em <- tryCatch(
  eigen_metrics(las_tls, r = 0.1, k = 10L, ncpu = 2L),
  error = function(e) { cat("  (failed:", conditionMessage(e), ")\n"); NULL }
)
if (!is.null(em)) {
  cat("  eigen_metrics rows:", nrow(em),
      "  E1x/E1y/E1z present:", all(c("E1x", "E1y", "E1z") %in% names(em)), "\n")
  branch_metrics(em, min_angle = 45, max_angle = 90, score_quantile = 0.75)
  cat("  Branch candidates:", sum(em$IsBranchCandidate, na.rm = TRUE),
      "/", nrow(em), "\n")
}

# -- classify_lw_points() standalone (TLS) ------------------------------------
cat("\n--- classify_lw_points() standalone ---\n")
if (!is.null(res_tls)) {
  las_lw <- tryCatch(
    classify_lw_points(res_tls$las_t, res_tls$map,
                       method                = "eigen",
                       neigh_k               = 20L,
                       verticality_threshold = 0.75,
                       ncpu                  = 2L),
    error = function(e) { cat("  (failed:", conditionMessage(e), ")\n"); NULL }
  )
  if (!is.null(las_lw))
    cat("  Stem points:", sum(las_lw$Stem, na.rm = TRUE),
        "  Branch points:", sum(las_lw$Branch, na.rm = TRUE), "\n")
}

# -- refit_trees() (TLS) -------------------------------------------------------
cat("\n--- refit_trees() ---\n")
if (!is.null(res_tls) && !is.null(res_tls$qual)) {
  n_bad <- sum(res_tls$qual$quality_class %in% c("bad", "marginal"), na.rm = TRUE)
  cat("  bad/marginal before refit:", n_bad, "\n")
  if (n_bad > 0L) {
    refit_res <- tryCatch(
      refit_trees(res_tls$las_t, res_tls$map, res_tls$segs, res_tls$qual,
                  strategies     = c("relax_ransac", "tighten_ransac",
                                     "meanshift_stem"),
                  refit_marginal = TRUE),
      error = function(e) { cat("  (failed:", conditionMessage(e), ")\n"); NULL }
    )
    if (!is.null(refit_res)) {
      cat("  Quality after refit:\n")
      print(table(refit_res$qual_table$quality_class))
      print(refit_res$refit_log)
    }
  } else {
    cat("  No bad/marginal trees -- skipping refit.\n")
  }
}

cat("\nDone.\n")
