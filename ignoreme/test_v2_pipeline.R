# spanner v2 pipeline walkthrough
#
# This script is meant to be both a smoke test and a copy/paste starting point
# for new data. It runs the same workflow twice:
#   1. MLS example data, then refit weak trees.
#   2. TLS example data, then refit weak trees.
#
# The important tuning pattern is:
#   detect trees -> assign points -> find stem points -> fit stem segments
#   -> score quality -> refit only bad/marginal trees
#
# For your own data, start by changing the paths below, then tune the method
# arguments in the "MLS settings" or "TLS settings" sections.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", helpers = FALSE, attach_testthat = FALSE)
} else {
  library(spanner)
}
library(lidR)

# Set SPANNER_RUN_OPTIONAL=true to try slower alternative methods at the end.
run_optional <- tolower(Sys.getenv("SPANNER_RUN_OPTIONAL", "false")) %in%
  c("1", "true", "yes")

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

# Replace these with your own .las/.laz files when using this as a template.
mls_path <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
tls_path <- system.file("extdata", "TLS_Clip.laz", package = "spanner")

# Use a CRS if your file does not already carry one. The TLS example needs it.
tls_crs <- 26912

# ---------------------------------------------------------------------------
# Quality and refit settings
# ---------------------------------------------------------------------------

# These values decide which trees are "good", "marginal", or "bad".
# For v2 segmenters, RANSAC error/inlier diagnostics may be unavailable; the
# quality scorer will use the metrics that exist instead of penalizing NA.
quality_args <- list(
  max_mean_error       = 0.05,
  max_max_error        = 0.20,
  min_inlier_frac      = 0.70,
  min_stem_coverage    = 0.50,
  expected_stem_height = 3.0,
  dbh_ht_ratio_range   = c(0.004, 0.05),
  max_radius_cv        = 0.40,
  min_segments         = 3L
)

# Keep refit narrow and quick. Raise max_trees or max_iterations after the
# detection/segmentation settings look sane.
refit_settings <- list(
  max_trees             = 10L,
  strategies            = c("tighten_ransac"),
  max_iterations        = 1L,
  z_max                 = 2.5,
  n_ransac              = 8L,
  n_best                = 5L,
  max_ransac_iterations = 100L,
  ncpu                  = 4L,
  refit_marginal = FALSE
)

# ---------------------------------------------------------------------------
# MLS settings
# ---------------------------------------------------------------------------

# MLS is often sparser and more directional than TLS. Hough detection works
# well when stems are visible as circular slices near breast height.
#
# To improve detections:
#   - If trees are missed, widen min_h/max_h or increase max_d.
#   - If false trees appear, raise min_h, lower max_d, or increase pixel_size.
#   - If stems are broken/patchy, loosen verticality_threshold a little.
mls_methods <- list(
  map = use_hough(
    min_h = 0.5,
    max_h = 4,
    h_step = 0.5,
    pixel_size = 0.025,
    max_d = 0.7,
    min_density = 0.05,
    min_votes = 2L
  ),
  assign = assign_voronoi(max_dist = 2.0),
  stem = stem_eigen(
    neigh_k = 30L,
    linearity_threshold = 0.2,
    verticality_threshold = 0.5,
    axis_refine_stems = TRUE,
    axis_min_points = 5L,
    axis_search_radius = 0.6
  ),
  segment = seg_ransac_cylinder()
)

# ---------------------------------------------------------------------------
# TLS settings
# ---------------------------------------------------------------------------

# TLS normally has denser, cleaner stem structure. Raster/eigen detection is a
# good first choice when many stems are present and the scan is dense.
#
# To improve detections:
#   - If small stems are missed, lower dens_threshold or minimum_polygon_area.
#   - If noisy patches become trees, raise dens_threshold/eigen_threshold.
#   - If diameters are too large, lower max_dia.
tls_methods <- list(
  map = use_raster_eigen(
    res = 0.025,
    pt_spacing = 0.0254,
    dens_threshold = 0.20,
    neigh_sizes = c(0.25, 0.15, 0.66),
    eigen_threshold = 0.70,
    grid_slice_min = 0.5,
    grid_slice_max = 2.5,
    minimum_polygon_area = 0.005,
    cylinder_fit_type = "ransac",
    max_dia = 1,
    SDvert = 0.33,
    n_pts = 20,
    n_best = 25,
    inliers = 0.9,
    conf = 0.99,
    max_angle = 20
  ),
  assign = assign_voronoi(max_dist = 2.0),
  stem = stem_eigen(
    neigh_k = 20L,
    linearity_threshold = 0.02,
    verticality_threshold = 0.60,
    axis_refine_stems = TRUE,
    axis_min_points = 10L,
    axis_search_radius = 0.45
  ),
  segment = seg_irls_cylinder()
)

# ---------------------------------------------------------------------------
# Small reporting helpers
# ---------------------------------------------------------------------------

timed <- function(label, expr) {
  cat("\n", label, "\n", sep = "")
  flush.console()
  t <- system.time(value <- force(expr))
  cat(sprintf("  elapsed: %.2f sec\n", unname(t[["elapsed"]])))
  flush.console()
  attr(value, "elapsed_sec") <- unname(t[["elapsed"]])
  value
}

safe <- function(label, expr) {
  tryCatch(
    timed(label, expr),
    error = function(e) {
      cat("  failed: ", conditionMessage(e), "\n", sep = "")
      flush.console()
      NULL
    }
  )
}

print_rule <- function(label) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat(label, "\n", sep = "")
  cat(strrep("=", 72), "\n", sep = "")
}

print_quality <- function(qual) {
  if (is.null(qual)) return(invisible(NULL))
  print(table(qual$quality_class, useNA = "ifany"))
}

print_quality_hints <- function(qual) {
  if (is.null(qual) || nrow(qual) == 0L) return(invisible(NULL))

  bad <- qual[qual$quality_class %in% c("bad", "marginal"), , drop = FALSE]
  if (nrow(bad) == 0L) {
    cat("  all scored trees are good\n")
    return(invisible(NULL))
  }

  cat("  trees to inspect:", paste(head(bad$TreeID, 12L), collapse = ", "))
  if (nrow(bad) > 12L) cat(" ...")
  cat("\n")

  coverage_col <- "Stem_coverage"
  radius_col <- "Radius_CV"
  if (coverage_col %in% names(bad)) {
    low_cov <- bad$TreeID[which(bad[[coverage_col]] < quality_args$min_stem_coverage)]
    if (length(low_cov) > 0L) {
      cat("  low stem coverage:", paste(head(low_cov, 8L), collapse = ", "), "\n")
      cat("    try wider height slices, looser stem extraction, or more refit.\n")
    }
  }
  if (radius_col %in% names(bad)) {
    high_cv <- bad$TreeID[which(bad[[radius_col]] > quality_args$max_radius_cv)]
    if (length(high_cv) > 0L) {
      cat("  unstable radii:", paste(head(high_cv, 8L), collapse = ", "), "\n")
      cat("    try tighter stem points, smaller max_dia, or IRLS segmentation.\n")
    }
  }
}

prep_las <- function(path, crs = NULL) {
  las <- readLAS(path)
  if (!is.null(crs)) sf::st_crs(las) <- crs
  las <- classify_ground(las, ptd())
  las <- normalize_height(las, tin(classified = FALSE))
  filter_poi(las, Z > 0)
}

show_stage_note <- function(stage) {
  notes <- c(
    detect = "Tree detection controls how many candidate stems enter the rest of the pipeline.",
    assign = "Point assignment decides which points belong to each detected tree.",
    stem = "Stem extraction should keep vertical stem points and drop branches/crown clutter.",
    segment = "Segmentation fits stem slices; valid segment count and radius stability matter most.",
    quality = "Quality scoring tells you whether to trust a tree or send it to refit."
  )
  cat("  note:", notes[[stage]], "\n")
}

# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

run_pipeline <- function(las, label, methods) {
  print_rule(label)
  cat("Points after normalization:", npoints(las), "\n")

  show_stage_note("detect")
  map <- timed("[Stage 1] find_trees()", {
    find_trees(las, method = methods$map)
  })
  cat("  trees detected:", nrow(map), "\n")
  if (nrow(map) == 0L) {
    cat("  no trees found; tune the detection settings and rerun\n")
    return(NULL)
  }

  show_stage_note("assign")
  las <- timed("[Stage 2] tree_points()", {
    tree_points(las, map, method = methods$assign)
  })
  n_assigned <- sum(!is.na(las$treeID) & las$treeID > 0L)
  cat("  points assigned to a tree:", n_assigned, "\n")

  show_stage_note("stem")
  las <- timed("[Stage 3] stem_points()", {
    stem_points(las, map, method = methods$stem)
  })
  cat("  stem points:", sum(las$Stem, na.rm = TRUE), "\n")

  show_stage_note("segment")
  segs <- timed("[Stage 4] stem_segments()", {
    stem_segments(las, map, method = methods$segment)
  })
  cat("  total segments:", nrow(segs), "\n")
  cat("  valid segments:", sum(segs$Valid, na.rm = TRUE), "\n")
  if (any(segs$Valid, na.rm = TRUE)) {
    rr <- range(segs$Radius[segs$Valid], na.rm = TRUE)
    cat("  radius range:", paste(round(rr, 3), collapse = " - "), "m\n")
  }

  vol <- timed("[Stage 5] compute_stem_volume()", {
    compute_stem_volume(segs)
  })
  print(vol[, c("TreeID", "Volume_m3", "N_segments", "Mean_radius")])

  show_stage_note("quality")
  qual <- safe("[Stage 6] assess_tree_quality()", {
    do.call(assess_tree_quality, c(list(tree_locs = map, seg_table = segs),
                                   quality_args))
  })
  print_quality(qual)
  print_quality_hints(qual)

  inv <- safe("[Stage 7] tree_inventory()", {
    tree_inventory(map, las, seg_table = segs, vol_table = vol, qual_table = qual)
  })
  if (!is.null(inv)) {
    cols_show <- intersect(c("TreeID", "height", "crown_area",
                             "Volume_m3", "Mean_radius"), names(inv))
    print(sf::st_drop_geometry(inv)[, cols_show, drop = FALSE])
  }

  invisible(list(
    map = map,
    las = las,
    segs = segs,
    vol = vol,
    qual = qual,
    inv = inv
  ))
}

run_refit <- function(res, label, settings = refit_settings) {
  print_rule(paste0(label, " | refit_trees()"))

  if (is.null(res) || is.null(res$qual)) {
    cat("No pipeline result or quality table; refit skipped.\n")
    return(NULL)
  }

  refit_targets <- res$qual$quality_class %in% c("bad", "marginal")
  refit_ids <- head(res$qual$TreeID[refit_targets], settings$max_trees)
  n_refit <- length(refit_ids)

  cat("bad/marginal selected for refit:", n_refit, "\n")
  cat("refit tree IDs:", paste(refit_ids, collapse = ", "), "\n")
  cat("refit budget:", settings$max_iterations, "round(s),",
      settings$max_ransac_iterations, "RANSAC iterations per slice cap\n")

  if (n_refit == 0L) {
    cat("No bad/marginal trees; refit skipped.\n")
    return(NULL)
  }

  out <- safe(paste(label, "timed refit_trees()"), {
    refit_trees(
      res$las,
      res$map,
      res$segs,
      res$qual,
      strategies = settings$strategies,
      refit_marginal = TRUE,
      tree_ids = refit_ids,
      max_iterations = settings$max_iterations,
      z_max = settings$z_max,
      n_ransac = settings$n_ransac,
      n_best = settings$n_best,
      max_ransac_iterations = settings$max_ransac_iterations,
      quality_args = quality_args
    )
  })

  if (!is.null(out)) {
    cat("\nQuality after refit:\n")
    print_quality(out$qual_table)
    print_quality_hints(out$qual_table)
    cat("\nRefit log:\n")
    print(out$refit_log)
    cat(sprintf("\nRefit elapsed seconds: %.2f\n", attr(out, "elapsed_sec")))
  }

  out
}

# ---------------------------------------------------------------------------
# Run MLS, then refit MLS
# ---------------------------------------------------------------------------

las_mls <- timed("Load/prep MLS", prep_las(mls_path))
res_mls <- run_pipeline(
  las = las_mls,
  label = "MLS | Hough detection + graph assignment + eigen stems + RANSAC cylinders",
  methods = mls_methods
)
refit_mls <- run_refit(res_mls, "MLS")

# ---------------------------------------------------------------------------
# Run TLS, then refit TLS
# ---------------------------------------------------------------------------

las_tls <- timed("Load/prep TLS", prep_las(tls_path, crs = tls_crs))
res_tls <- run_pipeline(
  las = las_tls,
  label = "TLS | raster/eigen detection + Voronoi assignment + eigen stems + IRLS cylinders",
  methods = tls_methods
)
refit_tls <- run_refit(res_tls, "TLS")

# ---------------------------------------------------------------------------
# Optional alternatives to try when the main run is not good enough
# ---------------------------------------------------------------------------

if (run_optional && !is.null(res_mls)) {
  safe("MLS optional: stem_hough()", {
    las_sh <- stem_points(
      res_mls$las,
      res_mls$map,
      method = stem_hough(
        h_base_lo = 1, h_base_hi = 2.5, h_step = 0.5,
        max_d = 0.5, pixel_size = 0.05,
        min_density = 0.1, min_votes = 3L
      )
    )
    cat("  stem_hough stem points:", sum(las_sh$Stem, na.rm = TRUE), "\n")
    invisible(las_sh)
  })

  safe("MLS optional: seg_irls_cylinder()", {
    las_mls_s <- stem_points(
      res_mls$las,
      res_mls$map,
      stem_eigen(neigh_k = 30L, verticality_threshold = 0.65)
    )
    segs_ir <- stem_segments(las_mls_s, res_mls$map, seg_irls_cylinder())
    cat("  valid segments:", sum(segs_ir$Valid, na.rm = TRUE), "\n")
    invisible(segs_ir)
  })

  safe("MLS optional: assign_crop()", {
    las_crop <- tree_points(las_mls, res_mls$map, method = assign_crop(l = 3))
    n_crop <- sum(!is.na(las_crop$treeID) & las_crop$treeID > 0L)
    cat("  points assigned:", n_crop, "\n")
    invisible(las_crop)
  })
}

if (run_optional && !is.null(res_tls)) {
  safe("TLS optional: use_hough()", {
    map_tls_h <- find_trees(
      las_tls,
      use_hough(min_h = 1, max_h = 3, h_step = 0.5,
                pixel_size = 0.025, max_d = 0.6)
    )
    cat("  trees detected:", nrow(map_tls_h), "\n")
    invisible(map_tls_h)
  })

  safe("TLS optional: process_tree_data()", {
    ptd_res <- process_tree_data(
      res_tls$map,
      res_tls$las,
      seg_table = res_tls$segs,
      qual_table = res_tls$qual,
      return_sf = FALSE
    )
    cat("  rows:", nrow(ptd_res), "\n")
    invisible(ptd_res)
  })

  safe("TLS optional: eigen_metrics() + branch_metrics()", {
    em <- eigen_metrics(las_tls, r = 0.1, k = 10L, ncpu = 2L)
    branch_metrics(em, min_angle = 45, max_angle = 90, score_quantile = 0.75)
    cat("  rows:", nrow(em), "\n")
    cat("  branch candidates:", sum(em$IsBranchCandidate, na.rm = TRUE), "\n")
    invisible(em)
  })

  safe("TLS optional: classify_lw_points()", {
    las_lw <- classify_lw_points(
      res_tls$las,
      res_tls$map,
      method = "eigen",
      neigh_k = 20L,
      verticality_threshold = 0.75,
      ncpu = 2L
    )
    cat("  stem points:", sum(las_lw$Stem, na.rm = TRUE), "\n")
    cat("  branch points:", sum(las_lw$Branch, na.rm = TRUE), "\n")
    invisible(las_lw)
  })
}

cat("\nDone.\n")
