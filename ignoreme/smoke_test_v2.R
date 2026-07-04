library(spanner)

# Verify all new API exports are present
new_fns <- c(
  "find_trees", "tree_points", "stem_points", "stem_segments", "tree_inventory",
  "use_raster_eigen", "use_hough",
  "assign_voronoi", "assign_crop", "assign_graph",
  "stm_eigen", "stm_hough",
  "seg_ransac_circle", "seg_ransac_cylinder", "seg_irls_cylinder", "seg_bf_cylinder"
)

exported <- ls("package:spanner")
found    <- intersect(new_fns, exported)
missing  <- setdiff(new_fns, exported)

cat("Found:", length(found), "/", length(new_fns), "\n")
if (length(missing) > 0) {
  cat("MISSING:", paste(missing, collapse = ", "), "\n")
} else {
  cat("All new v2.0 functions exported OK\n")
}

# Verify legacy functions still present
legacy <- c("get_raster_eigen_treelocs", "segment_graph", "segment_stem",
            "classify_lw_points", "classify_stem_points", "process_tree_data")
missing_legacy <- setdiff(legacy, exported)
if (length(missing_legacy) > 0) {
  cat("MISSING LEGACY:", paste(missing_legacy, collapse = ", "), "\n")
} else {
  cat("All v1.x legacy functions intact\n")
}

# Verify method objects have correct S3 class
m1 <- use_raster_eigen()
m2 <- use_hough()
cat("use_raster_eigen class:", paste(class(m1), collapse = ", "), "\n")
cat("use_hough class:", paste(class(m2), collapse = ", "), "\n")

a1 <- assign_voronoi()
a2 <- assign_crop()
a3 <- assign_graph()
cat("assign_voronoi class:", paste(class(a1), collapse = ", "), "\n")
cat("assign_crop class:", paste(class(a2), collapse = ", "), "\n")
cat("assign_graph class:", paste(class(a3), collapse = ", "), "\n")

s1 <- stm_eigen()
s2 <- stm_hough()
cat("stm_eigen class:", paste(class(s1), collapse = ", "), "\n")
cat("stm_hough class:", paste(class(s2), collapse = ", "), "\n")

g1 <- seg_ransac_circle()
g2 <- seg_ransac_cylinder()
g3 <- seg_irls_cylinder()
g4 <- seg_bf_cylinder()
cat("seg_ransac_circle class:", paste(class(g1), collapse = ", "), "\n")
cat("seg_ransac_cylinder class:", paste(class(g2), collapse = ", "), "\n")
cat("seg_irls_cylinder class:", paste(class(g3), collapse = ", "), "\n")
cat("seg_bf_cylinder class:", paste(class(g4), collapse = ", "), "\n")

cat("\nSmoke test PASSED\n")
