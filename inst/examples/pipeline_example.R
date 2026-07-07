# spanner v2 TLS/MLS staged workflow example
#
# This example uses the MLS_Clip.laz file shipped in inst/extdata and runs:
# find_trees() -> tree_points() -> stem_points() -> stem_segments() ->
# tree_inventory()

library(spanner)
library(lidR)
library(sf)

mls_file <- system.file("extdata", "MLS_Clip.laz", package = "spanner")
las <- lidR::readLAS(mls_file, select = "xyzcr")
if (lidR::is.empty(las)) stop("MLS_Clip.laz could not be read.")

# Set the CRS used by the package example data when it is not already present.
sf::st_crs(las) <- 26912

# Stage 1: detect tree seed locations.
tree_map <- find_trees(
  las,
  method = use_hough(
    min_h = 1.0,
    max_h = 3.0,
    h_step = 0.5,
    max_d = 0.8,
    min_density = 0.05,
    min_votes = 3L
  )
)

# Stage 2: assign points to detected trees.
las_tree <- tree_points(
  las,
  tree_map,
  method = assign_voronoi(max_dist = 6)
)

# Stage 3: classify likely stem points.
las_stem <- stem_points(
  las_tree,
  tree_map,
  method = stem_hough(
    h_base_lo = 1.0,
    h_base_hi = 2.5,
    h_step = 0.5,
    max_d = 0.8,
    min_density = 0.05,
    min_votes = 3L
  )
)

# Stage 4: fit per-tree stem segments.
stem_seg <- stem_segments(
  las_stem,
  tree_map,
  method = seg_ransac_cylinder(
    dz = 0.5,
    overlap = 0.1,
    max_radius = 0.8,
    min_pts = 10L
  )
)

# Optional summaries for quality assessment and inventory output.
stem_vol <- compute_stem_volume(stem_seg, method = "cylinder")
stem_quality <- assess_tree_quality(
  tree_map,
  stem_seg,
  vol_table = stem_vol,
  expected_stem_height = 3.0
)

# Stage 5: create a tree-level inventory.
inventory <- tree_inventory(
  tree_map,
  las_tree,
  seg_table = stem_seg,
  vol_table = stem_vol,
  qual_table = stem_quality
)

print(inventory)
print(table(stem_quality$quality_class))
