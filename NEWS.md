# spanner 2.0.0

## New modular pipeline

spanner 2.0 introduces a five-stage pipeline API that replaces the previous
collection of stand-alone functions.  Each stage dispatches on a swappable
*method object*, making it easy to mix and match detection and fitting
strategies without rewriting pipeline code.

```
find_trees()      → tree_points()  → stem_points()
stem_segments()   → tree_inventory()
```

### Stage 1 — `find_trees(las, method)`

Detects individual tree locations and returns an `sf` POINT table
(`TreeID`, `X`, `Y`, `Radius`, …).

* `use_raster_eigen()` – wraps the existing `get_raster_eigen_treelocs()`
  raster-eigenvalue detector.  Used as `find_trees(las, use_raster_eigen())`.
* `use_hough()` – new Hough-accumulator stem detector.  Stacks horizontal
  density slices across a user-defined height range, accumulates votes in
  a 2-D pixel grid, and clusters peaks into tree centres.  Works on both
  TLS and MLS without pre-classified ground.

### Stage 2 — `tree_points(las, map, method)`

Attributes every point to a tree ID.

* `assign_graph()` – wraps the existing `segment_graph()` minimum-spanning-
  tree segmentation.
* `assign_voronoi()` – fast Voronoi-based assignment; each point is claimed
  by its nearest tree seed within an optional maximum radius.
* `assign_crop()` – cylindrical crop around each tree seed.

### Stage 3 — `stem_points(las, map, method)`

Labels stem/bole points within the segmented cloud.

* `stem_eigen()` – two-scale eigenvalue leaf-wood classifier (previously
  `classify_lw_points()`); adds `Stem` and `StemProb` columns.
* `stem_hough()` – new per-tree Hough-based stem detector; labels points that
  lie on circular cross-sections fitted by the Hough accumulator.

### Stage 4 — `stem_segments(las, map, method)`

Fits 3-D cylinders across overlapping height slices for every tree and
returns a tidy segment table (`TreeID`, `Segment`, `Z_low`, `Z_high`,
`Z_mid`, `Radius`, `Diameter`, `Fit_method`, `N_pts`, `Valid`).

* `seg_ransac_cylinder()` – RANSAC 3-D cylinder fit.  All trees × all slices
  are processed in a single vectorised C++ call, making this 50–100× faster
  than the previous per-slice R loop.
* `seg_irls_cylinder()` – IRLS (iteratively reweighted least squares) 3-D
  cylinder fit with the same performance characteristics.

### Stage 5 — `tree_inventory(map, segmented_las, seg_table, …)`

Rolls up all pipeline outputs into a per-tree `sf` inventory table.
Computes crown geometry, merges segment and quality tables, and optionally
returns convex-hull crown polygons.

## New functions

* `find_trees()` – stage 1 dispatcher; see pipeline section above.
* `tree_points()` – stage 2 dispatcher.
* `stem_points()` – stage 3 dispatcher.
* `stem_segments()` – stage 4 dispatcher.
* `tree_inventory()` – stage 5 dispatcher.
* `use_raster_eigen()`, `use_hough()` – method factories for `find_trees()`.
* `assign_voronoi()`, `assign_crop()`, `assign_graph()` – method factories for `tree_points()`.
* `stem_eigen()`, `stem_hough()` – method factories for `stem_points()`.
* `seg_ransac_cylinder()`, `seg_irls_cylinder()` – method factories for `stem_segments()`.

* `compute_stem_volume()` – (renamed from `compute_bole_volume()`) aggregates
  a segment table from `stem_segments()` into per-tree outside-bark stem
  volume.  Accepts any table with `TreeID`, `Segment`, `Z_low`, `Z_high`,
  `Radius`, `Valid`; RANSAC diagnostic columns are optional.

* `allometric_li_geodesic()` – builds a `lidR::segment_trees()`-compatible
  algorithm using allometric crown envelopes and geodesic growing on a local
  kNN graph.

* `allometric_random_walker()` – seeded diffusion-style crown segmentation
  algorithm for `lidR::segment_trees()` with allometric penalties.

* `allometric_supervoxel_segment()` – supervoxel-graph segmentation for
  `lidR::segment_trees()` with backfilling and small-island cleanup.

* `compute_lai()` – estimates leaf area index (LAI) and leaf area density
  (LAD) per grid cell from any height-normalized lidar point cloud using the
  MacArthur-Horn gap-fraction Beer-Lambert inversion.  Accepts a single `LAS`
  or a `LAScatalog` (parallelised over chunks).  Returns a three-layer
  `SpatRaster`: `LAI`, `LAD_mean`, `LAD_max`.

* `compute_pad_voxels()` – estimates 3-D plant area density (PAD, m²/m³)
  by tracing lidar pulses downward through a regular voxel grid and applying
  the Beer-Lambert inversion.  Per-voxel output includes `directed`,
  `transmitted`, `intercepted`, `occluded`, `PAD`, and `VoxelClass`.

* `compute_transmittance_raster()` – converts the per-voxel PAD table from
  `compute_pad_voxels()` into a 2-D `SpatRaster` of Beer-Lambert canopy
  transmittance (τ ∈ [0, 1]) for each (X, Y) column.

* `branch_metrics()` – scores and flags branch-candidate points from
  `eigen_metrics()` output.  Adds `AxisAngle`, `BranchScore`, and
  `IsBranchCandidate` to the input `data.table` in-place.

* `classify_lw_points()` – two-scale eigenvalue leaf-wood classifier; adds
  `Bole`, `BoleProb`, `Branch`, `BranchProb`, and `Other` columns.  Used
  internally by `stem_eigen()` and available as a standalone function.

* `refit_trees()` – re-runs stem fitting with relaxed or tightened parameters
  for trees flagged as `"bad"` or `"marginal"` by `assess_tree_quality()`.
  Supports multiple sequential refitting strategies and handles `NA` tree IDs
  without crashing.

## Changes to existing functions

* `process_tree_data()` – when `seg_table` is supplied and `vol_table` is
  `NULL`, stem volume is now computed automatically via `compute_stem_volume()`
  and merged into the result.

* `assess_tree_quality()` – now accepts segment tables from `stem_segments()`
  in addition to `segment_bole()`.  `RANSAC_err` and `Inlier_frac` columns
  are optional; quality-score components that depend on them default to zero
  when absent.

* `eigen_metrics()` – now returns three additional columns: `E1x`, `E1y`,
  `E1z` (components of the dominant eigenvector λ₁), required by
  `branch_metrics()`.

## Deprecated functions

The following functions remain exported for backwards compatibility but are
superseded by the v2 pipeline and may be removed in a future release.

* `segment_bole()` – use `stem_segments()` instead.
* `classify_stem_points()` – use `stem_points(method = stem_eigen())` or
  `stem_points(method = stem_hough())` instead.

## Internal C++ improvements

* New Hough accumulator for stem detection (`C_stack_map`,
  `C_hough_stem_plot`): processes all height slices in a single pass over the
  point cloud with OpenMP parallelism.
* New batch cylinder fitters (`C_ransac_plot_cylinders`,
  `C_irls_plot_cylinders`): accept all trees and all slices in one call,
  eliminating the per-slice R dispatch overhead that dominated runtime in v1.
* New Voronoi and crop tree-ID assigners (`C_tree_ids_voronoi`,
  `C_tree_ids_crop`) implemented directly in C++.

# spanner 1.0.4

* Fixed non-API call to R (`R_UnboundValue`) in compiled code (flagged by CRAN r-devel checks on Linux and Windows).
* Fixed `create_rotation_gif()` example to use `if (interactive())` inside `\donttest{}` so they are skipped on headless CI environments.

# spanner 1.0.3

* Replaced bare `Rf_error()` call in `src/backports.h` with the parenthesized form `(Rf_error)(...)` to avoid interception by future Rcpp macro wrapping (see RcppCore/Rcpp#1247).

# spanner 1.0.2

* Added `screen_size` parameter to `create_rotation_gif()` to allow custom window dimensions (e.g., `c(800, 600)` or `c(1920, 1080)`)

***CRAN Release for v1.0.2 on Feb 3rd, 2026***
  
* Added `colorize_las()` function with four coloring methods: `attr` (attribute-based), `rgb` (raster extraction), `pcv` (true 3D ambient occlusion), and `ssao` (fast screen-space ambient occlusion)
* Added `download_naip_for_las()` function to automatically download NAIP imagery from Microsoft Planetary Computer STAC API
* Added C++ implementations with OpenMP parallelization for PCV and SSAO ambient occlusion methods
* Expanded comprehensive test suite for all new functionality
* Removed deprecated `sp` and `raster` package dependencies for CRAN compliance
* Replaced all `grid_metrics()` calls with `pixel_metrics()` for native terra support
* Updated CRS handling to use `sf::st_crs()` throughout codebase
* Fixed documentation examples to include proper namespace prefixes and run without deprecated packages
    * Includes edge case testing and error handling validation
    * Test suite designed for CRAN compliance (uses `skip_on_cran()` for intensive tests)

* **Bug fixes:**
    * Fixed `segment_graph()` scientific notation issue where large point indices (e.g., 800000) were converted to strings like '8e+05' by `cppRouting::makegraph()`, causing "not all nodes are in the graph" errors. Now explicitly formats indices without scientific notation.
    * Fixed `segment_graph()` null pointer error when no trees detected
    * Replaced `1:nrow()` with `seq_len(nrow())` to prevent zero-length errors
    * Added early return when tree.locations is NULL or empty

## Previous Changes in 1.0.2

* New eigen metrics added (to match CloudCompare)
    * Roughness: Distance from query point to fitted plane through neighborhood centroid
    * Mean Curvature: Differential geometry-based curvature using quadric surface fitting
    * Gaussian Curvature: Product of principal curvatures from quadric surface
    * PCA1: Eigenvector projection variance normalized by eigensum
    * PCA2: Second eigenvector projection variance normalized by eigensum
    * NumNeighbors: Count of points in sphere neighborhood
* Code optimization and cleanup:
    * Removed 20+ unused C++ exports to streamline the package interface
    * Fixed namespace issues by removing ::: calls in internal code
    * Reduced exported C++ functions to only those actively used (C_eigen_in_sphere, C_count_in_disc, C_count_in_sphere, cppCylinderFit)
    * Removed lidR from Imports field (kept in Depends and LinkingTo only)
* Updated documentation examples:
    * Fixed examples to use sf::st_coordinates() with proper namespace qualification
    * Updated get_raster_eigen_treelocs() examples with optimized parameters for better tree detection in forests with interlocking crowns (res=0.25, dens_threshold=0.25, eigen_threshold=0.75, minimum_polygon_area=0.005)
    * Corrected circle area formula in plotting examples (Radius^2*3.14)
* Bug fixes - replaced null checks with is.empty for LAS objects and stopped R from collapsing the one-row subset into a vector, so the circle fit still receives a 2-column input when there's on one tree in `get_raster_eigen_treelocs`.
* Removed depends on magrittr and removed all %>% in codebase
* Added `process_tree_data` funciton that takes the output of `get_raster_eigen_treelocs` and `segment_graph` to adds information about the height, crown area and volume, and diameter for each unique TreeID. It also has an optional parameter to return either points or hulls as an `sf` object for each tree.
* Added the citation for the package
* Added a couple default datasets and got rid of getExampleData()
* Added the xyz normals as returns for eigen_metrics()
* Added PatchMorph functions:
    * process_rasters_patchmorph: Processes an input raster by reclassifying it based on suitability levels and applying gap and spur distance transformations to generate a list of processed rasters.
    * plot_raster_by_name: Plots a raster from a list of rasters based on the provided raster name.
    * sum_rasters_by_suitability: Sums rasters from a list based on their suitability levels and returns a list of summed rasters for each suitability level.
* Added dependencies fopr sf and terra
* Modified get_raster_eigen_treelocs and segment_graph to use sf and not write any intermediate files to optput locations with parallel processing to make sure that all possible operations use
    * the available CPU cores
    *efficient data structures; used lapply for list operations and dplyr::bind_rows for combining data frames.
    *reduce redundant ralculations; stored intermediate results and reused them where possible.
    *removed unnecessary objects and used more efficient data structures.
* Added spanner_pal() which is a custom color palette


# spanner 1.0.1

* Dependencies for raster now replaced by terra and calls to sp replaced with sf (in spanner only).
* Consumed the ransac cylinder fitting code from TreeLS and removed it as a depends/import. This code had to be modified slightly to bring it in line with R for 4.1.2.
* Updated get_raster_eigen_treelocs to better handle tree cover when smaller slices are needed.
* Updated the readme file to include all relevant literature.
* Fixed a bug resulting in treeID's from segment_graph not matching those created in get_raster_eigen_treelocs. Added parameters to specify where the point cloud slice should be used for matching resulting IDs (so they correspond to those provided from get_raster_eigen_treelocs) 

# spanner 1.0.0

* Manuscript release version

# spanner 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
