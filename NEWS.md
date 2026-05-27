# spanner 1.1.0

## New functions

* `compute_lai()` – estimates leaf area index (LAI) and leaf area density
  (LAD) per grid cell from any height-normalized lidar point cloud using the
  MacArthur-Horn gap-fraction Beer-Lambert inversion ([`lidR::LAD()`] under
  the hood). Accepts a single `LAS` or a `LAScatalog` (parallelised over
  chunks). Returns a three-layer `SpatRaster`: `LAI`, `LAD_mean`, `LAD_max`.

* `compute_pad_voxels()` – estimates 3-D plant area density (PAD, m²/m³)
  by tracing lidar pulses downward through a regular voxel grid and applying
  the Beer-Lambert inversion following Grau et al. (2017) and the voxelmon
  framework (Tenny et al. 2025). For each voxel the function reports
  `directed`, `transmitted`, `intercepted`, `occluded`, `PAD`, and a
  `VoxelClass` label (`"foliage"`, `"wood"`, `"empty"`, `"occluded"`). LAI
  per column and vertical height profiles can be derived directly from the
  returned `data.table`.

* `compute_transmittance_raster()` – converts the per-voxel PAD table from
  `compute_pad_voxels()` into a 2-D `SpatRaster` of Beer-Lambert canopy
  transmittance (τ ∈ [0, 1]) for each (X, Y) column:
  τ = exp(−G × Σ min(PAD_i, max_pad) × vox_size). Per-voxel PAD is capped
  at `max_pad` before summation (same default as `compute_pad_voxels()`) to
  prevent numerically extreme wood-voxel values from collapsing transmittance
  to zero. The returned raster can be multiplied directly against a vostokR
  `solar_potential` surface to estimate below-canopy irradiance.

* `branch_metrics()` – scores and flags branch-candidate points from the
  output of `eigen_metrics()`. Adds `AxisAngle` (°, angle of the dominant
  local axis from vertical), `BranchScore` (min–max scaled composite), and
  `IsBranchCandidate` (logical) to the input `data.table` in-place.

* `classify_lw_points()` – combined leaf-wood classifier (bole + branch +
  other in a single two-scale eigen pass). Adds `Bole`, `BoleProb`,
  `Branch`, `BranchProb`, and `Other` columns to the returned LAS.

* `segment_bole()` – fits RANSAC circles to bole points per tree and returns
  a per-tree segment table with diameter, height, and fit-quality statistics.
  Accepts a `ncpu` parameter to parallelise fitting over trees via OpenMP.
  Prefers the `Bole` column produced by `classify_lw_points()` over the
  legacy `Stem` column from `classify_stem_points()` when both are present.

* `refit_trees()` – re-runs `segment_bole()` with relaxed or tightened
  parameters for trees flagged as `"bad"` or `"marginal"` by
  `assess_tree_quality()`. Handles `NA` tree IDs in the segmented point
  cloud without crashing.

## New features in existing functions

* `eigen_metrics()` – now returns three additional columns: `E1x`, `E1y`,
  `E1z` (components of the dominant eigenvector λ₁). Required by
  `branch_metrics()` and available for any orientation-based analysis.

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
