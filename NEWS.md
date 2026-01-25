# spanner 1.0.2

## CRAN Resubmission (2026-01-24)

* **Modernized spatial dependencies for CRAN compliance:**
    * Removed deprecated `sp` package completely from all dependencies
    * Eliminated `raster` package dependency by using modern lidR 4.x API
    * Replaced all `grid_metrics()` calls with `pixel_metrics()` for native terra support
    * All rasterization now returns `terra::SpatRaster` objects exclusively
    * Added `.onLoad()` hook to set `options(lidR.raster.default = "terra")`
    
* **Updated CRS handling to use sf ecosystem:**
    * Replaced `lidR::crs()` with `sf::st_crs()` throughout codebase
    * Updated all CRS assignments to use `sf::st_crs()` wrapper functions
    * Fixed CRS conversions to work without sp package

* **Enhanced documentation examples for CRAN:**
    * Added namespace prefixes (`lidR::`, `sf::`, `terra::`) to all function calls
    * Added safety checks for sparse datasets (null/empty data validation)
    * Fixed plot functions to use `terra::plot()` for SpatRaster objects
    * Fixed variable name typo in sum_rasters_by_suitability example
    * All examples now run successfully without sp or raster packages

* **Added comprehensive test suite:**
    * Implemented 60+ tests using testthat (edition 3)
    * Tests cover all core functionality: eigen metrics, cylinder fitting, tree processing, segmentation, PatchMorph algorithms
    * Includes edge case testing and error handling validation
    * Test suite designed for CRAN compliance (uses `skip_on_cran()` for intensive tests)

* **Bug fixes:**
    * Fixed `segment_graph()` null pointer error when no trees detected
    * Replaced `1:nrow()` with `seq_len(nrow())` to prevent zero-length errors
    * Added early return when tree.locations is NULL or empty

* **Performance validation:**
    * R CMD check --as-cran: Status OK (2 NOTEs only)
    * Total check time: 4.38 minutes (under 10-minute CRAN limit)
    * All unwrapped examples: < 1 second each (under 5-second CRAN limit)
    * Tests: 30 PASS, 16 SKIP on CRAN, 0 FAIL

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
