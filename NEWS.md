# spanner 1.1.0

## Key Updates

### Major Performance and Functionality Enhancements
- **process_rasters_patchmorph()** fully rewritten for speed and efficiency, with:
  - Kernel caching (40–60% faster)
  - Improved memory handling
  - Robust error recovery
  - Thorough parameter validation
- **getCircularKernel()** optimized with vectorized operations for 2–3× speed gains
- Reduced memory footprint and external dependencies
- Standardized coding style and improved error messages across functions

### New and Improved Visualization Tools
- **plot_raster_by_name()** redesigned with ggplot2 integration:
  - Multiple palettes and custom colors
  - Professional themes and flexible layouts
  - Fallback to base plotting if ggplot2 unavailable
- **plot_multiple_rasters()** for side-by-side raster comparisons with synchronized scales and consistent styling
- **plot_raster_interactive()** with Plotly for zooming, panning, hover info, and HTML export
- Removed rarely used contour functionality to streamline plotting

### Expanded Tree Detection Capabilities
- **comprehensive_tree_detection()**: End-to-end workflow combining raster-based detection, RANSAC cylinder fitting, and geometric filtering, returning sf outputs with quality classifications
- **validate_detected_trees()**: Filters detections by radius, error, and point count
- **plot_tree_detection()**: Visualizes detected trees over CHMs with size and color scaling
- **colorizeLAS()**: Colors LAS point clouds by any attribute with customizable palettes and NA handling

### Codebase and Documentation Improvements
- Cleaned repository of unused code, artifacts, and redundant examples
- Focused package scope on core tree detection
- Comprehensive roxygen2 documentation updates
- Expanded examples and added `simple_example.R` demonstrating complete workflow

### Bug Fixes
- Corrected NAMESPACE exports
- Improved cylinder fitting error handling
- Enhanced geometric feature integration
- Better handling of edge cases in detection pipeline


# spanner 1.0.2
* Added the citation for the package
* Added a couple default datasets and fixed the call for getExampleData()
* Added the xyz normals as returns for eigen_metrics()
* Added PatchMorph functions:
** process_rasters_patchmorph: Processes an input raster by reclassifying it based on suitability levels and applying gap and spur distance transformations to generate a list of processed rasters.
** plot_raster_by_name: Plots a raster from a list of rasters based on the provided raster name.
** sum_rasters_by_suitability: Sums rasters from a list based on their suitability levels and returns a list of summed rasters for each suitability level.
* Added dependencies fopr sf and terra
* Modified get_raster_eigen_treelocs and segment_graph
** to use sf and not write any intermediate files to optput locations
** wit parallel processing to make sure that all possible operations use the available CPU cores.
** efficient data structures; used lapply for list operations and dplyr::bind_rows for combining data frames.
** reduce redundant ralculations; stored intermediate results and reused them where possible.
** removed unnecessary objects and used more efficient data structures.
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
