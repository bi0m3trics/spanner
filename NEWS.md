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
