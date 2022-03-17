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
