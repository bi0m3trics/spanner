library(spanner)
library(lidR)
library(lidRviewer)

las <- readLAS("C:/Users/ajsm/Desktop/Temp/Setup49.las")

myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.025, pt_spacing = 0.0254,
                                       dens_threshold = 0.25,
                                       neigh_sizes = c(0.25, 0.15, 0.66),
                                       eigen_threshold = 0.75,
                                       grid_slice_min = 1,
                                       grid_slice_max = 2,
                                       minimum_polygon_area = 0.005,
                                       cylinder_fit_type = "ransac",
                                       max_dia = 1,
                                       SDvert = 0.33,
                                       n_pts = 20,
                                       n_best = 25,
                                       inliers = 0.9,
                                       conf = 0.99,
                                       max_angle = 20)


myTreeGraph = segment_graph(las = las, tree.locations = myTreeLocs, k = 50,
                            distance.threshold = 0.5,
                            use.metabolic.scale = TRUE,
                            metabolic.scale.function = '1/((2*x)^(1/8))',
                            ptcloud_slice_min = 1,
                            ptcloud_slice_max = 2,
                            subsample.graph = 0.1,
                            return.dense = TRUE)

bole_las <- classify_stem_points(myTreeGraph,
                                myTreeLocs,
                                method = "lewos",
                                z_min = 0.05,
                                z_max = 30,
                                search_radius = 1,
                                neigh_k = 30L,
                                linearity_threshold = 0, # ignores Linearity
                                verticality_threshold = 0.8,
                                n_propagation = 10L,
                                k_propagation = 30L,
                                hough_min_radius = 0.01,
                                hough_max_radius = 1,
                                hough_min_votes = 30L,
                                ncpu = 16L,
                                voxel_thin = NULL
)

# plot(bole_las, color = "WoodProb")

seg_tbl  <- segment_bole(bole_las,                         # Takes a minute!
                         myTreeLocs,
                         algorithm = "hough",
                         z_min = 0.1,
                         z_max = 30,
                         dz = 1.0,
                         overlap = 0.1,
                         use_stem_column = TRUE,
                         search_radius = NULL,
                         inlier_tol = 0.03,
                         n_ransac = 30L,
                         conf = 0.99,
                         inlier_frac = 0.85,
                         n_best = 10L,
                         max_radius = 1,
                         min_pts = 10L,
                         hough_min_radius = 0.02,
                         hough_max_radius = 1,
                         hough_min_votes = 15L,
                         voxel_thin = NULL,
                         ncpu = 16L  )
print(seg_tbl)

vol_tbl  <- compute_bole_volume(seg_table = seg_tbl,
                                method = "smalian")
print(vol_tbl)

qual_tbl <- assess_tree_quality(myTreeLocs,
                                seg_tbl,
                                vol_table = vol_tbl,
                                max_mean_error = 0.02,
                                max_max_error = 0.05,
                                min_inlier_frac = 0.7,
                                min_bole_coverage = 0.5,
                                dbh_ht_ratio_range = c(0.004, 0.05),
                                max_radius_cv = 0.4,
                                min_segments = 3L)
table(qual_tbl$quality_class)

refit_result <- refit_trees(myTreeGraph,
                            myTreeLocs,
                            seg_tbl,
                            qual_tbl,
                            strategies = "meanshift_bole",
                            refit_marginal = FALSE,
                            max_iterations = 3L,
                            z_min = 0.1,
                            z_max = NULL,
                            dz = 0.5,
                            overlap = 0.1,
                            search_radius = NULL,
                            inlier_tol = 0.03,
                            n_ransac = 20L,
                            conf = 0.99,
                            inlier_frac = 0.85,
                            n_best = 25L,
                            max_radius = 1,
                            min_pts = 10L,
                            relax_delta = 0.15,
                            dbscan_eps_factor = 2,
                            dbscan_min_pts = 5L,
                            ms_bandwidth = NULL,
                            ms_max_iter = 30L,
                            ms_tol = 0.001,
                            ncpu = 16L)

print(refit_result$refit_log)

myTreeGraph = segment_graph(las = las, tree.locations = refit_result$seg_table,
                            k = 50,
                            distance.threshold = 0.5,
                            use.metabolic.scale = TRUE,
                            metabolic.scale.function = '1/((2*x)^(1/8))',
                            ptcloud_slice_min = 1,
                            ptcloud_slice_max = 2,
                            subsample.graph = 0.1,
                            return.dense = TRUE)

processed_data <- process_tree_data(treeData = refit_result$seg_table,
                                    segmentedLAS = myTreeGraph,
                                    return_sf = TRUE,
                                    seg_table = seg_tbl,
                                    vol_table = vol_tbl,
                                    qual_table = qual_tbl)

myTreeGraph <- spanner::colorize_las(myTreeGraph, method = "attr",
                                     attribute_name = "treeID",
                                     palette = random.colors(500))
lidRviewer::view(myTreeGraph)
