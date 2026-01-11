library(spanner)
library(lidRviewer)
library(lidR)

# Load your lidar data and prep it dfor processing
las <- readLAS("E:/GithubRepos/spanner - Copy - Copy/inst/extdata/SparsePatchA.laz")
projection(las) = sp::CRS("+init=epsg:26912")
las = classify_ground(las, csf(sloop_smooth = FALSE,
                               class_threshold = 0.5,
                               cloth_resolution = 0.5, rigidness = 1L,
                               iterations = 500L, time_step = 0.65))
las = normalize_height(las, tin())
las = classify_noise(las, ivf(0.25, 3))
las = filter_poi(las, Classification != LASNOISE)

# Run comprehensive tree detection
trees <- comprehensive_tree_detection(las, neigh_sizes = c(0.5, 0.25, 0.75),
                                      res = 0.0254, grid_slice_min = 1.2,
                                      grid_slice_max = 1.8, dens_threshold = 0.15,
                                      vert_threshold = 0.6, minimum_area = 0.1,
                                      pt_spacing = 0.025, SDvert = 0.15,
                                      ransac_error_threshold = 0.15,
                                      cylinder_search_radius = 0.4,
                                      min_cylinder_points = 15,
                                      n_ransac_samples = 10, ransac_confidence = 0.95,
                                      ransac_inliers = 0.8, max_angle = 20,
                                      n_best = 20,
                                      second_pass_vert_threshold = 0.75,
                                      second_pass_aniso_threshold = 0.7,
                                      ncpu = parallel::detectCores(logical = FALSE),
                                      verbose = TRUE)

validate_detected_trees(trees, min_radius = 0.05, max_radius = 1,
                        max_error = 0.2, min_points = 10)

# trees is now an sf object with tree locations, radii, and quality metrics
plot_tree_detection(las, trees, color_by = "Quality", point_size = 3)

# Segement point cloud
treegraph = segment_graph(las = las, tree.locations = trees, k = 50,
                          distance.threshold = 0.5, use.metabolic.scale = FALSE,
                          ptcloud_slice_min = 0.6666, ptcloud_slice_max = 2.0,
                          subsample.graph = 0.1, return.dense = FALSE)

# View it
las_colored <- colorizeLAS(treegraph, attribute = "treeID", palette_func = lidR::random.colors)

view(las_colored)
