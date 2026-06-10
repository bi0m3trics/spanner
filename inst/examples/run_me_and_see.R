library(spanner)
library(lidR)
library(copc4R)

# LASfile <- system.file("extdata", "ALS_Clip.laz", package = "spanner")
# las <- readLAS(LASfile, select = "xyz")

# aoi <- get_aoi()
# tiles <- download_3dep_copc(aoi, dest_dir = "tiles/")
# las <- read_copc(tiles$file[1])
las <- read_copc("C:/GitHubRepos/spanner/tiles/3dep_copc_merge.copc.laz", progress = FALSE)
las <- as_las(las)

las <- classify_ground(las, ptd())
las <- normalize_height(las, tin())
# las <- filter_poi(las, Classification != 2 & Z >= 0.5)

library(future)
plan(multisession)

ctg <- readLAScatalog("C:/GitHubRepos/spanner/tiles/3dep_copc_merge.copc.laz")
opt_chunk_size(ctg) <- 300
opt_chunk_buffer(ctg) <- 20
opt_chunk_alignment(ctg) <- c(0, 0)
opt_progress(ctg) <- TRUE

my_catalog_segmentation <- function(chunk)
{
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)

  trees <- segment_trees(
    las,
    allometric_li_geodesic(
      seeds = NULL,
      hmin = 2,
      k = 12,
      max_jump = 1.5,
      bin_height = 1,
      crown_a = 0.35,
      crown_b = 0.65,
      max_crown_radius = NULL,
      center_drift = 1.5,
      density_radius = 1,
      gap_penalty = 3,
      geodesic_threshold = 5,
      use_connected_components = TRUE,
      min_component_size = 10,
      crown_profile = c("beta", "parabolic", "cone"),
      verbose = FALSE
    ))

  metrics <- crown_metrics(las = trees,
                           func = .stdtreemetrics,
                           geom = "point",
                           concaveman = c(3, 0),
                           attribute = "treeID")

  bbox  <- st_bbox(chunk)
  metrics <- sf::st_crop(metrics, bbox)
  return(metrics)
}

opt    <- list(need_buffer = TRUE,   # catalog_apply will throw an error if buffer = 0
               automerge   = TRUE)   # catalog_apply will merge the outputs into a single object
output <- catalog_apply(ctg, my_catalog_segmentation, .options = opt)


# allometric_li_geodesic
las_geodesic <- segment_trees(
  las,
  allometric_li_geodesic(
    seeds = NULL,
    hmin = 2,
    k = 12,
    max_jump = 1.5,
    bin_height = 1,
    crown_a = 0.35,
    crown_b = 0.65,
    max_crown_radius = NULL,
    center_drift = 1.5,
    density_radius = 1,
    gap_penalty = 3,
    geodesic_threshold = 5,
    use_connected_components = TRUE,
    min_component_size = 10,
    crown_profile = c("beta", "parabolic", "cone"),
    verbose = TRUE
  )
)
lidR::plot(las_geodesic, color = "treeID")

# allometric_random_walker
las_rw <- segment_trees(
  las,
  allometric_random_walker(
    seeds = NULL,
    hmin = 2,
    k = 12,
    alpha = 1,
    beta = 0.5,
    gamma = 0.5,
    delta = 1.5,
    eta = 1,
    crown_a = 0.35,
    crown_b = 0.65,
    eigen_radius = 1,
    density_radius = 1,
    probability_threshold = 0.35,
    max_iterations = 150,
    tolerance = 1e-05,
    crown_profile = c("beta", "parabolic", "cone"),
    verbose = TRUE
  )
)
lidR::plot(las_rw, color = "treeID")

# allometric_supervoxel_segment (tuned defaults)
las_sv <- segment_trees(
  las,
  allometric_supervoxel_segment(
    seeds = NULL,
    hmin = 2,
    voxel_size = 1.2,
    k_voxel = 18,
    min_voxel_points = 3,
    crown_a = 0.35,
    crown_b = 0.65,
    component_gap = 1.5,
    merge_threshold = 6,
    anisotropy_weight = 0.5,
    verticality_weight = 0.75,
    density_weight = 1,
    allometry_weight = 1.9,
    crown_profile = c("beta", "parabolic", "cone"),
    verbose = TRUE
  )
)
lidR::plot(las_sv, color = "treeID")


# las_sv <- spanner::colorize_las(las_sv, method = "attr", attribute_name = "treeID", palette = random.colors(500))
# lidRviewer::view(las_sv)

message("Done. All algorithms ran successfully.")
