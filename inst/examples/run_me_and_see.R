library(spanner)
library(lidR)
library(copc4R)

# LASfile <- system.file("extdata", "ALS_Clip.laz", package = "spanner")
# las <- readLAS(LASfile, select = "xyz")

# aoi <- get_aoi()
# tiles <- download_3dep_copc(aoi, dest_dir = "tiles/")
# las <- read_copc(tiles$file[1])
las <- read_copc("C:/GitHubRepos/spanner/tiles/3dep_copc_merge.copc.laz",
                 progress = FALSE)
las <- as_las(las)

las <- classify_ground(las, ptd())
las <- normalize_height(las, tin())
# las <- filter_poi(las, Classification != 2 & Z >= 0.5)
# writeLAS(las, "C:/GitHubRepos/spanner/tiles/3dep_copc_merge.laz")

# ==============================================================================
# ==============================================================================

library(future)
plan(multisession, workers = parallelly::availableCores(logical = FALSE) - 1)
set_lidr_threads(1L)

ctg <- readLAScatalog("C:/GitHubRepos/spanner/tiles/3dep_copc_merge.laz")
opt_chunk_size(ctg) <- 500
opt_chunk_buffer(ctg) <- 20
opt_chunk_alignment(ctg) <- c(0, 0)
opt_progress(ctg) <- TRUE

my_catalog_allometric_li_geodesic <- function(chunk)
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

  points <- crown_metrics(las = trees,
                           func = .stdtreemetrics,
                           geom = "point",
                           concaveman = c(3, 0),
                           attribute = "treeID")
  crowns <- crown_metrics(las = trees,
                          func = .stdtreemetrics,
                          geom = "convex",
                          concaveman = c(3, 0),
                          attribute = "treeID")
  chm <- rasterize_canopy(las = trees,
                          res = 0.5,
                          algorithm = spikefree(1.5))
  # chm = pitfill_stonge2008(chm) # Pits and spikes filling

  bbox  <- st_bbox(chunk)
  points <- sf::st_crop(points, bbox)

  crowns <- dplyr::semi_join(crowns,
                           dplyr::select(sf::st_drop_geometry(points), treeID),
                           by = "treeID")
  # chm <- terra::crop(chm, bbox)


  list(points, crowns, terra::wrap(chm))
}

opt    <- list(need_buffer = TRUE,   # catalog_apply will throw an error if buffer = 0
               automerge   = FALSE)  # keep as list of lists; merge manually after unwrapping terra objects
output_geodesic <- catalog_apply(ctg, my_catalog_allometric_li_geodesic, .options = opt)

extract_part <- function(x, idx) if (!is.null(x)) x[[idx]] else NULL
pts_list <- lapply(output_geodesic, extract_part, idx = 1)
crowns_list <- lapply(output_geodesic, extract_part, idx = 2)
chm_list <- lapply(output_geodesic, extract_part, idx = 3)
chm_list <- lapply(chm_list, terra::unwrap)

pts_list <- Filter(Negate(is.null), pts_list)
crowns_list <- Filter(Negate(is.null), crowns_list)
chm_list <- Filter(Negate(is.null), chm_list)

if (length(pts_list) == 0 || length(crowns_list) == 0 || length(chm_list) == 0) stop("No trees found. Check inputs and parameters.")

pts_all_geodesic <- do.call(rbind, pts_list)
crowns_all_geodesic <- do.call(rbind, crowns_list)
chm_geodesic <- do.call(terra::merge, chm_list)

# ==============================================================================

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

# ==============================================================================
# ==============================================================================

library(future)
plan(multisession, workers = parallelly::availableCores(logical = FALSE) - 1)
set_lidr_threads(1L)

ctg <- readLAScatalog("C:/GitHubRepos/spanner/tiles/3dep_copc_merge.laz")
opt_chunk_size(ctg) <- 500
opt_chunk_buffer(ctg) <- 20
opt_chunk_alignment(ctg) <- c(0, 0)
opt_progress(ctg) <- TRUE

my_catalog_segmentation_allometric_random_walker <- function(chunk)
{
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)

  seeds_rw <- lidR::locate_trees(
    las,
    lidR::lmf(ws = function(z) pmin(12, pmax(3, 0.1 * z + 3)), hmin = 2)
  )

  trees <- segment_trees(
    las,
    allometric_random_walker(
      seeds = seeds_rw,
      hmin = 2,
      k = 15,
      alpha = 1.2,
      beta = 0.5,
      gamma = 0.7,
      delta = 2.0,
      eta = 1,
      crown_a = 0.32,
      crown_b = 0.65,
      eigen_radius = 1,
      density_radius = 1,
      crown_profile = c("beta", "parabolic", "cone"),
      verbose = TRUE
    )
  )

  points <- crown_metrics(las = trees,
                          func = .stdtreemetrics,
                          geom = "point",
                          concaveman = c(3, 0),
                          attribute = "treeID")
  crowns <- crown_metrics(las = trees,
                          func = .stdtreemetrics,
                          geom = "convex",
                          concaveman = c(3, 0),
                          attribute = "treeID")
  chm <- rasterize_canopy(las = trees,
                          res = 0.5,
                          algorithm = spikefree(1.5))
  # chm = pitfill_stonge2008(chm) # Pits and spikes filling

  bbox  <- st_bbox(chunk)
  points <- sf::st_crop(points, bbox)

  crowns <- dplyr::semi_join(crowns,
                             dplyr::select(sf::st_drop_geometry(points), treeID),
                             by = "treeID")
  # chm <- terra::crop(chm, bbox)


  list(points, crowns, terra::wrap(chm))
}

opt    <- list(need_buffer = TRUE,   # catalog_apply will throw an error if buffer = 0
               automerge   = FALSE)  # keep as list of lists; merge manually after unwrapping terra objects
output_rw <- catalog_apply(ctg, my_catalog_segmentation_allometric_random_walker, .options = opt)

extract_part <- function(x, idx) if (!is.null(x)) x[[idx]] else NULL
pts_list <- lapply(output_rw, extract_part, idx = 1)
crowns_list <- lapply(output_rw, extract_part, idx = 2)
chm_list <- lapply(output_rw, extract_part, idx = 3)
chm_list <- lapply(chm_list, terra::unwrap)

pts_list <- Filter(Negate(is.null), pts_list)
crowns_list <- Filter(Negate(is.null), crowns_list)
chm_list <- Filter(Negate(is.null), chm_list)

if (length(pts_list) == 0 || length(crowns_list) == 0 || length(chm_list) == 0) stop("No trees found. Check inputs and parameters.")

pts_all_rw <- do.call(rbind, pts_list)
crowns_all_rw <- do.call(rbind, crowns_list)
chm_rw <- do.call(terra::merge, chm_list)

# ==============================================================================

# allometric_random_walker
# Pre-compute seeds with a wider LMF window to control seed count directly.
# Increase the coefficient (0.1) or the intercept (3) to get fewer seeds.
seeds_rw <- lidR::locate_trees(
  las,
  lidR::lmf(ws = function(z) pmin(12, pmax(2.7, 0.1 * z + 2.7)), hmin = 2)
)

las_rw <- segment_trees(
  las,
  allometric_random_walker(
    seeds = seeds_rw,
    hmin = 2,
    k = 15,
    alpha = 1.2,
    beta = 0.5,
    gamma = 0.7,
    delta = 2.0,
    eta = 1,
    crown_a = 0.32,
    crown_b = 0.65,
    eigen_radius = 1,
    density_radius = 1,
    crown_profile = c("beta", "parabolic", "cone"),
    verbose = TRUE
  )
)

lidR::plot(las_rw, color = "treeID")

las_rw <- spanner::colorize_las(las_rw, method = "attr",
                                attribute_name = "treeID",
                                palette = random.colors(500))
lidRviewer::view(las_rw)

# ==============================================================================
# ==============================================================================

library(future)
plan(multisession, workers = parallelly::availableCores(logical = FALSE) - 1)
set_lidr_threads(1L)

ctg <- readLAScatalog("C:/GitHubRepos/spanner/tiles/3dep_copc_merge.laz")
opt_chunk_size(ctg) <- 500
opt_chunk_buffer(ctg) <- 20
opt_chunk_alignment(ctg) <- c(0, 0)
opt_progress(ctg) <- TRUE

my_catalog_segmentation_allometric_supervoxel_segment <- function(chunk)
{
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)

  trees <- segment_trees(
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

  points <- crown_metrics(las = trees,
                          func = .stdtreemetrics,
                          geom = "point",
                          concaveman = c(3, 0),
                          attribute = "treeID")
  crowns <- crown_metrics(las = trees,
                          func = .stdtreemetrics,
                          geom = "convex",
                          concaveman = c(3, 0),
                          attribute = "treeID")
  chm <- rasterize_canopy(las = trees,
                          res = 0.5,
                          algorithm = spikefree(1.5))
  # chm = pitfill_stonge2008(chm) # Pits and spikes filling

  bbox  <- st_bbox(chunk)
  points <- sf::st_crop(points, bbox)

  crowns <- dplyr::semi_join(crowns,
                             dplyr::select(sf::st_drop_geometry(points), treeID),
                             by = "treeID")
  # chm <- terra::crop(chm, bbox)


  list(points, crowns, terra::wrap(chm))
}

opt    <- list(need_buffer = TRUE,   # catalog_apply will throw an error if buffer = 0
               automerge   = FALSE)  # keep as list of lists; merge manually after unwrapping terra objects
output_sv <- catalog_apply(ctg, my_catalog_segmentation_allometric_supervoxel_segment, .options = opt)

extract_part <- function(x, idx) if (!is.null(x)) x[[idx]] else NULL
pts_list <- lapply(output_sv, extract_part, idx = 1)
crowns_list <- lapply(output_sv, extract_part, idx = 2)
chm_list <- lapply(output_sv, extract_part, idx = 3)
chm_list <- lapply(chm_list, terra::unwrap)

pts_list <- Filter(Negate(is.null), pts_list)
crowns_list <- Filter(Negate(is.null), crowns_list)
chm_list <- Filter(Negate(is.null), chm_list)

if (length(pts_list) == 0 || length(crowns_list) == 0 || length(chm_list) == 0) stop("No trees found. Check inputs and parameters.")

pts_all_sv <- do.call(rbind, pts_list)
crowns_all_sv <- do.call(rbind, crowns_list)
chm_sv <- do.call(terra::merge, chm_list)

# ==============================================================================

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


las_sv <- spanner::colorize_las(las_sv, method = "attr",
                                attribute_name = "treeID",
                                palette = random.colors(500))
lidRviewer::view(las_sv)

# ==============================================================================
# ==============================================================================

message("Done. All algorithms ran successfully.")

library(ggplot2)
library(viridis)
library(sf)
library(terra)

# Convert CHM to data frame for ggplot
chm_df <- as.data.frame(chm, xy = TRUE, na.rm = TRUE)
names(chm_df)[3] <- "height"

ggplot() +
  geom_raster(
    data = chm_df,
    aes(x = x, y = y, fill = height)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Height (m)"
  ) +
  geom_sf(
    data = crowns_all,
    fill = NA,
    color = "white",
    linewidth = 0.4
  ) +
  # geom_sf(
  #   data = pts_all,
  #   aes(size = npoints),
  #   shape = 3,
  #   fill = NA,
  #   color = "white",
  #   stroke = 0.8
  # ) +
  scale_size_continuous(
    range = c(1, 6),
    name = "npoints"
  ) +
  coord_sf(expand = FALSE) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )
