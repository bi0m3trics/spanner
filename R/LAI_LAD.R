# ============================================================================
# spanner – LAI_LAD.R
#
# Leaf Area Index (LAI) and Plant Area Density (PAD) estimation from
# height-normalized lidar point clouds.
#
# Two approaches are implemented:
#
#   MacArthur-Horn  (compute_lai):
#     Applies the gap-fraction Beer-Lambert inversion from MacArthur &
#     Wilson (1966) as implemented in lidR::LAD(). Returns a SpatRaster of
#     LAI (and mean/max LAD) per grid cell. Works on any normalized single-
#     or multi-return lidar (ALS, UAS, TLS, MLS).
#
#   3-D voxel Beer-Lambert  (compute_pad_voxels):
#     Partitions the cloud into a 3-D voxel grid and, per vertical column,
#     counts directed / transmitted / intercepted pulses to estimate plant
#     area density (PAD, m²/m³) following Grau et al. (2017) and the
#     voxelmon framework (Tenny et al., 2025). Returns a data.table with
#     per-voxel metrics and an interpretive classification.
#
# Public functions:
#   compute_lai()                  – MacArthur-Horn LAI raster
#   compute_pad_voxels()           – Beer-Lambert 3-D voxel PAD
#   compute_transmittance_raster() – column transmittance (τ) raster from PAD
#
# Internal helpers:
#   .lai_cell()              – single-cell LAD profile + LAI summary
#   .compute_lai_las()       – pixel_metrics wrapper for a single LAS
#   .compute_lai_catalog()   – catalog_apply wrapper for a LAScatalog
# ============================================================================


# ----------------------------------------------------------------------------
#' Compute leaf area index from a normalized lidar point cloud
#'
#' `compute_lai` estimates leaf area index (LAI) and leaf area density (LAD)
#' per 2-D grid cell using the MacArthur-Horn gap-fraction Beer-Lambert
#' inversion as implemented in [lidR::LAD()]. The point cloud must be
#' height-normalized (`Z ≥ 0` above ground). Both a single `LAS` object and a
#' `LAScatalog` are accepted; catalog processing is parallelised over chunks.
#'
#' @details
#' The MacArthur-Horn method models the vertical gap-fraction profile as a
#' Beer-Lambert process:
#' \deqn{
#'   \mathrm{LAD}(z) = -\frac{\ln\!\bigl(\mathrm{gap\;fraction}(z)\bigr)}{k \cdot dz}
#' }
#' where \eqn{k} is the extinction coefficient (default 0.5 for a spherical
#' leaf angle distribution). LAI is the vertical integral:
#' \deqn{
#'   \mathrm{LAI} = \int \mathrm{LAD}(z)\,\mathrm{d}z
#'               \approx \sum_i \mathrm{LAD}(z_i) \cdot dz
#' }
#'
#' For `LAScatalog` input, set `opt_chunk_size()`, `opt_chunk_buffer()`, and
#' a `future::plan()` on the catalog before calling to control parallelism.
#' A buffer of at least `res_xy` is recommended to avoid edge artefacts.
#'
#' @param las A height-normalized `LAS` or `LAScatalog`. `Z` must represent
#'   height above ground (not raw elevation). The `ReturnNumber` attribute must
#'   be present (load with `select` containing `"r"`) so that the gap-fraction
#'   computation uses first returns only. Single-return scanners where every
#'   point has `ReturnNumber == 0` are handled automatically.
#' @param res_xy numeric. Grid cell size in the XY plane (m). Default `10`.
#' @param dz numeric. Vertical layer thickness for LAD profile integration
#'   (m). Default `1`.
#' @param k numeric. Beer-Lambert extinction coefficient. Default `0.5`
#'   (spherical leaf angle distribution).
#' @param min_height numeric. Minimum height above ground (m) included in
#'   the LAD profile. Points below this threshold are ignored. Default `0`.
#'
#' @return A [terra::SpatRaster] with three layers:
#'   \describe{
#'     \item{`LAI`}{Leaf area index (dimensionless). Vertical integral of
#'       LAD over all height bins.}
#'     \item{`LAD_mean`}{Mean LAD (m²/m³) across all positive-valued
#'       profile layers.}
#'     \item{`LAD_max`}{Maximum LAD (m²/m³) in the vertical profile.}
#'   }
#'   Grid cells with no returns return `NA`.
#'
#' @references
#' MacArthur, R.H. & Horn, J.W. (1969). Foliage profiles by vertical 
#' measurements. Ecology 50(5): 802–804. https://doi.org/10.2307/1933693
#'
#' Tenny JT, Sankey TT, Munson SM, Sánchez Meador AJ, Goetz SJ. (2025).
#' Canopy and surface fuels measurement using terrestrial lidar single-scan
#' approach in the Mogollon Highlands of Arizona. International Journal 
#' of Wildland Fire 34, WF24221. https://doi.org/10.1071/WF24221
#'
#' @seealso [compute_pad_voxels()], [lidR::LAD()]
#'
#' @examples
#' \donttest{
#' LASfile <- system.file("extdata", "MixedConifer.laz", package = "lidR")
#' las <- lidR::readLAS(LASfile)
#' las <- lidR::normalize_height(las, lidR::tin())
#' lai <- compute_lai(las, res_xy = 20, dz = 1, k = 0.5)
#' terra::plot(lai[["LAI"]], main = "MacArthur-Horn LAI")
#' terra::plot(lai[["LAD_mean"]], main = "Mean LAD (m2/m3)")
#' }
#'
#' @export
compute_lai <- function(las,
                        res_xy     = 10,
                        dz         = 1,
                        k          = 0.5,
                        min_height = 0) {

  if (!inherits(las, c("LAS", "LAScatalog")))
    stop("'las' must be a LAS or LAScatalog object.")
  if (!is.numeric(res_xy) || length(res_xy) != 1L || res_xy <= 0)
    stop("'res_xy' must be a single positive number.")
  if (!is.numeric(dz) || length(dz) != 1L || dz <= 0)
    stop("'dz' must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1L || k <= 0)
    stop("'k' must be a single positive number.")

  if (inherits(las, "LAScatalog"))
    return(.compute_lai_catalog(las, res_xy, dz, k, min_height))

  .compute_lai_las(las, res_xy, dz, k, min_height)
}


# ---- Internal helpers -------------------------------------------------------

#' @keywords internal
.lai_cell <- function(Z, ReturnNumber, dz, k, min_height) {
  # Use first returns only: multi-return pulses would make the intercepted
  # count exceed directed pulses in some bins, giving gap_frac < 0 and
  # log(negative) = NaN.  Fall back to all returns for single-return scanners
  # (ReturnNumber all 0 / all NA).
  fr <- !is.na(ReturnNumber) & ReturnNumber == 1L
  Z_use <- if (any(fr)) Z[fr] else Z[!is.na(Z)]
  Z_use <- Z_use[!is.na(Z_use)]

  if (length(Z_use) == 0L)
    return(list(LAI = NA_real_, LAD_mean = NA_real_, LAD_max = NA_real_))

  lad_obj  <- lidR::LAD(Z_use, dz = dz, k = k, z0 = min_height)
  lad_prof <- if (is.data.frame(lad_obj)) lad_obj[["lad"]] else lad_obj

  if (is.null(lad_prof) || length(lad_prof) == 0L)
    return(list(LAI = NA_real_, LAD_mean = NA_real_, LAD_max = NA_real_))

  # Guard against any residual NaN / Inf before summarising.
  lad_prof[!is.finite(lad_prof)] <- NA_real_

  LAI      <- sum(lad_prof, na.rm = TRUE) * dz
  finite   <- lad_prof[!is.na(lad_prof) & lad_prof > 0]
  LAD_mean <- if (length(finite) > 0L) mean(finite)   else 0
  LAD_max  <- if (length(finite) > 0L) max(finite)    else NA_real_

  list(LAI      = as.numeric(LAI),
       LAD_mean = as.numeric(LAD_mean),
       LAD_max  = as.numeric(LAD_max))
}


#' @keywords internal
.compute_lai_las <- function(las, res_xy, dz, k, min_height) {
  # lidR::pixel_metrics() deparsed the formula to a string and re-evaluates
  # it via eval(parse(text = cmd)), which drops all closure variables.
  # Fix: use bquote() to substitute literal parameter values into the formula
  # and reference the helper by its package-qualified name so it survives
  # the round-trip through deparse/parse.
  # ReturnNumber is included so .lai_cell can filter to first returns.
  func <- eval(bquote(
    ~ spanner:::.lai_cell(Z, ReturnNumber, .(dz), .(k), .(min_height))
  ))

  result <- lidR::pixel_metrics(las, func = func, res = res_xy)
  # lidR >= 4.1 already returns a SpatRaster; terra::rast(SpatRaster) copies
  # only the metadata, not the values, causing an empty raster.  Use result
  # directly when possible; only convert from stars/RasterBrick if needed.
  out <- if (inherits(result, "SpatRaster")) result else terra::rast(result)
  names(out) <- c("LAI", "LAD_mean", "LAD_max")
  out
}


#' @keywords internal
.compute_lai_catalog <- function(ctg, res_xy, dz, k, min_height) {
  # catalog_apply worker: each chunk is processed independently, then merged.
  # Bake parameter values into the formula with bquote for the same reason as
  # .compute_lai_las: lidR re-parses the formula string inside the worker,
  # so closure variables are not available.
  # ReturnNumber is included so .lai_cell can filter to first returns.
  func <- eval(bquote(
    ~ spanner:::.lai_cell(Z, ReturnNumber, .(dz), .(k), .(min_height))
  ))
  res_ <- res_xy

  .process_chunk <- function(cluster) {
    las <- lidR::readLAS(cluster)
    if (is.null(las) || lidR::is.empty(las)) return(NULL)

    result <- lidR::pixel_metrics(las, func = func, res = res_)
    if (is.null(result)) return(NULL)

    r <- if (inherits(result, "SpatRaster")) result else terra::rast(result)
    names(r) <- c("LAI", "LAD_mean", "LAD_max")

    # Crop to the chunk core to prevent seam duplication from the buffer
    bbox <- lidR::ext(cluster)
    if (!is.null(bbox)) terra::crop(r, bbox) else r
  }

  chunks <- lidR::catalog_apply(ctg, .process_chunk)
  chunks <- Filter(Negate(is.null), chunks)
  if (length(chunks) == 0L)
    stop("No valid chunks returned — check that the catalog is not empty.")
  do.call(terra::merge, chunks)
}


# ----------------------------------------------------------------------------
#' Estimate 3-D plant area density from a normalized lidar point cloud
#'
#' `compute_pad_voxels` partitions a height-normalized point cloud into a
#' regular 3-D voxel grid and, for each vertical column of voxels, traces the
#' Beer-Lambert attenuation of lidar pulses downward through the canopy. The
#' approach closely follows Grau et al. (2017) and the voxelmon framework
#' (Tenny et al., 2025).
#'
#' @details
#' **Pulse tracking (vertical column approximation)**
#'
#' First returns (or all returns for single-return scanners) are used as a
#' proxy for lidar pulses. For each voxel at grid position `(vx, vy, vz)`:
#'
#' \describe{
#'   \item{`directed`}{Total first returns in the same XY column at height
#'     ≥ `vz − vox_size/2`. These pulses were directed toward this voxel from
#'     above.}
#'   \item{`intercepted`}{First returns whose height falls inside this
#'     voxel's Z range `[vz − vox_size/2, vz + vox_size/2)`.}
#'   \item{`transmitted`}{Pulses that passed through this voxel without
#'     interception: `directed − intercepted`.}
#'   \item{`occluded`}{Ratio of pulses intercepted before reaching this
#'     voxel:
#'     \deqn{\mathrm{occluded} = 1 - \frac{\mathrm{transmitted} +
#'       \mathrm{intercepted}}{\mathrm{directed}}}
#'     Always `0` in the vertical column approximation (no horizontal
#'     shielding). Non-zero values arise with true 3-D ray tracing
#'     from a known scanner position (e.g., terrestrial single-scan mode).}
#' }
#'
#' **Plant Area Density**
#'
#' PAD is estimated by inverting the Beer-Lambert law for each voxel:
#' \deqn{
#'   \mathrm{PAD} = \frac{-\ln\!\left(1 - \dfrac{\mathrm{intercepted}}
#'     {\mathrm{intercepted} + \mathrm{transmitted}}\right)}
#'     {G \cdot L_{\mathrm{path}}}
#' }
#' where \eqn{G = 0.5} assumes a spherical leaf angle distribution and
#' \eqn{L_{\mathrm{path}} = f_{\mathrm{path}} \times \mathrm{vox\_size}}
#' with \eqn{f_{\mathrm{path}} = 0.843} following Grau et al. (2017), who
#' derived the expected mean path length through a cubic voxel for an
#' isotropic pulse distribution. For purely vertical (nadir-looking) ALS,
#' set `mean_path_factor = 1.0`.
#'
#' **Voxel classification** (after Tenny et al. 2025):
#' \describe{
#'   \item{`"foliage"`}{`PAD ∈ [min_pad, max_pad]` and `occluded ≤ max_occlusion`.}
#'   \item{`"wood"`}{`PAD > max_pad` and `occluded ≤ max_occlusion`.
#'     Solid returns (e.g., trunks, branches) produce very high PAD.}
#'   \item{`"empty"`}{`PAD < min_pad` and `occluded ≤ max_occlusion`.
#'     Open air with no intercepted pulses.}
#'   \item{`"occluded"`}{`occluded > max_occlusion` or `PAD` is non-finite.
#'     Voxel cannot be reliably characterised.}
#' }
#'
#' **LAI from PAD voxels**
#'
#' Column LAI can be derived by summing PAD × vox_size over foliage voxels:
#' ```r
#' pad_vox[VoxelClass == "foliage",
#'         .(LAI_pad = sum(PAD * vox_size, na.rm = TRUE)),
#'         by = .(X, Y)]
#' ```
#'
#' **Height profile**
#'
#' A vertical profile (similar to the voxelmon Profile output) can be
#' computed with:
#' ```r
#' pad_vox[, .(PAD_mean    = mean(PAD[VoxelClass == "foliage"], na.rm = TRUE),
#'             foliage_frac = mean(VoxelClass == "foliage"),
#'             empty_frac   = mean(VoxelClass == "empty"),
#'             occluded_frac = mean(VoxelClass == "occluded")),
#'         by = .(HAG_bin = floor(Z / vox_size) * vox_size)]
#' ```
#'
#' @param las A height-normalized `LAS` object. `Z` must represent height
#'   above ground. For catalogs, apply this function per chunk via
#'   [lidR::catalog_apply()].
#' @param vox_size numeric. Voxel edge length (m) in all three dimensions.
#'   Default `0.5`.
#' @param G numeric. Leaf projection coefficient. Default `0.5` (spherical
#'   leaf angle distribution).
#' @param mean_path_factor numeric. Multiplier applied to `vox_size` to
#'   obtain the mean pulse path length through a voxel. Default `0.843`
#'   (Grau et al. 2017; isotropic pulse directions). Use `1.0` for purely
#'   vertical (nadir) ALS pulses.
#' @param min_height numeric. Minimum height above ground (m) for included
#'   points. Default `0`.
#' @param max_occlusion numeric \[0, 1\]. Occlusion ratio above which a voxel
#'   is classified `"occluded"`. Default `0.8`.
#' @param min_pad numeric. Minimum PAD (m²/m³) for a voxel to be classified
#'   as `"foliage"`. Default `0.01`.
#' @param max_pad numeric. Maximum PAD (m²/m³) for a `"foliage"` voxel;
#'   voxels above this are classified as `"wood"`. Default `6.0`.
#' @param ncpu integer. Number of OpenMP threads for parallel voxelization and
#'   per-column cumulative sum. Default `1L`. Values greater than
#'   `parallel::detectCores()` are silently capped by the runtime. Set to
#'   `parallel::detectCores()` to use all available cores.
#'
#' @return A [data.table::data.table] with one row per occupied voxel and
#'   columns:
#'   \describe{
#'     \item{`X`, `Y`, `Z`}{numeric. Voxel centre coordinates (m). `Z` is
#'       height above ground for a normalized input.}
#'     \item{`HAG`}{numeric. Height above ground of the voxel centre (m).
#'       Equal to `Z` for a normalized point cloud.}
#'     \item{`directed`}{integer. Pulses directed toward this voxel from
#'       above (cumulative first returns in the column at ≥ this voxel's
#'       lower face).}
#'     \item{`transmitted`}{integer. Pulses that passed through this voxel
#'       without interception.}
#'     \item{`intercepted`}{integer. Pulses intercepted (returned) within
#'       this voxel.}
#'     \item{`occluded`}{numeric. Fraction of directed pulses intercepted
#'       before reaching this voxel. Always `0` under the vertical column
#'       approximation.}
#'     \item{`PAD`}{numeric. Plant area density (m²/m³).}
#'     \item{`VoxelClass`}{character. One of `"foliage"`, `"wood"`,
#'       `"empty"`, or `"occluded"`.}
#'   }
#'   Only voxels that contain at least one first-return point are included.
#'   Voxels with zero directed pulses (i.e., below the lowest return in a
#'   column) are not represented.
#'
#' @references
#' Grau, E., Durrieu, S., Fournier, R., Gastellu-Etchegorry, J.-P., &
#' Yin, T. (2017). Estimation of 3D vegetation density with Terrestrial
#' Laser Scanning data using voxels. A sensitivity analysis of influencing
#' parameters. *Remote Sensing of Environment*, 191, 373–388.
#' \doi{10.1016/j.rse.2017.01.032}
#'
#' Tenny JT, Sankey TT, Munson SM, Sánchez Meador AJ, Goetz SJ. (2025).
#' Canopy and surface fuels measurement using terrestrial lidar single-scan
#' approach in the Mogollon Highlands of Arizona. International Journal 
#' of Wildland Fire 34, WF24221. https://doi.org/10.1071/WF24221
#'
#' @seealso [compute_lai()], [lidR::LAD()]
#'
#' @examples
#' \donttest{
#' LASfile <- system.file("extdata", "TLS_Clip.laz", package = "spanner")
#' las <- lidR::readTLSLAS(LASfile, select = "xyzcr",
#'                         "-filter_with_voxel 0.01")
#' sf::st_crs(las) <- 26912
#' las <- lidR::normalize_height(las, lidR::tin())
#'
#' pad <- compute_pad_voxels(las, vox_size = 0.5)
#' print(pad)
#'
#' # Height profile of mean foliage PAD
#' profile <- pad[, .(PAD_mean = mean(PAD[VoxelClass == "foliage"],
#'                                    na.rm = TRUE)),
#'                by = .(HAG_bin = floor(Z / 0.5) * 0.5)]
#' print(profile[order(HAG_bin)])
#'
#' # LAI per XY column from foliage voxels
#' lai_pad <- pad[VoxelClass == "foliage",
#'                .(LAI_pad = sum(PAD * 0.5, na.rm = TRUE)),
#'                by = .(X, Y)]
#' print(lai_pad)
#' }
#'
#' @export
compute_pad_voxels <- function(las,
                               vox_size         = 0.5,
                               G                = 0.5,
                               mean_path_factor = 0.843,
                               min_height       = 0,
                               max_occlusion    = 0.8,
                               min_pad          = 0.01,
                               max_pad          = 6.0,
                               ncpu             = 1L) {

  if (!inherits(las, "LAS"))
    stop("'las' must be a LAS object. For catalogs, apply per-chunk via catalog_apply().")
  if (!is.numeric(vox_size) || length(vox_size) != 1L || vox_size <= 0)
    stop("'vox_size' must be a single positive number.")
  if (!is.numeric(G) || length(G) != 1L || G <= 0 || G > 1)
    stop("'G' must be in (0, 1]; 0.5 assumes a spherical leaf angle distribution.")
  if (!is.numeric(mean_path_factor) || length(mean_path_factor) != 1L ||
      mean_path_factor <= 0)
    stop("'mean_path_factor' must be a single positive number.")
  ncpu <- max(1L, as.integer(ncpu))

  dt <- las@data

  # -------------------------------------------------------------------
  # 1. Select first returns as pulse proxies.
  #    For single-return scanners (ReturnNumber absent), treat all returns
  #    as first returns (each point represents one complete pulse).
  # -------------------------------------------------------------------
  if ("ReturnNumber" %in% names(dt)) {
    pulse_mask <- !is.na(dt$ReturnNumber) & dt$ReturnNumber == 1L &
                  !is.na(dt$Z) & dt$Z >= min_height
  } else {
    pulse_mask <- !is.na(dt$Z) & dt$Z >= min_height
  }

  pulses <- dt[pulse_mask, .(X, Y, Z)]
  if (nrow(pulses) == 0L)
    stop("No valid first-return points after filtering. Check min_height and the LAS data.")

  # -------------------------------------------------------------------
  # 2–4. Voxelization + pulse counting via C++ (C_voxelize_pad).
  #
  #   Memory advantage: reads X/Y/Z by pointer with no intermediate copy;
  #   accumulates into a compact hash map (size = unique voxels << n_points).
  #   Speed advantage: OpenMP over both the voxelization pass (per-thread
  #   maps, serial merge) and the per-column cumulative sum.
  # -------------------------------------------------------------------
  cpp_out <- C_voxelize_pad(pulses$X, pulses$Y, pulses$Z,
                             vox_size, ncpu = ncpu)

  vox <- data.table::data.table(
    vx          = cpp_out$vx,
    vy          = cpp_out$vy,
    vz          = cpp_out$vz,
    intercepted = cpp_out$P_Intercepted,
    directed    = cpp_out$P_Directed,
    transmitted = cpp_out$P_Transmitted
  )

  # occluded: always 0 in the vertical column model (no lateral shadowing).
  vox[, occluded := 0]

  # -------------------------------------------------------------------
  # 3. Plant area density via Beer-Lambert inversion.
  #    PAD = -log(gap_frac) / (G * path_length)
  #    gap_frac = P_Transmitted / P_Directed
  #                = 1 - P_Intercepted / (P_Intercepted + P_Transmitted)
  #    path_length = mean_path_factor * vox_size
  # -------------------------------------------------------------------
  path_len <- mean_path_factor * vox_size
  gap_frac <- vox$transmitted / vox$directed

  # Clamp to (eps, 1]: log(0) = -Inf (fully intercepted voxels → very high
  # PAD → classified as "wood"); log > 0 impossible by construction.
  gap_frac <- pmin(pmax(gap_frac, .Machine$double.eps), 1.0)
  vox[, PAD := -log(gap_frac) / (G * path_len)]

  # -------------------------------------------------------------------
  # 4. Classify voxels.
  # -------------------------------------------------------------------
  vox[, VoxelClass := "empty"]
  vox[!is.finite(PAD) | occluded > max_occlusion,
      VoxelClass := "occluded"]
  vox[is.finite(PAD) & PAD >= min_pad & PAD <= max_pad &
        occluded <= max_occlusion,
      VoxelClass := "foliage"]
  vox[is.finite(PAD) & PAD > max_pad &
        occluded <= max_occlusion,
      VoxelClass := "wood"]

  # -------------------------------------------------------------------
  # 5. Final column order matching the documented output.
  # -------------------------------------------------------------------
  data.table::setnames(vox, c("vx", "vy", "vz"), c("X", "Y", "Z"))
  vox[, HAG := Z]   # equal to Z for a normalized point cloud

  data.table::setcolorder(vox, c("X", "Y", "Z", "HAG",
                                 "directed", "transmitted", "intercepted",
                                 "occluded", "PAD", "VoxelClass"))
  vox[]
}


# ----------------------------------------------------------------------------
#' Compute a column transmittance raster from PAD voxels
#'
#' `compute_transmittance_raster` converts the per-voxel PAD table produced by
#' [compute_pad_voxels()] into a 2-D raster of Beer-Lambert canopy
#' transmittance (τ ∈ \[0, 1\]) — the fraction of light from above that reaches
#' the forest floor for each (X, Y) column. The result can be used directly to
#' correct a vostokR `solar_potential` raster for canopy interception.
#'
#' @details
#' Per-column transmittance is derived by integrating PAD vertically through
#' all non-occluded voxels, with each voxel's PAD capped at `max_pad`:
#'
#' \deqn{
#'   \tau(x, y) = \exp\!\left(-G \sum_{i}
#'     \min(\mathrm{PAD}_{i},\, \mathrm{max\_pad}) \cdot \mathrm{vox\_size}\right)
#' }
#'
#' The cap is critical: wood voxels (class `"wood"`) have PAD values computed
#' from a gap fraction clamped at `.Machine$double.eps`, yielding PAD ≈ 171
#' m²/m³. Without a cap, even a single wood voxel per column drives τ to
#' essentially zero everywhere. Capping at `max_pad` (same value used in
#' [compute_pad_voxels()]) gives physically plausible column transmittance.
#'
#' Occluded voxels (class `"occluded"`) are excluded from the sum and treated
#' as transparent (conservative — avoids inflating occlusion where pulse
#' sampling is unreliable). `NA` PAD values are also skipped (`na.rm = TRUE`).
#'
#' **Correcting a vostokR solar potential raster:**
#' ```r
#' # solar_potential is in Wh m⁻² day⁻¹ (vostokR output)
#' # Resample tau to the same grid as solar_potential, then multiply:
#' tau_r      <- compute_transmittance_raster(pad, G = 0.5, vox_size = 0.5,
#'                                            max_pad = 6.0,
#'                                            crs = sf::st_crs(las)$wkt)
#' tau_rs      <- terra::resample(tau_r, solar_potential, method = "bilinear")
#' solar_floor <- solar_potential * tau_rs
#' ```
#'
#' @param pad     data.table returned by [compute_pad_voxels()].
#' @param G       numeric (length 1). Projection coefficient; must match the
#'   `G` argument used in [compute_pad_voxels()]. Default `0.5` (spherical
#'   leaf angle distribution).
#' @param vox_size numeric (length 1, > 0). Voxel edge length in metres; must
#'   match the `vox_size` argument used in [compute_pad_voxels()].
#' @param max_pad numeric (length 1, > 0). Cap applied to each voxel's PAD
#'   before summing. Should match the `max_pad` argument used in
#'   [compute_pad_voxels()] (default `6.0`). This prevents numerically extreme
#'   wood-voxel PAD values from collapsing transmittance to zero.
#' @param crs character or CRS object (optional). Coordinate reference system
#'   for the output raster, e.g. `sf::st_crs(las)$wkt`. If `NULL` the raster
#'   is returned without a CRS.
#'
#' @return A single-layer [terra::SpatRaster] named `"tau"` with resolution
#'   equal to `vox_size`. Cell values are transmittance fractions in \[0, 1\]:
#'   \describe{
#'     \item{1.0}{fully open — no canopy above that column}
#'     \item{0.0}{fully opaque — all light intercepted}
#'   }
#'
#' @seealso [compute_pad_voxels()], [compute_lai()]
#'
#' @examples
#' \dontrun{
#' las <- lidR::readLAS(system.file("extdata", "ALS_Clip.laz", package = "spanner"))
#' pad <- compute_pad_voxels(las, vox_size = 0.5, max_pad = 6.0)
#' tau_r <- compute_transmittance_raster(pad, G = 0.5, vox_size = 0.5,
#'                                       max_pad = 6.0,
#'                                       crs = sf::st_crs(las)$wkt)
#' terra::plot(tau_r, main = "Column transmittance (tau)", range = c(0, 1))
#'
#' # Correct a vostokR solar_potential raster (must be computed separately):
#' # tau_rs      <- terra::resample(tau_r, solar_potential, method = "bilinear")
#' # solar_floor <- solar_potential * tau_rs
#' }
#'
#' @export
compute_transmittance_raster <- function(pad,
                                         G        = 0.5,
                                         vox_size,
                                         max_pad  = 6.0,
                                         crs      = NULL) {

  if (!data.table::is.data.table(pad))
    stop("'pad' must be a data.table returned by compute_pad_voxels().")
  required_cols <- c("X", "Y", "PAD", "VoxelClass")
  missing_cols  <- setdiff(required_cols, names(pad))
  if (length(missing_cols) > 0L)
    stop("'pad' is missing columns: ", paste(missing_cols, collapse = ", "),
         ". Pass the data.table returned by compute_pad_voxels().")
  if (!is.numeric(G) || length(G) != 1L || G <= 0 || G > 1)
    stop("'G' must be a single number in (0, 1]; use the same value as compute_pad_voxels().")
  if (missing(vox_size) || !is.numeric(vox_size) || length(vox_size) != 1L ||
      vox_size <= 0)
    stop("'vox_size' must be a single positive number matching the compute_pad_voxels() call.")
  if (!is.numeric(max_pad) || length(max_pad) != 1L || max_pad <= 0)
    stop("'max_pad' must be a single positive number; use the same value as compute_pad_voxels().")

  # Sum PAD (capped at max_pad) over all non-occluded voxels per (X, Y) column.
  # The cap is essential: wood voxels have uncapped PAD ≈ 171 m²/m³ (derived
  # from gap_frac clamped at .Machine$double.eps). Without the cap, even one
  # wood voxel per column collapses τ = exp(-G × 85) → 0 everywhere.
  # Occluded voxels are excluded (treated as transparent — conservative).
  col_lai <- pad[VoxelClass != "occluded",
                 .(col_LAI = sum(pmin(PAD, max_pad) * vox_size, na.rm = TRUE)),
                 by = .(X, Y)]

  # Beer-Lambert column transmittance: τ = exp(-G × Σ min(PADᵢ, max_pad) × vox_size)
  col_lai[, tau := exp(-G * col_LAI)]

  # Build SpatRaster from the regular (X, Y) grid.
  # terra::rast(df, type="xyz") infers resolution from point spacing, which
  # equals vox_size for the regular grid produced by compute_pad_voxels().
  crs_str <- if (!is.null(crs)) as.character(crs) else ""
  r <- terra::rast(col_lai[, .(x = X, y = Y, tau)],
                   type = "xyz",
                   crs  = crs_str)
  names(r) <- "tau"
  r
}
