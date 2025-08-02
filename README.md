# spanner <img src="https://github.com/bi0m3trics/spanner/blob/master/img/spanner_hex_logo.png" width="150" align="right"/>
![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
[![](https://www.r-pkg.org/badges/version/spanner)](https://cran.r-project.org/package=spanner)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4624277.svg)](https://doi.org/10.5281/zenodo.4624277)

Definition of spanner
<br/>1 (chiefly British): WRENCH
<br/>2: a wrench that has a hole, projection, or hook at one or both ends of the head for engaging with a corresponding device on the object that is to be turned
<br/>3: utilities to support landscape-, forest-, and tree-related data collection, manipulation, analysis, modelling, and visualization. 

<img src="https://github.com/bi0m3trics/spanner/blob/master/img/tshirt3.png" width="100%"  height="auto" align="center"/>

# Install `spanner`

Get the latest released version of spanner from github.

```r
remotes::install_github('bi0m3trics/spanner')
```

# Workflows
## Terrestial Lidar Tree Attributes and Segmentation

<img align="right" height="240" src="./img/output.gif">

The following is the full processing pipeline described in <a href="https://doi.org/10.3390/rs13122297">Donager et al. (2021)</a>, and provides an example from downloading an example dataset, preprocesing it using lidR's functionality, estimating tree locations and DBH by rasterizing individual point cloud values of relative neighborhood density (at 0.3 and 1 m radius) and verticality within a slice of the normalized point cloud around breast height to 
(1.37 m), to individual tree segmentation following ecological principles for “growing” trees based on input locations in a graph-theory approach. Relies heavily on work of <a href = "https://www.sciencedirect.com/science/article/pii/S0034425720304314">Roussel et al (2020)</a>, <a href="https://www.sciencedirect.com/science/article/abs/pii/S0924271615002373?via%3Dihub">Tao and others (2015)</a>, and <a href="https://www.sciencedirect.com/science/article/abs/pii/S0168169917301114?via%3Dihub">de Conto et al. (2017)</a>.<br/><br/>

Citation: Donager, Jonathon J., Andrew J. Sánchez Meador, and Ryan C. Blackburn 2021. Adjudicating Perspectives on Forest Structure: How Do Airborne, Terrestrial, and Mobile Lidar-Derived Estimates Compare? Remote Sensing 13, no. 12: 2297. https://doi.org/10.3390/rs13122297

## PatchMorph: Patch Delineation Algorithm

<img align="right" height="240" src="./img/pm_output.png">

The patchwoRk function implements a patch delineation algorithm [at present it only implements <a href="https://link.springer.com/article/10.1007/s10980-007-9104-8">'PatchMorph' (Girvetz & Greco 2007)</a>, which can delineate patches across a range of spatial scales based on three organism-specific thresholds - (1) land cover density threshold, (2) habitat gap maximum thickness (gap threshold), and (3) habitat patch minimum thickness (spur threshold)].

Citation: Girvetz EH, and Greco SE. 2007. How to define a patch: a spatial model for hierarchically delineating organism-specific habitat patches. Landscape Ecology 22: 1131-1142. http://dx.doi.org/10.1007/s10980-007-9104-8

## Key Features

### Optimized Tree Detection

The package provides a comprehensive tree detection workflow that combines:

1. **Optimized geometric feature computation** using `C_geometric_features_simple`
2. **Raster-based initial detection** with `get_raster_eigen_treelocs_optimized`
3. **Two-phase RANSAC cylinder fitting** for improved accuracy
4. **Enhanced filtering** for poor initial fits

### Main Functions

#### `comprehensive_tree_detection()`

The flagship function that provides a complete tree detection workflow:

```r
library(spanner)
library(lidR)

# Load your LAS data
las <- readLAS("your_file.laz")

# Run comprehensive tree detection
trees <- comprehensive_tree_detection(
  las = las,
  neigh_sizes = c(0.5, 0.25, 0.75),  # Radii for geometric analysis
  res = 0.05,                         # Grid resolution
  ransac_error_threshold = 0.15,      # Error threshold for cylinder fitting
  verbose = TRUE
)

# View results
print(trees)
plot_tree_detection(las, trees)
```

#### `get_raster_eigen_treelocs_optimized()`

Memory-efficient tree detection using direct raster operations:

```r
# Optimized tree detection
tree_candidates <- get_raster_eigen_treelocs_optimized(
  las = las,
  neigh_sizes = c(2, 1, 3),
  res = 0.02,
  dens_threshold = 0.1,
  vert_threshold = 0.5
)
```

#### `C_geometric_features_simple()`

Fast computation of geometric features with OpenMP parallelization:

```r
# Compute geometric features
features <- C_geometric_features_simple(
  las = las,
  radius = 0.5,
  max_neighbors = 50,
  ncpu = 4
)
```

## Workflow Overview

The comprehensive tree detection follows this three-phase approach:

### Phase 1: Initial Detection
- Downsample point cloud using systematic voxel grid
- Compute geometric features (verticality, density metrics)
- Create density and verticality rasters
- Identify candidate tree regions using raster operations
- Cluster candidate points to get initial tree locations

### Phase 2: RANSAC Cylinder Fitting
- Extract points around each tree candidate
- Fit cylinders using RANSAC algorithm
- Classify fits as "Good" or "Poor" based on error threshold
- Store tree parameters (location, radius, error)

### Phase 3: Enhanced Processing for Poor Fits
- Apply stricter geometric filtering (higher verticality, anisotropy thresholds)
- Refit cylinders with filtered point sets
- Reclassify improved fits as "Enhanced"

## Output

The `comprehensive_tree_detection()` function returns an `sf` object containing:

- **TreeID**: Unique identifier for each detected tree
- **X, Y, Z**: Tree location coordinates
- **Radius**: Estimated tree radius (DBH/2)
- **Error**: RANSAC fitting error
- **PointCount**: Number of points used for cylinder fitting
- **Quality**: Fit quality classification ("Good", "Enhanced", "Poor")
- **geometry**: Spatial coordinates as sf geometry

## Example Usage

### Basic Tree Detection

```r
library(spanner)
library(lidR)

# Load your data
las <- readLAS("forest_plot.laz")

# Run detection with default parameters
trees <- comprehensive_tree_detection(las, verbose = TRUE)

# Summary
summary(trees)
table(trees$Quality)
```

### Custom Parameters

```r
# Customize detection parameters
trees <- comprehensive_tree_detection(
  las = las,
  neigh_sizes = c(0.75, 0.3, 1.0),     # Larger neighborhoods
  res = 0.03,                          # Finer resolution
  grid_slice_min = 1.0,                # Lower slice
  grid_slice_max = 2.5,                # Higher slice
  dens_threshold = 0.12,               # Lower density threshold
  vert_threshold = 0.65,               # Higher verticality threshold
  ransac_error_threshold = 0.12,       # Stricter error threshold
  cylinder_search_radius = 0.5,        # Larger search radius
  verbose = TRUE
)
```

### Validation and Filtering

```r
# Validate detected trees
validated_trees <- validate_detected_trees(
  trees,
  min_radius = 0.08,      # Minimum 8cm radius
  max_radius = 0.8,       # Maximum 80cm radius
  max_error = 0.2,        # Maximum RANSAC error
  min_points = 15         # Minimum points per tree
)

# Plot results
plot_tree_detection(las, validated_trees, color_by = "Quality")
```

## Performance Notes

- The optimized geometric feature computation uses OpenMP for parallelization
- Memory usage is minimized through direct raster operations
- Processing time scales approximately linearly with point cloud size
- Recommended for point clouds up to ~50M points on standard hardware

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This package is licensed under GPL-3.

## Citation

If you use this package in your research, please cite:

```
Sanchez Meador, A., Donager, J., Blackburn, R., Cannon, J. (2025). spanner: Utilities for Forest and Tree Analysis. 
R package version 1.1.0. https://github.com/bi0m3trics/spanner
```
