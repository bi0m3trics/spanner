//  ===============================================================================
//
//  SSAO (Screen Space Ambient Occlusion) Implementation
//
//  Developers:
//
//  Andrew J. Sanchez Meador - andrew.sanchezmeador@nau.edu
//
//  COPYRIGHT: Andrew J. Sanchez Meador, 2026
//
//  This piece of software is open and free to use, redistribution and modifications
//  should be done in accordance to the GNU General Public License >= 3
//
//  Use this software as you wish, but no warranty is provided whatsoever.
//
//  Fast ambient occlusion using screen-space depth buffer technique.
//  Much faster than full 3D PCV as it works on projected 2D depth maps.
//
//  ===============================================================================

#ifndef SSAO_HPP
#define SSAO_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <Rcpp.h>
#include "SpatialIndex.h"
#include "Progress.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace lidR;

/**
 * Compute SSAO (Screen Space Ambient Occlusion) for a point cloud
 * 
 * @param las S4 LAS object with point cloud data
 * @param kernel_size Sampling kernel size in pixels
 * @param pixel_size Size of each pixel in world coordinates
 * @param num_samples Number of random samples per point
 * @param ncpu Number of CPU cores to use
 * @return Vector of SSAO values (one per point, 0=occluded, 1=visible)
 */
vector<double> computeSSAO(Rcpp::S4 las,
                           int kernel_size,
                           double pixel_size,
                           int num_samples,
                           int ncpu);

/**
 * Project 3D point cloud to 2D depth map
 * 
 * @param X X coordinates
 * @param Y Y coordinates  
 * @param Z Z coordinates
 * @param pixel_size Pixel size in world units
 * @param min_x Pre-computed minimum X
 * @param min_y Pre-computed minimum Y
 * @param depth_map Output depth map (min Z per pixel)
 * @param width Output map width
 * @param height Output map height
 * @param point_map Output map of point indices
 */
void projectToDepthMap(const Rcpp::NumericVector& X,
                      const Rcpp::NumericVector& Y,
                      const Rcpp::NumericVector& Z,
                      double pixel_size,
                      double min_x,
                      double min_y,
                      vector<vector<double>>& depth_map,
                      int& width,
                      int& height,
                      vector<vector<int>>& point_map);

/**
 * Compute SSAO for a single point using depth map
 * 
 * @param grid_x X position in depth map
 * @param grid_y Y position in depth map
 * @param depth Current point depth
 * @param depth_map The depth map
 * @param width Depth map width
 * @param height Depth map height
 * @param kernel_size Sampling radius
 * @param num_samples Number of samples
 * @return SSAO value
 */
double computeSSAOForPoint(int grid_x,
                          int grid_y,
                          double depth,
                          const vector<vector<double>>& depth_map,
                          int width,
                          int height,
                          int kernel_size,
                          int num_samples);

#endif // SSAO_HPP
