//  ===============================================================================
//
//  PCV (Portion de Ciel Visible) Implementation
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
//  Based on:
//  Duguet, Florent & Girardeau-Montaut, Daniel. (2004). 
//  Rendu en Portion de Ciel Visible de Gros Nuages de Points 3D.
//
//  This implementation computes ambient occlusion for point clouds by calculating
//  the visible portion of the sky from each point.
//
//  ===============================================================================

#ifndef PCV_HPP
#define PCV_HPP

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
 * Compute PCV (sky visibility) for a point cloud
 * 
 * @param las S4 LAS object with point cloud data
 * @param radius Search radius for neighborhood
 * @param num_directions Number of directional rays to cast
 * @param ncpu Number of CPU cores to use
 * @return Vector of PCV values (one per point)
 */
vector<double> computePCV(Rcpp::S4 las,
                          double radius, 
                          int num_directions,
                          int ncpu);

/**
 * Compute PCV for a single point using spatial index
 * 
 * @param point Current point coordinates and index
 * @param tree Spatial index for neighbor lookup
 * @param radius Search radius
 * @param num_directions Number of rays
 * @return PCV value for the point
 */
double computePCVForPoint(const PointXYZ& point,
                          SpatialIndex& tree,
                          double radius,
                          int num_directions);

/**
 * Compute maximum elevation angle in a given direction
 * 
 * @param point Query point
 * @param neighbors Vector of neighbor points
 * @param dir_x Direction X component
 * @param dir_y Direction Y component
 * @return Maximum elevation angle in radians
 */
double computeMaxAngleInDirection(const vector<double>& point,
                                  const vector<vector<double>>& neighbors,
                                  double dir_x,
                                  double dir_y);

#endif // PCV_HPP
