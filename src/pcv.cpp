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
//  ===============================================================================

#include "pcv.h"

const double PI = 3.14159265358979323846;
const double HALF_PI = PI / 2.0;

double computeMaxAngleInDirection(const PointXYZ& point,
                                  const vector<PointXYZ>& neighbors,
                                  double dir_x,
                                  double dir_y) {
    double max_angle = -std::numeric_limits<double>::infinity();
    
    const double cos_threshold = 0.7; // ~45 degrees
    
    for (const auto& neighbor : neighbors) {
        // Compute horizontal distance and direction
        double dx = neighbor.x - point.x;
        double dy = neighbor.y - point.y;
        double dist_2d = sqrt(dx * dx + dy * dy);
        
        if (dist_2d < 1e-10) continue;
        
        // Normalize direction vector
        double ndx = dx / dist_2d;
        double ndy = dy / dist_2d;
        
        // Check if neighbor is in this direction (dot product)
        double dot = ndx * dir_x + ndy * dir_y;
        
        if (dot > cos_threshold) {
            // Compute elevation angle
            double dz = neighbor.z - point.z;
            double angle = atan2(dz, dist_2d);
            max_angle = std::max(max_angle, angle);
        }
    }
    
    return max_angle;
}

double computePCVForPoint(const PointXYZ& point,
                          SpatialIndex& tree,
                          double radius,
                          int num_directions) {
    
    // Find neighbors using spatial index
    vector<PointXYZ> neighbors;
    Sphere sphere(point.x, point.y, point.z, radius);
    tree.lookup(sphere, neighbors);
    
    if (neighbors.empty()) {
        return HALF_PI; // Maximum visibility (no occlusion)
    }
    
    // Compute sky visibility for each direction
    double total_angle = 0.0;
    
    for (int i = 0; i < num_directions; ++i) {
        double angle = 2.0 * PI * i / num_directions;
        double dir_x = cos(angle);
        double dir_y = sin(angle);
        
        double max_angle = computeMaxAngleInDirection(point, neighbors, dir_x, dir_y);
        
        // Sky portion in this direction
        if (std::isfinite(max_angle)) {
            total_angle += (HALF_PI - max_angle);
        } else {
            total_angle += HALF_PI;
        }
    }
    
    return total_angle / num_directions;
}

vector<double> computePCV(Rcpp::S4 las,
                          double radius,
                          int num_directions,
                          int ncpu) {
    
    // Extract point cloud data
    Rcpp::DataFrame data = Rcpp::as<Rcpp::DataFrame>(las.slot("data"));
    Rcpp::NumericVector X = data["X"];
    Rcpp::NumericVector Y = data["Y"];
    Rcpp::NumericVector Z = data["Z"];
    int n_points = X.size();
    
    vector<double> pcv_values(n_points);
    Progress pb(n_points, "Computing PCV ambient occlusion");
    
    // Build spatial index (automatically selects GridPartition or SparsePartition)
    SpatialIndex tree(las);
    
    #ifdef _OPENMP
    omp_set_num_threads(ncpu);
    #endif
    
    bool abort = false;
    
    #pragma omp parallel for
    for (int i = 0; i < n_points; ++i) {
        if (abort) continue;
        if (pb.check_interrupt()) abort = true;
        pb.increment();
        
        PointXYZ point(X[i], Y[i], Z[i], i);
        pcv_values[i] = computePCVForPoint(point, tree, radius, num_directions);
    }
    
    return pcv_values;
}
