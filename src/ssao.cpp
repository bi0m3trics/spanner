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
//  ===============================================================================

#include "ssao.h"
#include <Rcpp.h>
#include <random>

const double PI = 3.14159265358979323846;

void projectToDepthMap(const Rcpp::NumericVector& X,
                      const Rcpp::NumericVector& Y,
                      const Rcpp::NumericVector& Z,
                      double pixel_size,
                      double min_x,
                      double min_y,
                      vector<vector<double>>& depth_map,
                      int& width,
                      int& height,
                      vector<vector<int>>& point_map) {
    
    int n_points = X.size();
    
    // Find max bounds
    double max_x = *std::max_element(X.begin(), X.end());
    double max_y = *std::max_element(Y.begin(), Y.end());
    
    // Calculate grid dimensions
    width = static_cast<int>((max_x - min_x) / pixel_size) + 1;
    height = static_cast<int>((max_y - min_y) / pixel_size) + 1;
    
    // Initialize depth map with max values (far plane)
    depth_map.resize(height);
    point_map.resize(height);
    for (int i = 0; i < height; ++i) {
        depth_map[i].resize(width, std::numeric_limits<double>::max());
        point_map[i].resize(width, -1);
    }
    
    // Project points to grid and keep closest (minimum Z)
    for (int i = 0; i < n_points; ++i) {
        int grid_x = static_cast<int>((X[i] - min_x) / pixel_size);
        int grid_y = static_cast<int>((Y[i] - min_y) / pixel_size);
        
        if (grid_x >= 0 && grid_x < width && grid_y >= 0 && grid_y < height) {
            if (Z[i] < depth_map[grid_y][grid_x]) {
                depth_map[grid_y][grid_x] = Z[i];
                point_map[grid_y][grid_x] = i;
            }
        }
    }
}

double computeSSAOForPoint(int grid_x,
                          int grid_y,
                          double depth,
                          const vector<vector<double>>& depth_map,
                          int width,
                          int height,
                          int kernel_size,
                          int num_samples) {
    
    double occlusion = 0.0;
    int valid_samples = 0;
    
    // Random number generator for sampling
    static std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    
    for (int s = 0; s < num_samples; ++s) {
        // Random offset within kernel
        int offset_x = static_cast<int>(dis(gen) * kernel_size);
        int offset_y = static_cast<int>(dis(gen) * kernel_size);
        
        int sample_x = grid_x + offset_x;
        int sample_y = grid_y + offset_y;
        
        // Check bounds
        if (sample_x < 0 || sample_x >= width || sample_y < 0 || sample_y >= height) {
            continue;
        }
        
        double sample_depth = depth_map[sample_y][sample_x];
        
        // Skip if no data at this location
        if (sample_depth == std::numeric_limits<double>::max()) {
            continue;
        }
        
        // Calculate depth difference
        double depth_diff = sample_depth - depth;
        
        // If sample is higher (lower Z in depth map = closer to camera = occluder)
        if (depth_diff < 0) {
            // Distance-based falloff
            double dist = sqrt(offset_x * offset_x + offset_y * offset_y);
            double falloff = 1.0 - (dist / kernel_size);
            
            // Angle-based occlusion (steeper angle = more occlusion)
            double angle_factor = std::min(1.0, abs(depth_diff) / (dist + 0.1));
            
            occlusion += falloff * angle_factor;
        }
        
        valid_samples++;
    }
    
    if (valid_samples > 0) {
        occlusion /= valid_samples;
    }
    
    // Return 1.0 - occlusion (so high values = less occluded = brighter)
    return 1.0 - std::min(1.0, occlusion);
}

vector<double> computeSSAO(Rcpp::S4 las,
                           int kernel_size,
                           double pixel_size,
                           int num_samples,
                           int ncpu) {
    
    // Extract point cloud data
    Rcpp::DataFrame data = Rcpp::as<Rcpp::DataFrame>(las.slot("data"));
    Rcpp::NumericVector X = data["X"];
    Rcpp::NumericVector Y = data["Y"];
    Rcpp::NumericVector Z = data["Z"];
    int n_points = X.size();
    
    vector<double> ssao_values(n_points, 1.0); // Default to fully lit
    
    Progress progress(n_points, "Computing SSAO (Screen Space Ambient Occlusion)");
    
    // Pre-compute bounds once
    double min_x = *std::min_element(X.begin(), X.end());
    double min_y = *std::min_element(Y.begin(), Y.end());
    
    // Project to depth map
    vector<vector<double>> depth_map;
    vector<vector<int>> point_map;
    int width, height;
    
    Rcpp::Rcout << "Projecting to depth map..." << std::endl;
    projectToDepthMap(X, Y, Z, pixel_size, min_x, min_y, depth_map, width, height, point_map);
    Rcpp::Rcout << "Depth map size: " << width << " x " << height << std::endl;
    
    // Compute SSAO for each point
    #ifdef _OPENMP
    omp_set_num_threads(ncpu);
    #endif
    
    bool abort = false;
    
    #pragma omp parallel for
    for (int i = 0; i < n_points; ++i) {
        if (abort) continue;
        if (progress.check_interrupt()) abort = true;
        progress.increment();
        
        int grid_x = static_cast<int>((X[i] - min_x) / pixel_size);
        int grid_y = static_cast<int>((Y[i] - min_y) / pixel_size);
        
        if (grid_x >= 0 && grid_x < width && grid_y >= 0 && grid_y < height) {
            ssao_values[i] = computeSSAOForPoint(
                grid_x, grid_y, Z[i],
                depth_map, width, height,
                kernel_size, num_samples
            );
        }
    }
    
    return ssao_values;
}
