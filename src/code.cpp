/*
 * ======= CALC EIGEN METRICS IN SPHERE FUNCTION =========
 */
// [[Rcpp::depends(lidR)]]

#include <RcppArmadillo.h>
#include <limits>
#include <cmath>
#include "SpatialIndex.h"
#include "Progress.h"

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace lidR;
using namespace Rcpp;

// [[Rcpp::export]]
List C_eigen_in_sphere(S4 las, double radius, int ncpu)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  int n = X.size();
  int n_metrics = 27;  // Increased to include all CloudCompare metrics

  //https://ethz.ch/content/dam/ethz/special-interest/baug/igp/photogrammetry-remote-sensing-dam/documents/pdf/timo-jan-cvpr2016.pdf
  NumericVector eigenlar_sph(n);
  NumericVector eigenmed_sph(n);
  NumericVector eigensmall_sph(n);
  NumericVector eigensum_sph(n);
  NumericVector curvature_sph(n);
  NumericVector omnivariance_sph(n);
  NumericVector anisotropy_sph(n);
  NumericVector eigentropy_sph(n);
  NumericVector linearity_sph(n);
  NumericVector verticality_sph(n);
  NumericVector planarity_sph(n);
  NumericVector sphericity_sph(n);
  NumericVector nx_sph(n);
  NumericVector ny_sph(n);
  NumericVector nz_sph(n);
  
  // Additional CloudCompare metrics
  NumericVector surface_variation_sph(n);
  NumericVector change_curvature_sph(n);
  NumericVector surface_density_sph(n);
  NumericVector volume_density_sph(n);
  NumericVector moment_order1_sph(n);
  NumericVector normal_change_rate_sph(n);
  
  // Missing CloudCompare features
  NumericVector roughness_sph(n);
  NumericVector mean_curvature_sph(n);
  NumericVector gaussian_curvature_sph(n);
  NumericVector pca1_sph(n);
  NumericVector pca2_sph(n);
  NumericVector num_neighbors_sph(n);

  List out(n_metrics);

  // One copy for the index. The index automatically detects the TLS tag in the LAS object
  lidR::GridPartition tree(las);

  // Initiate the progress bar
  Progress pb(n, "Calculating eigen metrics in sphere: ");

  bool abort = false;

#pragma omp parallel for num_threads(ncpu)
  for(unsigned int i = 0 ; i < n ; i++)
  {

    if (abort) continue;
    if (pb.check_interrupt()) abort = true;
    pb.increment();

    std::vector<PointXYZ> sphpts;                   // creation of an STL container of points for the sphere neighborhood object
    Sphere sphere(X[i], Y[i], Z[i], radius);        // creation of a sphere object
    tree.lookup(sphere, sphpts);                    // lookup the points in the sphere neighborhood

    int k = sphpts.size();  // Determine the neighborhood size of the sphere

    // Initialize default values for edge cases
    double eigen_largest = 0.0, eigen_medium = 0.0, eigen_smallest = 0.0;
    double eigensum = 0.0, curvature = 0.0, omnivariance = 0.0;
    double anisotropy = 0.0, eigentropy = 0.0, linearity = 0.0;
    double verticality = 0.0, planarity = 0.0, sphericity = 0.0;
    double nx = 0.0, ny = 0.0, nz = 0.0;
    
    // Additional CloudCompare metrics initialization
    double surface_variation = 0.0, change_curvature = 0.0;
    double surface_density = 0.0, volume_density = 0.0;
    double moment_order1 = 0.0, normal_change_rate = 0.0;
    double roughness = 0.0, mean_curvature = 0.0, gaussian_curvature = 0.0;
    double pca1 = 0.0, pca2 = 0.0;
    int num_neighbors = k;
    
    // Compute centroid for all points (needed for some metrics)
    double cx = 0.0, cy = 0.0, cz = 0.0;
    if (k > 0) {
      for (unsigned int j = 0; j < sphpts.size(); j++) {
        cx += sphpts[j].x;
        cy += sphpts[j].y;
        cz += sphpts[j].z;
      }
      cx /= k; cy /= k; cz /= k;
    }

    // Only perform eigenvalue decomposition if we have enough points
    if (k >= 3) {
      // Compute covariance matrix manually (like CloudCompare)
      arma::mat cov_matrix(3, 3, arma::fill::zeros);
      for (unsigned int j = 0; j < sphpts.size(); j++) {
        double dx = sphpts[j].x - cx;
        double dy = sphpts[j].y - cy;
        double dz = sphpts[j].z - cz;
        
        cov_matrix(0,0) += dx * dx;
        cov_matrix(0,1) += dx * dy;
        cov_matrix(0,2) += dx * dz;
        cov_matrix(1,0) += dy * dx;
        cov_matrix(1,1) += dy * dy;
        cov_matrix(1,2) += dy * dz;
        cov_matrix(2,0) += dz * dx;
        cov_matrix(2,1) += dz * dy;
        cov_matrix(2,2) += dz * dz;
      }
      
      // Normalize covariance matrix by k (not k-1) to match CloudCompare
      cov_matrix /= k;
      
      // Compute eigenvalues and eigenvectors
      arma::vec eigenvalues;
      arma::mat eigenvectors;
      arma::eig_sym(eigenvalues, eigenvectors, cov_matrix);
      
      // Sort eigenvalues in descending order (CloudCompare convention)
      arma::uvec sorted_indices = arma::sort_index(eigenvalues, "descend");
      
#pragma omp critical
{
  // Calculate decomposed eigen values

  eigen_largest  = eigenvalues[sorted_indices[0]];
  eigen_medium   = eigenvalues[sorted_indices[1]];
  eigen_smallest = eigenvalues[sorted_indices[2]];
  eigensum       = eigen_largest + eigen_medium + eigen_smallest;
  curvature      = eigen_smallest/eigensum;
  
  // CloudCompare-style omnivariance calculation using raw eigenvalues
  if (eigen_largest <= 0.0 || eigen_medium <= 0.0 || eigen_smallest <= 0.0) {
    omnivariance = 0.0;
  } else {
    // Use raw eigenvalues (not normalized) like CloudCompare
    omnivariance = pow((eigen_largest * eigen_medium * eigen_smallest), (1.0/3.0));
  }
  
  anisotropy     = (eigen_largest - eigen_smallest) / eigen_largest;
  
  // CloudCompare: eigentropy uses RAW eigenvalues (no normalization)
  // From CCCoreLib Neighbourhood.cpp line 796
  if (eigen_largest > 0.0 && eigen_medium > 0.0 && eigen_smallest > 0.0) {
    eigentropy = -((eigen_largest * log(eigen_largest)) + (eigen_medium * log(eigen_medium)) + (eigen_smallest * log(eigen_smallest)));
  } else {
    eigentropy = 0.0;
  }
  
  linearity      = (eigen_largest - eigen_medium) / eigen_largest;
  
  // Get the eigenvector corresponding to the smallest eigenvalue (normal vector)
  arma::vec normal = eigenvectors.col(sorted_indices[2]);
  
  verticality    = 1 - std::abs(normal[2]);
  planarity      = (eigen_medium - eigen_smallest) / eigen_largest;
  sphericity     = eigen_smallest / eigen_largest;
  nx = normal[0];
  ny = normal[1];
  nz = normal[2];
  
  // Additional CloudCompare metrics calculations
  
  // Surface variation (change of curvature) - λ3 / (λ1 + λ2 + λ3)
  surface_variation = curvature;  // This is the same as curvature in CloudCompare
  
  // Change of curvature (alternative name for surface variation)
  change_curvature = surface_variation;
  
  // Surface density - CloudCompare uses DENSITY_2D: divides by circle area (πR²), not sphere surface area
  // From GeometricalAnalysisTools.cpp line 143: dimensionalCoef = M_PI * pow(kernelRadius, 2.0)
  double circle_area = M_PI * radius * radius;
  surface_density = k / circle_area;
  
  // Volume density - number of neighbors per unit volume  
  double sphere_volume = (4.0/3.0) * M_PI * radius * radius * radius;
  volume_density = k / sphere_volume;
  
  // 1st order moment - CloudCompare formula from Neighbourhood.cpp line 726-757
  // Uses projection onto 2nd eigenvector from QUERY POINT (not centroid)
  arma::vec e2 = eigenvectors.col(sorted_indices[1]);
  double m1 = 0.0;
  double m2 = 0.0;
  for (unsigned int j = 0; j < sphpts.size(); ++j) {
    // Compute vector from query point to neighbor (not from centroid!)
    double px = sphpts[j].x - X[i];
    double py = sphpts[j].y - Y[i];
    double pz = sphpts[j].z - Z[i];
    double dotProd = px * e2[0] + py * e2[1] + pz * e2[2];
    m1 += dotProd;
    m2 += dotProd * dotProd;
  }
  moment_order1 = (m2 < 1e-10) ? 0.0 : (m1 * m1) / m2;
  
  // Normal change rate - CloudCompare uses surface variation (curvature)
  // From CCCoreLib: normal_change_rate = l3 / (l1 + l2 + l3)
  normal_change_rate = curvature;  // curvature is already l3/sum
  
  // Roughness - distance from query point to fitted plane
  // Plane passes through centroid with normal as the smallest eigenvector
  double dx_query = X[i] - cx;
  double dy_query = Y[i] - cy;
  double dz_query = Z[i] - cz;
  roughness = fabs(dx_query * normal[0] + dy_query * normal[1] + dz_query * normal[2]);
  
  // Mean and Gaussian curvature - CloudCompare quadric surface fitting
  // From CCCoreLib/Neighbourhood.cpp lines 298-468 and 886-940
  // Fits a 2.5D quadric: Z = a + b*X + c*Y + d*X^2 + e*X*Y + f*Y^2
  if (k >= 5) {  // Need at least 5 points for quadric fitting
    // Step 1: Create local coordinate system where LS plane normal becomes Z axis
    // We already have the normal vector and eigenvectors
    arma::vec N = normal;  // Smallest eigenvector (normal to plane)
    
    // Create orthonormal basis for local coordinate system
    // Y vector: orthogonal to N, stable orientation
    arma::vec Y_local(3);
    Y_local(0) = -N(1);
    Y_local(1) = N(0);
    Y_local(2) = 0.0;
    double Y_norm = arma::norm(Y_local);
    
    if (Y_norm < 0.1) {  // N is nearly vertical, use different Y
      Y_local(0) = 0.0;
      Y_local(1) = N(2);
      Y_local(2) = -N(1);
      Y_norm = arma::norm(Y_local);
    }
    
    if (Y_norm > 1e-10) {
      Y_local = Y_local / Y_norm;  // Normalize
      
      // X = Y × N (cross product)
      arma::vec X_local = arma::cross(Y_local, N);
      X_local = X_local / arma::norm(X_local);  // Normalize
      
      // Step 2: Build the least squares system for quadric fitting
      // We solve: A*h = b where h = [a, b, c, d, e, f]^T
      // Each point contributes: [1, x, y, x^2, x*y, y^2] * h = z
      arma::mat A(k, 6, arma::fill::zeros);
      arma::vec b(k);
      
      double max_dim_sq = 0.0;  // Track max squared dimension for numerical stability
      
      for (unsigned int j = 0; j < k; ++j) {
        // Transform point to local coordinate system
        double dx = sphpts[j].x - cx;
        double dy = sphpts[j].y - cy;
        double dz = sphpts[j].z - cz;
        
        double lx = dx * X_local(0) + dy * X_local(1) + dz * X_local(2);
        double ly = dx * Y_local(0) + dy * Y_local(1) + dz * Y_local(2);
        double lz = dx * N(0) + dy * N(1) + dz * N(2);
        
        // Build matrix row: [1, x, y, x^2, x*y, y^2]
        A(j, 0) = 1.0;
        A(j, 1) = lx;
        A(j, 2) = ly;
        A(j, 3) = lx * lx;
        A(j, 4) = lx * ly;
        A(j, 5) = ly * ly;
        
        b(j) = lz;
        
        // Track maximum dimension
        max_dim_sq = std::max(max_dim_sq, lx * lx);
        max_dim_sq = std::max(max_dim_sq, ly * ly);
        max_dim_sq = std::max(max_dim_sq, lz * lz);
      }
      
      // Step 3: Solve least squares system using Armadillo
      // Check if system is well-conditioned
      if (max_dim_sq > 1e-10) {
        arma::vec h;
        bool solve_ok = arma::solve(h, A, b, arma::solve_opts::fast);
        
        if (solve_ok && h.n_elem == 6) {
          // Quadric coefficients: Z = a + b*X + c*Y + d*X^2 + e*X*Y + f*Y^2
          double a = h(0);
          double b_coef = h(1);
          double c_coef = h(2);
          double d_coef = h(3);
          double e_coef = h(4);
          double f_coef = h(5);
          
          // Step 4: Transform query point to local coordinates
          double dx_q = X[i] - cx;
          double dy_q = Y[i] - cy;
          double dz_q = Z[i] - cz;
          
          double qx = dx_q * X_local(0) + dy_q * X_local(1) + dz_q * X_local(2);
          double qy = dx_q * Y_local(0) + dy_q * Y_local(1) + dz_q * Y_local(2);
          
          // Step 5: Compute partial derivatives at query point
          // See "CURVATURE OF CURVES AND SURFACES – A PARABOLIC APPROACH" by ZVI HAR'EL
          double fx = b_coef + 2.0 * d_coef * qx + e_coef * qy;
          double fy = c_coef + e_coef * qx + 2.0 * f_coef * qy;
          double fxx = 2.0 * d_coef;
          double fyy = 2.0 * f_coef;
          double fxy = e_coef;
          
          double fx2 = fx * fx;
          double fy2 = fy * fy;
          double q = 1.0 + fx2 + fy2;
          
          // Step 6: Compute principal curvatures using differential geometry
          // For a surface z = f(x,y), use the Weingarten map (shape operator)
          
          // Direct computation using the formulas for principal curvatures
          // from differential geometry without matrix eigenvalue decomposition
          
          // The principal curvatures k1, k2 satisfy:
          // Mean curvature: H = (k1 + k2)/2
          // Gaussian curvature: K = k1 * k2
          
          // For a Monge patch z = f(x,y), the standard formulas are:
          // H = ((1+fy²)fxx - 2fx·fy·fxy + (1+fx²)fyy) / (2(1+fx²+fy²)^(3/2))
          // K = (fxx·fyy - fxy²) / (1+fx²+fy²)²
          
          double sqrt_q = std::sqrt(q);
          double q_cubed = q * sqrt_q;  // q^(3/2)
          
          // Mean curvature (signed)
          double H_numerator = (1.0 + fy2) * fxx - 2.0 * fx * fy * fxy + (1.0 + fx2) * fyy;
          double H = H_numerator / (2.0 * q_cubed);
          
          // Gaussian curvature
          double K = (fxx * fyy - fxy * fxy) / (q * q);
          
          // Take absolute values as CloudCompare does
          mean_curvature = std::abs(H);
          gaussian_curvature = std::abs(K);
        }
      }
    }
  }
  
  // PCA1 and PCA2 - Full eigenvector transformation approach
  // Project neighborhood onto principal components and measure spread
  if (eigensum > std::numeric_limits<double>::epsilon()) {
    // Get the first two principal component vectors
    arma::vec pc1 = eigenvectors.col(sorted_indices[0]);  // First principal component
    arma::vec pc2 = eigenvectors.col(sorted_indices[1]);  // Second principal component
    
    // Project all points onto PC1 and PC2 and compute variance of projections
    double var_pc1 = 0.0, var_pc2 = 0.0;
    for (unsigned int j = 0; j < k; ++j) {
      double dx = sphpts[j].x - cx;
      double dy = sphpts[j].y - cy;
      double dz = sphpts[j].z - cz;
      
      // Project onto PC1 and PC2
      double proj1 = dx * pc1(0) + dy * pc1(1) + dz * pc1(2);
      double proj2 = dx * pc2(0) + dy * pc2(1) + dz * pc2(2);
      
      var_pc1 += proj1 * proj1;
      var_pc2 += proj2 * proj2;
    }
    
    // Normalize by number of points (variance of projections)
    var_pc1 /= k;
    var_pc2 /= k;
    
    // PCA1 and PCA2 as normalized projection variances
    pca1 = var_pc1 / eigensum;
    pca2 = var_pc2 / eigensum;
  } else {
    pca1 = 0.0;
    pca2 = 0.0;
  }
}
    } // end if k >= 3
    
#pragma omp critical
{

  eigenlar_sph[i]   = eigen_largest;
  eigenmed_sph[i]   = eigen_medium;
  eigensmall_sph[i] = eigen_smallest;
  eigensum_sph[i]   = eigensum;
  curvature_sph[i]  = curvature;
  omnivariance_sph[i] = omnivariance;
  anisotropy_sph[i] = anisotropy;
  eigentropy_sph[i] = eigentropy;
  linearity_sph[i]  = linearity;
  verticality_sph[i] = verticality;
  planarity_sph[i]  = planarity;
  sphericity_sph[i] = sphericity;
  nx_sph[i] = nx;
  ny_sph[i] = ny;
  nz_sph[i] = nz;
  
  // Additional CloudCompare metrics assignments
  surface_variation_sph[i] = surface_variation;
  change_curvature_sph[i] = change_curvature;
  surface_density_sph[i] = surface_density;
  volume_density_sph[i] = volume_density;
  moment_order1_sph[i] = moment_order1;
  normal_change_rate_sph[i] = normal_change_rate;
  
  // Missing CloudCompare features assignments
  roughness_sph[i] = roughness;
  mean_curvature_sph[i] = mean_curvature;
  gaussian_curvature_sph[i] = gaussian_curvature;
  pca1_sph[i] = pca1;
  pca2_sph[i] = pca2;
  num_neighbors_sph[i] = num_neighbors;

}
  }
  if (abort) throw Rcpp::internal::InterruptedException();

  // Build and return the list of numericvectors
  out[0] = eigenlar_sph;
  out[1] = eigenmed_sph;
  out[2] = eigensmall_sph;
  out[3] = eigensum_sph;
  out[4] = curvature_sph;
  out[5] = omnivariance_sph;
  out[6] = anisotropy_sph;
  out[7] = eigentropy_sph;
  out[8] = linearity_sph;
  out[9] = verticality_sph;
  out[10] = planarity_sph;
  out[11] = sphericity_sph;
  out[12] = nx_sph;
  out[13] = ny_sph;
  out[14] = nz_sph;
  out[15] = surface_variation_sph;
  out[16] = change_curvature_sph;
  out[17] = surface_density_sph;
  out[18] = volume_density_sph;
  out[19] = moment_order1_sph;
  out[20] = normal_change_rate_sph;
  out[21] = roughness_sph;
  out[22] = mean_curvature_sph;
  out[23] = gaussian_curvature_sph;
  out[24] = pca1_sph;
  out[25] = pca2_sph;
  out[26] = num_neighbors_sph;

  return out;
}


/*
 * ======= LIDRPLUGINS' COUNT IN DISC FUNCTION =========
 */

// [[Rcpp::depends(lidR)]]

#include <RcppArmadillo.h>
#include "SpatialIndex.h"
#include "Progress.h"

using namespace lidR;
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector C_count_in_disc(NumericVector X, NumericVector Y, NumericVector x, NumericVector y, double radius, int ncpu)
{
  unsigned int n = x.length();
  IntegerVector output(n);

  lidR::GridPartition tree(X,Y);

  Progress pb(n, "Counting points in disc: ");
  bool abort= false;

#pragma omp parallel for num_threads(ncpu)
  for(unsigned int i = 0 ; i < n ; i++)
  {
    if (abort) continue;
    if (pb.check_interrupt()) abort = true;
    pb.increment();

    lidR::Circle disc(x[i], y[i], radius);
    std::vector<lidR::PointXYZ> pts;
    tree.lookup(disc, pts);

#pragma omp critical
{
  output[i] = pts.size();
}
  }
  if (abort) throw Rcpp::internal::InterruptedException();
  return output;
}

/*
 * ======= COUNT IN SPHERE FUNCTION =========
 */

// [[Rcpp::depends(lidR)]]

#include <RcppArmadillo.h>
#include "SpatialIndex.h"
#include "Progress.h"

using namespace lidR;
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector C_count_in_sphere(S4 las, double radius, int ncpu)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  unsigned int n = X.length();
  IntegerVector ptDen_sph(n);

  lidR::SpatialIndex tree(las);                  // the index for the las object

  Progress pb(n, "Counting points in sphere: ");
  bool abort= false;

#pragma omp parallel for num_threads(ncpu)
  for(unsigned int i = 0 ; i < n ; i++)
  {
    if (abort) continue;
    if (pb.check_interrupt()) abort = true;
    pb.increment();

    std::vector<PointXYZ> sphpts;               // creation of an STL container of points for the sphere neighborhood object
    Sphere sphere(X[i], Y[i], Z[i], radius);    // creation of a sphere object
    tree.lookup(sphere, sphpts);                // lookup the points in the sphere neighborhood

#pragma omp critical
{
  ptDen_sph[i] = sphpts.size();               // count the points in the sphere neighborhood
}
  }
  if (abort) throw Rcpp::internal::InterruptedException();
  return ptDen_sph;
}
