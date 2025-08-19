/*
 * ======= CALC EIGEN METRICS IN SPHERE FUNCTION =========
 */
// [[Rcpp::depends(lidR)]]

#include <RcppArmadillo.h>
#include <limits>
#include "SpatialIndex.h"
#include "Progress.h"

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
  int n_metrics = 15;

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

    // Only perform eigenvalue decomposition if we have enough points
    if (k >= 3) {
      // Compute centroid
      double cx = 0.0, cy = 0.0, cz = 0.0;
      for (unsigned int j = 0; j < sphpts.size(); j++) {
        cx += sphpts[j].x;
        cy += sphpts[j].y;
        cz += sphpts[j].z;
      }
      cx /= k; cy /= k; cz /= k;
      
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
  
  // CloudCompare-style eigentropy calculation using raw eigenvalues (exact match)
  if (eigen_largest > 0.0 && eigen_medium > 0.0 && eigen_smallest > 0.0) {
    eigentropy = -((eigen_largest * log(eigen_largest)) + (eigen_medium * log(eigen_medium)) + (eigen_smallest * log(eigen_smallest)));
  } else {
    eigentropy = 0.0;
  }
  
  linearity      = (eigen_largest - eigen_medium) / eigen_largest;
  
  // Get the eigenvector corresponding to the smallest eigenvalue (normal vector)
  arma::vec normal = eigenvectors.col(sorted_indices[2]);
  verticality    = 1 - abs(normal[2]);
  planarity      = (eigen_medium - eigen_smallest) / eigen_largest;
  sphericity     = eigen_smallest / eigen_largest;
  nx = normal[0];
  ny = normal[1];
  nz = normal[2];
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
