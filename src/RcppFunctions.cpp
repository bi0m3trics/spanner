/*
 * ======= CALC VERTICALITY IN SPHERE FUNCTION =========
 */
// [[Rcpp::depends(lidR)]]

#include <RcppArmadillo.h>
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
  int n_metrics = 12;

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
  
  List out(n_metrics);
  
  // One copy for the index. The index automatically detects the TLS tag in the LAS object
  lidR::GridPartition tree(las);                  
  
  // Initiate the progress bar
  Progress pb(n, "Calculating eigen metrics: ");
  
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
    
    // Memory Allocations
    arma::mat A(k,3);
    arma::mat coeff;        // Principle component matrix
    arma::mat score;
    arma::vec latent;       // Eigenvalues in descending order
    
    // Eigen decomposition is the real cost of the function
    for (unsigned int j = 0 ; j < sphpts.size() ; j++)
    {
      A(j,0) = sphpts[j].x;
      A(j,1) = sphpts[j].y;
      A(j,2) = sphpts[j].z;
    }
    arma::princomp(coeff, score, latent, A);
    
    #pragma omp critical
    {
      // Calculate decomposed eigen values
    
      double eigen_largest  = latent[0];
      double eigen_medium   = latent[1];
      double eigen_smallest = latent[2];
      double eigensum       = latent[0] + latent[1] + latent[2];
      double curvature      = eigen_smallest/eigensum;
      double omnivariance   = pow((latent[0] * latent[1] * latent[2]), (1.0/3.0));
      double anisotropy     = (latent[0] - latent[2]) / latent[0];
      double eigentropy     = -((latent[0] * log(latent[0])) + (latent[1] * log(latent[1])) + (latent[2] * log(latent[2])));
      double linearity      = (latent[0] - latent[1]) / latent[0];
      double verticality    = 1-abs(coeff(2,2));
      double planarity      = (latent[1] - latent[2])/latent[0];
      double sphericity     = latent[2]/latent[0];
      
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
  
  Progress pb(n, "Counting points in disc: ");
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