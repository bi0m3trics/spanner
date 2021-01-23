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
NumericVector C_vert_in_sphere(S4 las, double radius, int ncpu)
{
  DataFrame data = as<DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  int n = X.size();
  NumericVector vert_sph(n);
  
  // One copy for the index. The index automatically detects the TLS tag in the LAS object
  lidR::GridPartition tree(las);                  
  
  // Initiate the progress bar
  Progress pb(n, "Calculating Verticality: ");
  
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
    
    // Memory Allocaitons
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
      // Calculate verticality
      // double linearity = (latent[0] - latent[1]) / latent[0];
      double verticality = 1-abs(coeff(2,2));
      vert_sph[i] = verticality;
    }
  }
  
  if (abort) throw Rcpp::internal::InterruptedException();
  
  // return the numericvector
  return vert_sph;
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