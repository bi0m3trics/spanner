// Updated geometric_features_simple.cpp
// Computes full feature set (matching optimized version), supports OpenMP parallelism

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "nanoflann.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <atomic>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

struct PointCloud {
  std::vector<double> x, y, z;

  inline size_t kdtree_get_point_count() const { return x.size(); }

  inline double kdtree_distance(const double* p1, const size_t idx_p2, size_t) const {
    const double d0 = p1[0] - x[idx_p2];
    const double d1 = p1[1] - y[idx_p2];
    const double d2 = p1[2] - z[idx_p2];
    return d0*d0 + d1*d1 + d2*d2;
  }

  inline double kdtree_get_pt(const size_t idx, int dim) const {
    if (dim == 0) return x[idx];
    else if (dim == 1) return y[idx];
    else return z[idx];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
  nanoflann::L2_Simple_Adaptor<double, PointCloud>,
  PointCloud, 3> KDTree;

struct GeometricFeatures {
  double linearity = NA_REAL, planarity = NA_REAL, sphericity = NA_REAL, omnivariance = NA_REAL;
  double anisotropy = NA_REAL, eigenentropy = NA_REAL, surface_variation = NA_REAL, verticality = NA_REAL;
  double e_largest = NA_REAL, e_medium = NA_REAL, e_smallest = NA_REAL;
  double nx = NA_REAL, ny = NA_REAL, nz = NA_REAL;
  double roughness = NA_REAL, signed_roughness = NA_REAL;
  double num_neighbors = NA_REAL, surface_density = NA_REAL, volume_density = NA_REAL;
  double first_order_moment = NA_REAL;
  double height_above_ground = NA_REAL, relative_height = NA_REAL;
};

GeometricFeatures computeGeometricFeatures(
    size_t point_idx,
    const std::vector<nanoflann::ResultItem<unsigned int, double>>& indices_dists,
    const PointCloud& cloud,
    double radius) {

  GeometricFeatures f;
  size_t k = indices_dists.size();
  if (k < 3) return f;

  arma::mat pts(k, 3);
  double min_z = std::numeric_limits<double>::max(), max_z = -min_z;
  arma::vec3 centroid = arma::zeros<arma::vec>(3);

  for (size_t i = 0; i < k; i++) {
    size_t idx = indices_dists[i].first;
    double xi = cloud.x[idx], yi = cloud.y[idx], zi = cloud.z[idx];
    pts(i, 0) = xi; pts(i, 1) = yi; pts(i, 2) = zi;
    centroid(0) += xi; centroid(1) += yi; centroid(2) += zi;
    if (zi < min_z) min_z = zi;
    if (zi > max_z) max_z = zi;
  }
  centroid /= k;
  arma::mat centered = pts.each_row() - centroid.t();
  arma::mat cov = centered.t() * centered / (k - 1);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, cov);
  eigval = arma::sort(eigval, "descend");
  arma::uvec idx = arma::sort_index(eigval, "descend");
  eigvec = eigvec.cols(idx);

  double e1 = eigval(0), e2 = eigval(1), e3 = eigval(2);
  double sum_eigs = e1 + e2 + e3;
  if (sum_eigs < 1e-10) return f;

  f.e_largest = e1; f.e_medium = e2; f.e_smallest = e3;
  f.linearity = (e1 - e2) / e1;
  f.planarity = (e2 - e3) / e1;
  f.sphericity = e3 / e1;
  f.omnivariance = std::pow(e1 * e2 * e3, 1.0 / 3.0);
  f.anisotropy = (e1 - e3) / e1;
  f.eigenentropy = -(e1 * std::log(e1) + e2 * std::log(e2) + e3 * std::log(e3));
  f.surface_variation = e3 / sum_eigs;

  arma::vec normal = eigvec.col(2);
  f.nx = normal(0); f.ny = normal(1); f.nz = normal(2);
  f.verticality = 1.0 - std::abs(normal(2));

  arma::vec3 query = {cloud.x[point_idx], cloud.y[point_idx], cloud.z[point_idx]};
  arma::vec3 to_point = query - centroid;
  f.roughness = std::abs(arma::dot(to_point, normal));
  f.signed_roughness = arma::dot(to_point, normal);

  f.num_neighbors = k;
  double pi = M_PI;
  f.surface_density = k / (pi * radius * radius);
  f.volume_density = k / ((4.0 / 3.0) * pi * radius * radius * radius);

  double moment_sum = 0.0;
  for (size_t i = 0; i < k; i++) {
    size_t idx = indices_dists[i].first;
    arma::vec3 pt = {cloud.x[idx], cloud.y[idx], cloud.z[idx]};
    moment_sum += arma::norm(pt - centroid);
  }
  f.first_order_moment = moment_sum / k;

  double z = cloud.z[point_idx];
  f.height_above_ground = z - min_z;
  f.relative_height = (max_z > min_z) ? (z - min_z) / (max_z - min_z) : 0.0;

  return f;
}

// [[Rcpp::export]]
Rcpp::List C_geometric_features_simple(Rcpp::S4 las, double radius = 0.1, int max_neighbors = 50, int ncpu = 1, Rcpp::Nullable<Rcpp::CharacterVector> features = R_NilValue) {
  Rcpp::DataFrame data = Rcpp::as<Rcpp::DataFrame>(las.slot("data"));
  Rcpp::NumericVector x = data["X"], y = data["Y"], z = data["Z"];
  size_t n = x.size();

  PointCloud cloud;
  cloud.x.assign(x.begin(), x.end());
  cloud.y.assign(y.begin(), y.end());
  cloud.z.assign(z.begin(), z.end());

  KDTree index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
  index.buildIndex();

  NumericVector linearity(n), planarity(n), sphericity(n), omnivariance(n), anisotropy(n);
  NumericVector eigenentropy(n), surface_variation(n), verticality(n);
  NumericVector eLargest(n), eMedium(n), eSmallest(n);
  NumericVector nx(n), ny(n), nz(n);
  NumericVector roughness(n), signed_roughness(n);
  NumericVector num_neighbors(n), surface_density(n), volume_density(n);
  NumericVector first_order_moment(n), height_above_ground(n), relative_height(n);

#ifdef _OPENMP
#pragma omp parallel for num_threads(ncpu) schedule(dynamic, 100)
#endif
  for (int i = 0; i < static_cast<int>(n); ++i) {
    double query[3] = {cloud.x[i], cloud.y[i], cloud.z[i]};

    std::vector<nanoflann::ResultItem<unsigned int, double>> indices_dists;
    nanoflann::SearchParameters params;
    size_t found = index.radiusSearch(query, radius * radius, indices_dists, params);

    if (indices_dists.size() > (size_t)max_neighbors)
      indices_dists.resize(max_neighbors);

    GeometricFeatures f = computeGeometricFeatures(i, indices_dists, cloud, radius);

    linearity[i] = f.linearity;
    planarity[i] = f.planarity;
    sphericity[i] = f.sphericity;
    omnivariance[i] = f.omnivariance;
    anisotropy[i] = f.anisotropy;
    eigenentropy[i] = f.eigenentropy;
    surface_variation[i] = f.surface_variation;
    verticality[i] = f.verticality;
    eLargest[i] = f.e_largest;
    eMedium[i] = f.e_medium;
    eSmallest[i] = f.e_smallest;
    nx[i] = f.nx;
    ny[i] = f.ny;
    nz[i] = f.nz;
    roughness[i] = f.roughness;
    signed_roughness[i] = f.signed_roughness;
    num_neighbors[i] = f.num_neighbors;
    surface_density[i] = f.surface_density;
    volume_density[i] = f.volume_density;
    first_order_moment[i] = f.first_order_moment;
    height_above_ground[i] = f.height_above_ground;
    relative_height[i] = f.relative_height;
  }

  std::unordered_set<std::string> feature_set;
  if (features.isNotNull()) {
    CharacterVector fv(features);
    for (const auto& f : fv) {
      feature_set.insert(std::string(f));
    }
  }

  List out;
  auto add_if_selected = [&](std::string name, NumericVector& vec) {
    if (feature_set.empty() || feature_set.count(name)) out[name] = vec;
  };

  add_if_selected("linearity", linearity);
  add_if_selected("planarity", planarity);
  add_if_selected("sphericity", sphericity);
  add_if_selected("omnivariance", omnivariance);
  add_if_selected("anisotropy", anisotropy);
  add_if_selected("eigenentropy", eigenentropy);
  add_if_selected("surface_variation", surface_variation);
  add_if_selected("verticality", verticality);
  add_if_selected("eLargest", eLargest);
  add_if_selected("eMedium", eMedium);
  add_if_selected("eSmallest", eSmallest);
  add_if_selected("nx", nx);
  add_if_selected("ny", ny);
  add_if_selected("nz", nz);
  add_if_selected("roughness", roughness);
  add_if_selected("signed_roughness", signed_roughness);
  add_if_selected("num_neighbors", num_neighbors);
  add_if_selected("surface_density", surface_density);
  add_if_selected("volume_density", volume_density);
  add_if_selected("first_order_moment", first_order_moment);
  add_if_selected("height_above_ground", height_above_ground);
  add_if_selected("relative_height", relative_height);

  return out;
}
