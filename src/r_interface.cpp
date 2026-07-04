//  ===============================================================================
//
//  Developers:
//
//  Tiago de Conto - tdc.florestal@gmail.com -  https://github.com/tiagodc/
//
//  COPYRIGHT: Tiago de Conto, 2020
//
//  This piece of software is open and free to use, redistribution and modifications
//  should be done in accordance to the GNU General Public License >= 3
//
//  Use this software as you wish, but no warranty is provided whatsoever. For any
//  comments or questions on TreeLS, please contact the developer (prefereably through my github account)
//
//  If publishing any work/study/research that used the tools in TreeLS,
//  please don't forget to cite the proper sources!
//
//  Enjoy!
//
//  ===============================================================================

#include "methods.h"
#include "pcv.h"
#include "ssao.h"

// export tree positions point stack
List exportTreeMap(vector<HoughCenters>& coordinates){

  vector<double> xout;
  vector<double> yout;
  vector<double> zout;
  vector<double> radii;
  vector<bool> keyFlag;
  vector<bool> treeFlag;
  vector<unsigned short int> votes;
  vector<unsigned int> treeId;
  vector<unsigned int> discId;
  // vector<unsigned int> nPoints;

  unsigned int diskCounter = 1;
  unsigned int maxId = 0;

  for(auto& point : coordinates){

    if(point.tree_id == 0)
      continue;

    if(point.tree_id > maxId)
      maxId = point.tree_id;

    double z = (point.low_z + point.up_z)/2;

    // main point
    xout.push_back(point.main_circle.x_center);
    yout.push_back(point.main_circle.y_center);
    zout.push_back(z);

    radii.push_back(point.main_circle.radius);
    votes.push_back(point.main_circle.n_votes);
    treeId.push_back(point.tree_id);
    // nPoints.push_back(point.circles.size());
    discId.push_back(diskCounter);
    keyFlag.push_back(true);
    treeFlag.push_back(false);

    // other candidates
    for(auto& c_point : point.circles){

      xout.push_back(c_point.x_center);
      yout.push_back(c_point.y_center);
      zout.push_back(z);

      radii.push_back(c_point.radius);
      votes.push_back(c_point.n_votes);
      treeId.push_back(point.tree_id);
      // nPoints.push_back(point.circles.size());
      discId.push_back(diskCounter);
      keyFlag.push_back(false);
      treeFlag.push_back(false);
    }
    diskCounter++;
  }

  vector<double> xSums(maxId, 0);
  vector<double> ySums(maxId, 0);
  vector<unsigned int> counters(maxId, 0);

  for(auto& point : coordinates){

    if(point.tree_id == 0)
      continue;

    xSums[point.tree_id-1] += point.main_circle.x_center;
    ySums[point.tree_id-1] += point.main_circle.y_center;
    counters[point.tree_id-1]++;
  }

  for(unsigned int i = 0; i < maxId; ++i){

    if(counters[i] == 0)
      continue;

    double mainX = xSums[i] / counters[i];
    double mainY = ySums[i] / counters[i];

    xout.push_back(mainX);
    yout.push_back(mainY);
    zout.push_back(0);

    radii.push_back(0);
    votes.push_back(0);
    treeId.push_back(i+1);
    // nPoints.push_back(0);
    discId.push_back(0);
    keyFlag.push_back(true);
    treeFlag.push_back(true);
  }

  bool hasIds = maxId > 0;

  List out;
  out["X"] = xout;
  xout.clear();
  xout.shrink_to_fit();

  out["Y"] = yout;
  yout.clear();
  yout.shrink_to_fit();

  out["Z"] = zout;
  zout.clear();
  zout.shrink_to_fit();

  out["Intensity"] = votes;
  votes.clear();
  votes.shrink_to_fit();

  out["Radius"] = radii;
  radii.clear();
  radii.shrink_to_fit();

  out["DiscID"] = discId;
  discId.clear();
  discId.shrink_to_fit();

  out["Keypoint"] = keyFlag;
  keyFlag.clear();
  keyFlag.shrink_to_fit();

  out["TreePosition"] = treeFlag;
  treeFlag.clear();
  treeFlag.shrink_to_fit();

  if(hasIds){
    out["TreeID"] = treeId;
    treeId.clear();
    treeId.shrink_to_fit();
  }

  return out;

}

vector<vector<vector<double> > > getChunks(vector<vector<double> >& las, vector<unsigned int> ids){

  unordered_set<unsigned int> treeIds(ids.begin(), ids.end());
  unsigned int nTrees = treeIds.size();

  vector<vector<vector<double> > > output(nTrees);
  for(auto& trees : output){
    trees.resize(3);
  }

  unordered_map<unsigned int, unsigned int> idMap;
  unsigned int counter = 0;
  for(auto& i : treeIds){
    idMap[i] = counter;
    counter++;
  }

  for(unsigned int i = 0; i < las[0].size(); ++i){

    unsigned int& id = ids[i];

    if(id == 0)
      continue;

    unsigned int tempId = idMap[id];

    output[tempId][0].push_back( las[0][i] );
    output[tempId][1].push_back( las[1][i] );
    output[tempId][2].push_back( las[2][i] );
  }

  return output;

}

// [[Rcpp::export]]
SEXP cppCylinderFit(NumericMatrix& las, std::string method = "nm", unsigned int n = 10, double p = 0.95, double inliers = 0.9, double max_angle = 30, unsigned int n_best = 20){
  vector<vector<double> > cloud = rmatrix2cpp(las);
  vector<double> pars;

  double nmax = 100;
  if(method != "ransac" && method != "bf" && cloud[0].size() > nmax){
    double prop = nmax / (double)cloud[0].size();
    cloud = randomPoints(cloud, prop);
  }

  if(method == "irls"){
    vector<double> initPars = {0, M_PI/2, 0, 0, 0};
    pars = irlsCylinder(cloud, initPars);
  }else if(method == "nm"){
    pars = nmCylinderFit(cloud);
  }else if(method == "ransac"){
    pars = ransacCylinder(cloud, n, p, inliers);
  }else if(method == "bf"){
    pars = bruteForceRansacCylinder(cloud, n, p, inliers, n_best, max_angle, true)[0];
  }

  return wrap( pars );
}

// [[Rcpp::export]]
SEXP cppComputePCV(S4 las, double radius = 1.0, int num_directions = 60, int ncpu = 1){
  // Call PCV computation with S4 LAS object
  vector<double> pcv_values = computePCV(las, radius, num_directions, ncpu);
  
  return wrap(pcv_values);
}

// [[Rcpp::export]]
NumericVector cppComputeSSAO(S4 las, int kernel_size = 5, double pixel_size = 0.1,
                             int num_samples = 16, int ncpu = 4) {
  // Call SSAO computation with S4 LAS object
  vector<double> ssao_values = computeSSAO(las, kernel_size, pixel_size, num_samples, ncpu);
  
  return wrap(ssao_values);
}

// ============================================================================
//  spanner v2.0 additions
//  New C++ entry points for the explicit staged TLS/MLS pipeline.
// ============================================================================

// ----------------------------------------------------------------------------
// C_stack_map – Plot-level Hough tree mapping across stacked height slices.
//
// Iterates over non-overlapping horizontal slices [h_min, h_max) and runs the
// circular Hough Transform on each.  Disks that are vertically persistent
// across >= 75 % of slices (min_layers) are assigned sequential tree IDs.
//
// Parameters
//   las          NumericMatrix  N x 3 (X, Y, Z)
//   h_min        double         Lower slice boundary (m)
//   h_max        double         Upper slice boundary (m)
//   h_step       double         Slice thickness (m)
//   pixel_size   double         Hough accumulator pixel size (m)
//   max_radius   double         Maximum circle radius searched (m)
//   min_density  double [0,1]   Min fraction of max-count for a pixel to vote
//   min_votes    uint           Min Hough votes to accept a circle centre
//
// Returns an R List with the same structure as exportTreeMap():
//   X, Y, Z, Intensity (votes), Radius, DiscID, Keypoint, TreePosition, TreeID
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List C_stack_map(NumericMatrix las,
                 double       h_min       = 1.0,
                 double       h_max       = 3.0,
                 double       h_step      = 0.5,
                 double       pixel_size  = 0.025,
                 double       max_radius  = 0.25,
                 double       min_density = 0.1,
                 unsigned int min_votes   = 3) {

  if (h_max <= h_min)
    stop("C_stack_map: h_max must be greater than h_min");
  if (h_step <= 0 || pixel_size <= 0 || max_radius <= 0)
    stop("C_stack_map: h_step, pixel_size and max_radius must be positive");
  if (min_density <= 0 || min_density > 1)
    stop("C_stack_map: min_density must be in (0, 1]");

  vector<vector<double>> cloud = rmatrix2cpp(las);

  unsigned int n_slices = static_cast<unsigned int>(
      std::max(1.0, std::ceil((h_max - h_min) / h_step)));

  vector<vector<vector<double>>> slices = getSlices(cloud, h_min, h_max, h_step);
  cloud.clear();
  cloud.shrink_to_fit();

  vector<HoughCenters> all_disks;
  all_disks.reserve(slices.size() * 4);  // rough upper bound

  for (unsigned int i = 0; i < slices.size(); ++i) {
    if (slices[i][0].empty()) continue;
    Raster raster = getCounts(slices[i], pixel_size);
    // getCenters returns all candidate circle centres in this slice
    vector<HoughCenters> slice_disks = getCenters(&raster, max_radius,
                                                   min_density, min_votes);
    double z_lo = h_min + i * h_step;
    double z_hi = z_lo + h_step;
    for (auto& disk : slice_disks) {
      disk.low_z = z_lo;
      disk.up_z  = z_hi;
      all_disks.push_back(std::move(disk));
    }
  }

  if (all_disks.empty())
    return List::create();

  // Require vertical persistence across >= 75 % of height slices
  unsigned int min_layers = std::max(1u,
      static_cast<unsigned int>(std::ceil(n_slices * 0.75)));
  assignTreeId(all_disks, max_radius * 2.0, min_density, min_layers);

  return exportTreeMap(all_disks);
}


// ----------------------------------------------------------------------------
// C_hough_stem_points – Per-point stem labels via slice Hough (single tree).
//
// Seeds a Hough circle in the [h_base_lo, h_base_hi] band, then uses
// treeHough() to track the cylinder axis upward.  Each point is labelled
// as a stem point if it falls within (radius + pixel_size) of its slice's
// tracked Hough centre.
//
// Returns a List: Stem (logical), Segment (uint), Radius (double), Votes (uint)
// – one element per row of 'las'.
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List C_hough_stem_points(NumericMatrix las,
                          double       h_base_lo   = 1.0,
                          double       h_base_hi   = 2.5,
                          double       h_step      = 0.5,
                          double       max_radius  = 0.25,
                          double       pixel_size  = 0.025,
                          double       min_density = 0.1,
                          unsigned int min_votes   = 3) {

  unsigned int n = static_cast<unsigned int>(las.nrow());
  vector<bool>         stem_flag(n, false);
  vector<unsigned int> seg_out(n, 0);
  vector<double>       rad_out(n, 0.0);
  vector<unsigned int> vot_out(n, 0);

  if (n == 0)
    return List::create(Named("Stem")    = stem_flag,
                        Named("Segment") = seg_out,
                        Named("Radius")  = rad_out,
                        Named("Votes")   = vot_out);

  vector<vector<double>> cloud = rmatrix2cpp(las);

  // Seed: best circle in the base height range
  vector<vector<vector<double>>> base_slices =
      getSlices(cloud, h_base_lo, h_base_hi, h_base_hi - h_base_lo);

  if (base_slices.empty() || base_slices[0][0].empty())
    return List::create(Named("Stem")    = stem_flag,
                        Named("Segment") = seg_out,
                        Named("Radius")  = rad_out,
                        Named("Votes")   = vot_out);

  Raster base_raster = getCounts(base_slices[0], pixel_size);
  HoughCenters base_hc = getSingleCenter(&base_raster, max_radius,
                                          min_density, min_votes);
  if (base_hc.main_circle.n_votes < min_votes)
    return List::create(Named("Stem")    = stem_flag,
                        Named("Segment") = seg_out,
                        Named("Radius")  = rad_out,
                        Named("Votes")   = vot_out);

  // Track circle upward from base through full cloud height
  vector<HoughCenters> hc_stack = treeHough(cloud,
      h_base_lo, h_base_hi, h_step, max_radius, pixel_size,
      min_density, min_votes);

  // Label each point according to the tracked Hough circle at its height
  for (unsigned int i = 0; i < n; ++i) {
    double x = cloud[0][i];
    double y = cloud[1][i];
    double z = cloud[2][i];
    if (z < 0.0) continue;

    unsigned int seg_idx = static_cast<unsigned int>(std::floor(z / h_step));
    if (seg_idx >= hc_stack.size()) continue;

    const HoughCircle& hc = hc_stack[seg_idx].main_circle;
    if (hc.n_votes < min_votes) continue;

    double dist = std::sqrt(std::pow(x - hc.x_center, 2) +
                            std::pow(y - hc.y_center, 2));
    if (dist <= hc.radius + pixel_size) {
      stem_flag[i] = true;
      seg_out[i]   = seg_idx + 1u;
      rad_out[i]   = hc.radius;
      vot_out[i]   = hc.n_votes;
    }
  }

  return List::create(Named("Stem")    = stem_flag,
                      Named("Segment") = seg_out,
                      Named("Radius")  = rad_out,
                      Named("Votes")   = vot_out);
}


// ----------------------------------------------------------------------------
// C_hough_stem_plot – Per-point stem labels for a multi-tree plot.
//
// Iterates over unique non-zero TreeIDs in 'tree_ids', extracts the per-tree
// sub-cloud, and calls the single-tree Hough stem detector.  Results are
// assembled into plot-level output vectors aligned with the rows of 'las'.
//
// Parameters
//   las       NumericMatrix   N x 3 (X, Y, Z)
//   tree_ids  IntegerVector   TreeID per point (0 = unassigned / ground)
//   ... same scalar parameters as C_hough_stem_points
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List C_hough_stem_plot(NumericMatrix las,
                        IntegerVector tree_ids,
                        double       h_base_lo   = 1.0,
                        double       h_base_hi   = 2.5,
                        double       h_step      = 0.5,
                        double       max_radius  = 0.25,
                        double       pixel_size  = 0.025,
                        double       min_density = 0.1,
                        unsigned int min_votes   = 3) {

  unsigned int n = static_cast<unsigned int>(las.nrow());
  vector<bool>         stem_flag(n, false);
  vector<unsigned int> seg_out(n, 0);
  vector<double>       rad_out(n, 0.0);
  vector<unsigned int> vot_out(n, 0);

  if (n == 0)
    return List::create(Named("Stem")    = stem_flag,
                        Named("Segment") = seg_out,
                        Named("Radius")  = rad_out,
                        Named("Votes")   = vot_out);

  // Collect unique non-zero tree IDs
  std::unordered_set<int> id_set;
  for (int i = 0; i < static_cast<int>(n); ++i)
    if (tree_ids[i] > 0) id_set.insert(tree_ids[i]);

  if (id_set.empty())
    return List::create(Named("Stem")    = stem_flag,
                        Named("Segment") = seg_out,
                        Named("Radius")  = rad_out,
                        Named("Votes")   = vot_out);

  vector<vector<double>> full_cloud = rmatrix2cpp(las);

  for (int tid : id_set) {
    // Collect indices of points belonging to this tree
    vector<unsigned int> pt_idx;
    for (unsigned int i = 0; i < n; ++i)
      if (tree_ids[i] == tid) pt_idx.push_back(i);
    if (pt_idx.empty()) continue;

    // Extract per-tree sub-cloud (avoid full copy of full_cloud per tree)
    unsigned int m = static_cast<unsigned int>(pt_idx.size());
    vector<vector<double>> tree_cloud(3, vector<double>(m));
    for (unsigned int k = 0; k < m; ++k) {
      tree_cloud[0][k] = full_cloud[0][pt_idx[k]];
      tree_cloud[1][k] = full_cloud[1][pt_idx[k]];
      tree_cloud[2][k] = full_cloud[2][pt_idx[k]];
    }

    // Seed in base band
    vector<vector<vector<double>>> base_sl =
        getSlices(tree_cloud, h_base_lo, h_base_hi, h_base_hi - h_base_lo);
    if (base_sl.empty() || base_sl[0][0].empty()) continue;

    Raster bras = getCounts(base_sl[0], pixel_size);
    HoughCenters bhc = getSingleCenter(&bras, max_radius, min_density, min_votes);
    if (bhc.main_circle.n_votes < min_votes) continue;

    // Track upward
    vector<HoughCenters> hc_stack = treeHough(tree_cloud,
        h_base_lo, h_base_hi, h_step, max_radius, pixel_size,
        min_density, min_votes);

    // Label
    for (unsigned int k = 0; k < m; ++k) {
      double x = tree_cloud[0][k];
      double y = tree_cloud[1][k];
      double z = tree_cloud[2][k];
      if (z < 0.0) continue;

      unsigned int seg_idx = static_cast<unsigned int>(std::floor(z / h_step));
      if (seg_idx >= hc_stack.size()) continue;

      const HoughCircle& hc = hc_stack[seg_idx].main_circle;
      if (hc.n_votes < min_votes) continue;

      double dist = std::sqrt(std::pow(x - hc.x_center, 2) +
                              std::pow(y - hc.y_center, 2));
      if (dist <= hc.radius + pixel_size) {
        unsigned int orig = pt_idx[k];
        stem_flag[orig] = true;
        seg_out[orig]   = seg_idx + 1u;
        rad_out[orig]   = hc.radius;
        vot_out[orig]   = hc.n_votes;
      }
    }
  }

  return List::create(Named("Stem")    = stem_flag,
                      Named("Segment") = seg_out,
                      Named("Radius")  = rad_out,
                      Named("Votes")   = vot_out);
}


// ----------------------------------------------------------------------------
// C_tree_ids_voronoi – Nearest-tree (Voronoi) point-to-tree assignment.
//
// Each point receives the ID of its closest tree centre in 2-D Euclidean
// space.  Points farther than max_dist from every centre are assigned 0.
// Pass max_dist = -1 to assign ALL points (pure Voronoi partition).
//
// Returns an IntegerVector of tree IDs, length == pt_x.size().
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector C_tree_ids_voronoi(NumericVector pt_x,
                                  NumericVector pt_y,
                                  NumericVector tree_x,
                                  NumericVector tree_y,
                                  IntegerVector tree_ids,
                                  double max_dist = -1.0) {

  int n_pts   = pt_x.size();
  int n_trees = tree_x.size();
  IntegerVector result(n_pts, 0);
  if (n_trees == 0) return result;

  double max_d2 = (max_dist > 0.0) ? max_dist * max_dist
                                    : std::numeric_limits<double>::max();

  for (int i = 0; i < n_pts; ++i) {
    double px = pt_x[i], py = pt_y[i];
    double best_d2 = max_d2;
    int    best_j  = -1;
    for (int j = 0; j < n_trees; ++j) {
      double dx = px - tree_x[j];
      double dy = py - tree_y[j];
      double d2 = dx * dx + dy * dy;
      if (d2 < best_d2) { best_d2 = d2; best_j = j; }
    }
    if (best_j >= 0) result[i] = tree_ids[best_j];
  }
  return result;
}


// ----------------------------------------------------------------------------
// C_tree_ids_crop – Fixed-radius (or square) crop point-to-tree assignment.
//
// A point is assigned to the first tree centre whose crop contains it.
// Where crops overlap, the tree with the lowest index in 'tree_ids' wins.
// Points outside every crop receive 0.
//
// Parameters
//   length   radius (circle=TRUE) or half-side (circle=FALSE)
//   circle   TRUE = circular crop, FALSE = square crop
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector C_tree_ids_crop(NumericVector pt_x,
                               NumericVector pt_y,
                               NumericVector tree_x,
                               NumericVector tree_y,
                               IntegerVector tree_ids,
                               double length = 1.0,
                               bool   circle = true) {

  int n_pts   = pt_x.size();
  int n_trees = tree_x.size();
  IntegerVector result(n_pts, 0);
  if (n_trees == 0) return result;

  double r2 = length * length;

  for (int i = 0; i < n_pts; ++i) {
    double px = pt_x[i], py = pt_y[i];
    for (int j = 0; j < n_trees; ++j) {
      double dx = px - tree_x[j];
      double dy = py - tree_y[j];
      bool inside = circle ? (dx * dx + dy * dy < r2)
                           : (std::abs(dx) < length && std::abs(dy) < length);
      if (inside) { result[i] = tree_ids[j]; break; }
    }
  }
  return result;
}


// ----------------------------------------------------------------------------
// C_ransac_plot_cylinders – Batch RANSAC cylinder fitting for all trees x slices.
//
// Thin wrapper around ransacPlotCylinders() (same algorithm used by TreeLS).
// All trees and height-slices are processed in a single C++ call, eliminating
// the R-level double loop.
//
// Parameters
//   las          NumericMatrix   N x 3 (X, Y, Z)
//   tId          NumericVector   per-point tree ID (integer-valued)
//   segs         NumericVector   per-point segment index (0-based integer)
//   rads         NumericVector   per-point initial radius estimate (m)
//   n_samples    uint            RANSAC samples per iteration
//   conf         double          RANSAC confidence level
//   inlier_frac  double          expected inlier fraction
//   tol          double          inlier distance tolerance (m)
//
// Returns a nested R List: [[tree_idx]][[seg_idx]] = numeric cylinder params.
// Cylinder vector: [phi, theta, x_ctr, y_ctr, radius, ..., seg_id, tree_id]
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List C_ransac_plot_cylinders(NumericMatrix las,
                              NumericVector tId,
                              NumericVector segs,
                              NumericVector rads,
                              unsigned int n_samples   = 10,
                              double       conf        = 0.95,
                              double       inlier_frac = 0.9,
                              double       tol         = 0.05) {
  vector<vector<double>>  cloud(rmatrix2cpp(las));
  vector<unsigned int>    treeId(tId.begin(),  tId.end());
  vector<unsigned int>    segments(segs.begin(), segs.end());
  vector<double>          radii(rads.begin(),  rads.end());
  return wrap(ransacPlotCylinders(cloud, treeId, segments, radii,
                                  n_samples, conf, inlier_frac, tol));
}


// ----------------------------------------------------------------------------
// C_irls_plot_cylinders – Batch IRLS cylinder fitting for all trees x slices.
//
// Thin wrapper around irlsPlotCylinders(). Faster than RANSAC for clean data;
// less robust to heavy outliers.  All trees x slices in one C++ call.
//
// Parameters
//   las       NumericMatrix   N x 3 (X, Y, Z)
//   tId       NumericVector   per-point tree ID
//   segs      NumericVector   per-point segment index (0-based)
//   rads      NumericVector   per-point initial radius estimate (m)
//   n_points  uint            max pts per slice passed to IRLS (0 = all)
//   tol       double          inlier distance tolerance (m)
// ----------------------------------------------------------------------------
// [[Rcpp::export]]
List C_irls_plot_cylinders(NumericMatrix las,
                            NumericVector tId,
                            NumericVector segs,
                            NumericVector rads,
                            unsigned int n_points = 100,
                            double       tol      = 0.05) {
  vector<vector<double>>  cloud(rmatrix2cpp(las));
  vector<unsigned int>    treeId(tId.begin(),  tId.end());
  vector<unsigned int>    segments(segs.begin(), segs.end());
  vector<double>          radii(rads.begin(),  rads.end());
  return wrap(irlsPlotCylinders(cloud, treeId, segments, radii, n_points, tol));
}
