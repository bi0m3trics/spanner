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
