#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

inline long long key2d(int ix, int iy) {
  // Combine signed 32-bit grid coordinates without shifting negative values.
  const std::uint64_t ux = static_cast<std::uint32_t>(ix);
  const std::uint64_t uy = static_cast<std::uint32_t>(iy);
  return static_cast<long long>((ux << 32) | uy);
}

inline long long key3d(int ix, int iy, int iz) {
  // Keep hash arithmetic in unsigned space to avoid signed-overflow UB.
  const std::uint64_t ux = static_cast<std::uint32_t>(ix);
  const std::uint64_t uy = static_cast<std::uint32_t>(iy);
  const std::uint64_t uz = static_cast<std::uint32_t>(iz);
  std::uint64_t k = ux * 73856093ULL;
  k ^= uy * 19349663ULL;
  k ^= uz * 83492791ULL;
  return static_cast<long long>(k);
}

inline double sqr(double x) {
  return x * x;
}

inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

double allometry_profile(double q, int profile_id) {
  if (!std::isfinite(q)) return 0.05;
  if (q <= 0.0 || q >= 1.0) return 0.05;

  if (profile_id == 2) {
    double v = 1.0 - sqr((q - 0.62) / 0.45);
    return std::max(0.05, std::min(1.0, v));
  }

  if (profile_id == 3) {
    double v = std::min(1.0, std::max(0.05, q));
    return v;
  }

  // Default: smooth beta-like bell profile.
  double v = std::exp(-sqr((q - 0.62) / 0.30));
  return std::max(0.05, std::min(1.0, v));
}

double allometry_radius(double z,
                        double H,
                        double crown_a,
                        double crown_b,
                        double max_crown_radius,
                        int profile_id) {
  if (!std::isfinite(H) || H <= 0.0) return 0.0;
  double q = z / H;
  double base = crown_a * std::pow(H, crown_b);
  double r = base * allometry_profile(q, profile_id);
  if (std::isfinite(max_crown_radius) && max_crown_radius > 0.0) {
    r = std::min(r, max_crown_radius);
  }
  return std::max(0.0, r);
}

struct Grid2D {
  double cell;
  std::unordered_map<long long, std::vector<int> > bins;

  Grid2D(const NumericVector& x, const NumericVector& y, double cell_size)
    : cell(cell_size > 0 ? cell_size : 1.0) {
    const int n = x.size();
    bins.reserve(static_cast<size_t>(n * 1.3));
    for (int i = 0; i < n; ++i) {
      int ix = static_cast<int>(std::floor(x[i] / cell));
      int iy = static_cast<int>(std::floor(y[i] / cell));
      bins[key2d(ix, iy)].push_back(i);
    }
  }

  std::vector<int> query_radius(const NumericVector& x,
                                const NumericVector& y,
                                double x0,
                                double y0,
                                double radius) const {
    std::vector<int> out;
    if (radius <= 0) return out;

    int ix0 = static_cast<int>(std::floor(x0 / cell));
    int iy0 = static_cast<int>(std::floor(y0 / cell));
    int span = static_cast<int>(std::ceil(radius / cell));
    double r2 = radius * radius;

    for (int dx = -span; dx <= span; ++dx) {
      for (int dy = -span; dy <= span; ++dy) {
        auto it = bins.find(key2d(ix0 + dx, iy0 + dy));
        if (it == bins.end()) continue;
        const std::vector<int>& idx = it->second;
        for (size_t k = 0; k < idx.size(); ++k) {
          int j = idx[k];
          double d2 = sqr(x[j] - x0) + sqr(y[j] - y0);
          if (d2 <= r2) out.push_back(j);
        }
      }
    }

    return out;
  }
};

std::vector< std::vector<int> > build_knn(const NumericVector& x,
                                          const NumericVector& y,
                                          int k,
                                          double cell_size,
                                          double init_radius,
                                          int expand_steps = 4) {
  int n = x.size();
  std::vector< std::vector<int> > neigh(static_cast<size_t>(n));
  Grid2D grid(x, y, cell_size);

  for (int i = 0; i < n; ++i) {
    double radius = std::max(init_radius, cell_size);
    std::vector<int> cand;

    for (int step = 0; step < expand_steps; ++step) {
      cand = grid.query_radius(x, y, x[i], y[i], radius);
      if (static_cast<int>(cand.size()) > k) break;
      radius *= 1.75;
    }

    std::vector< std::pair<double, int> > dist_idx;
    dist_idx.reserve(cand.size());

    for (size_t t = 0; t < cand.size(); ++t) {
      int j = cand[t];
      if (j == i) continue;
      double d2 = sqr(x[j] - x[i]) + sqr(y[j] - y[i]);
      dist_idx.push_back(std::make_pair(d2, j));
    }

    if (dist_idx.empty()) {
      neigh[static_cast<size_t>(i)] = std::vector<int>();
      continue;
    }

    std::sort(dist_idx.begin(), dist_idx.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                if (a.first == b.first) return a.second < b.second;
                return a.first < b.first;
              });

    int take = std::min(k, static_cast<int>(dist_idx.size()));
    neigh[static_cast<size_t>(i)].reserve(static_cast<size_t>(take));
    for (int m = 0; m < take; ++m) {
      neigh[static_cast<size_t>(i)].push_back(dist_idx[static_cast<size_t>(m)].second);
    }
  }

  return neigh;
}

std::vector<int> local_density(const NumericVector& x,
                               const NumericVector& y,
                               double radius,
                               double cell_size) {
  int n = x.size();
  std::vector<int> density(static_cast<size_t>(n), 0);
  Grid2D grid(x, y, cell_size);

  for (int i = 0; i < n; ++i) {
    std::vector<int> cand = grid.query_radius(x, y, x[i], y[i], radius);
    density[static_cast<size_t>(i)] = static_cast<int>(cand.size()) - 1;
    if (density[static_cast<size_t>(i)] < 0) density[static_cast<size_t>(i)] = 0;
  }

  return density;
}

std::vector<int> map_seeds_to_points(const NumericVector& x,
                                     const NumericVector& y,
                                     const NumericVector& z,
                                     const NumericVector& seed_x,
                                     const NumericVector& seed_y,
                                     double hmin) {
  int n = x.size();
  int s = seed_x.size();
  std::vector<int> idx(static_cast<size_t>(s), -1);

  for (int i = 0; i < s; ++i) {
    double best = std::numeric_limits<double>::infinity();
    int best_j = -1;

    for (int j = 0; j < n; ++j) {
      if (z[j] < hmin) continue;
      double d2 = sqr(x[j] - seed_x[i]) + sqr(y[j] - seed_y[i]);
      if (d2 < best) {
        best = d2;
        best_j = j;
      }
    }

    if (best_j < 0) {
      for (int j = 0; j < n; ++j) {
        double d2 = sqr(x[j] - seed_x[i]) + sqr(y[j] - seed_y[i]);
        if (d2 < best) {
          best = d2;
          best_j = j;
        }
      }
    }

    idx[static_cast<size_t>(i)] = best_j;
  }

  return idx;
}

// Grid-based nearest-neighbor mapping from seeds to a pre-filtered point set.
// x, y contain only the above-hmin points; returned indices are in [0, x.size()).
// O(s * cells_checked) instead of O(n * s).
std::vector<int> map_seeds_to_points_grid(const NumericVector& x,
                                          const NumericVector& y,
                                          const NumericVector& seed_x,
                                          const NumericVector& seed_y) {
  const int nh = x.size();
  const int s  = seed_x.size();
  std::vector<int> idx(static_cast<size_t>(s), -1);
  if (nh == 0) return idx;

  // Coarse grid for fast proximity lookup.
  const double cell = 10.0;
  std::unordered_map<long long, std::vector<int> > bins;
  bins.reserve(static_cast<size_t>(nh / 2 + 1));
  for (int j = 0; j < nh; ++j) {
    int ix = static_cast<int>(std::floor(x[j] / cell));
    int iy = static_cast<int>(std::floor(y[j] / cell));
    bins[key2d(ix, iy)].push_back(j);
  }

  for (int i = 0; i < s; ++i) {
    double best = std::numeric_limits<double>::infinity();
    int   best_j = -1;
    const int ix0 = static_cast<int>(std::floor(seed_x[i] / cell));
    const int iy0 = static_cast<int>(std::floor(seed_y[i] / cell));

    for (int span = 0; span <= 100; ++span) {
      // Points in the next ring (span+1) are at Euclidean distance >= span*cell.
      // If we already hold a closer candidate, no further ring can beat it.
      const double next_min = static_cast<double>(span) * cell;
      if (best_j >= 0 && next_min * next_min >= best) break;

      for (int dx = -span; dx <= span; ++dx) {
        for (int dy = -span; dy <= span; ++dy) {
          // Visit only the border cells of the current ring.
          if (span > 0 && std::abs(dx) < span && std::abs(dy) < span) continue;
          const auto it = bins.find(key2d(ix0 + dx, iy0 + dy));
          if (it == bins.end()) continue;
          const std::vector<int>& pts = it->second;
          for (size_t t = 0; t < pts.size(); ++t) {
            const int jj = pts[t];
            const double d2 = sqr(x[jj] - seed_x[i]) + sqr(y[jj] - seed_y[i]);
            if (d2 < best) { best = d2; best_j = jj; }
          }
        }
      }
    }
    idx[static_cast<size_t>(i)] = best_j;
  }
  return idx;
}

std::vector<int> connected_components_points(const std::vector< std::vector<int> >& neigh,
                                             const NumericVector& x,
                                             const NumericVector& y,
                                             const NumericVector& z,
                                             double hmin,
                                             double max_jump) {
  int n = x.size();
  std::vector<int> comp(static_cast<size_t>(n), -1);
  int cid = 0;
  const double max2 = sqr(max_jump * 1.2);

  for (int i = 0; i < n; ++i) {
    if (z[i] < hmin || comp[static_cast<size_t>(i)] >= 0) continue;

    std::vector<int> stack(1, i);
    comp[static_cast<size_t>(i)] = cid;

    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();

      const std::vector<int>& nu = neigh[static_cast<size_t>(u)];
      for (size_t t = 0; t < nu.size(); ++t) {
        int v = nu[t];
        if (z[v] < hmin || comp[static_cast<size_t>(v)] >= 0) continue;

        double d2 = sqr(x[u] - x[v]) + sqr(y[u] - y[v]);
        if (d2 > max2) continue;

        comp[static_cast<size_t>(v)] = cid;
        stack.push_back(v);
      }
    }

    ++cid;
  }

  return comp;
}

void compute_local_eigenfeatures(const NumericVector& x,
                                 const NumericVector& y,
                                 const NumericVector& z,
                                 double radius,
                                 std::vector<double>& anisotropy,
                                 std::vector<double>& verticality) {
  int n = x.size();
  anisotropy.assign(static_cast<size_t>(n), 0.5);
  verticality.assign(static_cast<size_t>(n), 0.7);

  Grid2D grid(x, y, radius > 0 ? radius : 1.0);

  for (int i = 0; i < n; ++i) {
    std::vector<int> cand = grid.query_radius(x, y, x[i], y[i], radius);
    if (cand.size() < 4) continue;

    double mx = 0.0, my = 0.0, mz = 0.0;
    for (size_t t = 0; t < cand.size(); ++t) {
      int j = cand[t];
      mx += x[j];
      my += y[j];
      mz += z[j];
    }

    double inv = 1.0 / static_cast<double>(cand.size());
    mx *= inv;
    my *= inv;
    mz *= inv;

    arma::mat cov(3, 3, arma::fill::zeros);
    for (size_t t = 0; t < cand.size(); ++t) {
      int j = cand[t];
      arma::vec v(3);
      v(0) = x[j] - mx;
      v(1) = y[j] - my;
      v(2) = z[j] - mz;
      cov += v * v.t();
    }

    cov /= std::max(1.0, static_cast<double>(cand.size()) - 1.0);

    arma::vec eigval;
    arma::mat eigvec;
    bool ok = arma::eig_sym(eigval, eigvec, cov);
    if (!ok || eigval.n_elem != 3) continue;

    double l1 = std::max(1e-12, static_cast<double>(eigval(2)));
    double l3 = std::max(0.0, static_cast<double>(eigval(0)));

    double an = (l1 - l3) / l1;
    anisotropy[static_cast<size_t>(i)] = clamp01(an);

    arma::vec principal = eigvec.col(2);
    double vert = std::abs(principal(2));
    verticality[static_cast<size_t>(i)] = clamp01(vert);
  }
}

IntegerVector labels_to_integer_ids(const std::vector<int>& label,
                                    const NumericVector& z,
                                    double hmin) {
  int n = z.size();
  IntegerVector out(n, NA_INTEGER);
  for (int i = 0; i < n; ++i) {
    if (z[i] < hmin) {
      out[i] = NA_INTEGER;
    } else if (label[static_cast<size_t>(i)] > 0) {
      out[i] = label[static_cast<size_t>(i)];
    } else {
      out[i] = NA_INTEGER;
    }
  }
  return out;
}

} // namespace

extern "C" SEXP cpp_allometric_li_geodesic(SEXP xSEXP,
                                             SEXP ySEXP,
                                             SEXP zSEXP,
                                             SEXP seedXSEXP,
                                             SEXP seedYSEXP,
                                             SEXP seedZSEXP,
                                             SEXP hminSEXP,
                                             SEXP kSEXP,
                                             SEXP maxJumpSEXP,
                                             SEXP binHeightSEXP,
                                             SEXP crownASEXP,
                                             SEXP crownBSEXP,
                                             SEXP maxCrownRadiusSEXP,
                                             SEXP centerDriftSEXP,
                                             SEXP densityRadiusSEXP,
                                             SEXP gapPenaltySEXP,
                                             SEXP geodesicThresholdSEXP,
                                             SEXP useConnectedSEXP,
                                             SEXP minComponentSizeSEXP,
                                             SEXP profileIdSEXP) {
  NumericVector x(xSEXP), y(ySEXP), z(zSEXP);
  NumericVector seed_x(seedXSEXP), seed_y(seedYSEXP), seed_z(seedZSEXP);

  const int n = x.size();
  if (y.size() != n || z.size() != n) {
    stop("X, Y, Z must share length.");
  }

  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(z[i])) {
      stop("X, Y, Z contain non-finite values; clean the LAS before segmentation.");
    }
  }

  if (n == 0) return IntegerVector();

  const double hmin = as<double>(hminSEXP);
  const int k = std::max(2, as<int>(kSEXP));
  const double max_jump = std::max(1e-6, as<double>(maxJumpSEXP));
  const double bin_height = std::max(1e-6, as<double>(binHeightSEXP));
  const double crown_a = as<double>(crownASEXP);
  const double crown_b = as<double>(crownBSEXP);
  const double max_crown_radius = as<double>(maxCrownRadiusSEXP);
  const double center_drift = std::max(0.0, as<double>(centerDriftSEXP));
  const double density_radius = std::max(1e-6, as<double>(densityRadiusSEXP));
  const double gap_penalty = std::max(0.0, as<double>(gapPenaltySEXP));
  const double geodesic_threshold = std::max(1e-6, as<double>(geodesicThresholdSEXP));
  const bool use_connected_components = as<int>(useConnectedSEXP) == 1;
  const int min_component_size = std::max(1, as<int>(minComponentSizeSEXP));
  const int profile_id = as<int>(profileIdSEXP);

  const int s = seed_x.size();
  if (seed_y.size() != s || seed_z.size() != s || s == 0) {
    return IntegerVector(n, NA_INTEGER);
  }

  for (int i = 0; i < s; ++i) {
    if (!std::isfinite(seed_x[i]) || !std::isfinite(seed_y[i]) || !std::isfinite(seed_z[i])) {
      stop("Seed coordinates contain non-finite values.");
    }
  }

  // --- Height filter: build compact subset (z >= hmin) for KNN and propagation ---
  // This dramatically reduces memory for the KNN graph and avoids iterating over
  // ground/low points that would receive NA anyway.
  std::vector<int> hi_idx;
  hi_idx.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (z[i] >= hmin) hi_idx.push_back(i);
  }
  const int nh = static_cast<int>(hi_idx.size());

  if (nh == 0) return IntegerVector(n, NA_INTEGER);

  NumericVector hx(static_cast<size_t>(nh));
  NumericVector hy(static_cast<size_t>(nh));
  NumericVector hz(static_cast<size_t>(nh));
  for (int i = 0; i < nh; ++i) {
    const int j = hi_idx[static_cast<size_t>(i)];
    hx[i] = x[j];
    hy[i] = y[j];
    hz[i] = z[j];
  }

  // Build KNN, density, and connected components on the compact subset only.
  std::vector< std::vector<int> > neigh = build_knn(hx, hy, k, max_jump, max_jump * 2.0, 5);
  std::vector<int> density = local_density(hx, hy, density_radius, density_radius);
  std::vector<int> comp;
  if (use_connected_components) {
    // All hz values are >= hmin, so the hmin guard inside never fires.
    comp = connected_components_points(neigh, hx, hy, hz, hmin, max_jump);
  } else {
    comp.assign(static_cast<size_t>(nh), 0);
  }

  // Grid-based seed-to-point mapping — O(s) instead of O(n*s).
  std::vector<int> seed_idx = map_seeds_to_points_grid(hx, hy, seed_x, seed_y);

  std::vector<int> label(static_cast<size_t>(nh), 0);
  std::vector<double> dist(static_cast<size_t>(nh), std::numeric_limits<double>::infinity());

  std::vector<double> center_x(seed_x.begin(), seed_x.end());
  std::vector<double> center_y(seed_y.begin(), seed_y.end());
  std::vector<double> seed_h(seed_z.begin(), seed_z.end());
  std::vector<int> assigned_count(static_cast<size_t>(s), 0);

  typedef std::tuple<double, int, int> Node; // cost, point, label(1-based)
  struct Cmp {
    bool operator()(const Node& a, const Node& b) const {
      if (std::get<0>(a) == std::get<0>(b)) {
        if (std::get<2>(a) == std::get<2>(b)) return std::get<1>(a) > std::get<1>(b);
        return std::get<2>(a) > std::get<2>(b);
      }
      return std::get<0>(a) > std::get<0>(b);
    }
  };

  std::priority_queue<Node, std::vector<Node>, Cmp> pq;

  for (int i = 0; i < s; ++i) {
    int idx = seed_idx[static_cast<size_t>(i)];
    if (idx < 0 || idx >= nh) continue;
    int lid = i + 1;

    if (dist[static_cast<size_t>(idx)] > 0.0) {
      dist[static_cast<size_t>(idx)] = 0.0;
      label[static_cast<size_t>(idx)] = lid;
      assigned_count[static_cast<size_t>(i)] = 1;
      pq.push(std::make_tuple(0.0, idx, lid));
    }
  }

  while (!pq.empty()) {
    Node node = pq.top();
    pq.pop();

    double cur_cost = std::get<0>(node);
    int u = std::get<1>(node);
    int lid = std::get<2>(node);
    int sid = lid - 1;

    if (sid < 0 || sid >= s) continue;
    if (label[static_cast<size_t>(u)] != lid) continue;
    if (cur_cost > dist[static_cast<size_t>(u)]) continue;

    const std::vector<int>& nu = neigh[static_cast<size_t>(u)];
    for (size_t it = 0; it < nu.size(); ++it) {
      int v = nu[it];
      // All hx/hy/hz points are already >= hmin; only skip excessive height jumps.
      if (hz[v] > hz[u] + 1.5 * bin_height) continue;

      double H = seed_h[static_cast<size_t>(sid)];
      if (!std::isfinite(H) || H <= hmin) H = std::max(seed_z[static_cast<size_t>(sid)], hz[u]);
      if (H <= 0) continue;

      double r_allowed = allometry_radius(hz[v], H, crown_a, crown_b, max_crown_radius, profile_id);
      double base_radius = std::max(0.0, crown_a * std::pow(H, crown_b));
      // Avoid over-constraining growth in broad crowns when profile taper is sharp.
      r_allowed = std::max(r_allowed, 0.65 * base_radius);
      double d_center = std::sqrt(sqr(hx[v] - center_x[static_cast<size_t>(sid)]) +
                                  sqr(hy[v] - center_y[static_cast<size_t>(sid)]));

      if (d_center > r_allowed + center_drift) continue;

      double dxy = std::sqrt(sqr(hx[v] - hx[u]) + sqr(hy[v] - hy[u]));
      if (dxy > max_jump * 2.5) continue;

      double dens_pen = 0.0;
      if (density[static_cast<size_t>(v)] < 2) {
        dens_pen = (2.0 - static_cast<double>(density[static_cast<size_t>(v)])) * gap_penalty * 0.5;
      }

      double allom_pen = 0.0;
      if (d_center > r_allowed) {
        allom_pen = (d_center - r_allowed) * 0.8;
      }

      double cc_pen = 0.0;
      if (use_connected_components &&
          comp[static_cast<size_t>(u)] >= 0 &&
          comp[static_cast<size_t>(v)] >= 0 &&
          comp[static_cast<size_t>(u)] != comp[static_cast<size_t>(v)]) {
        cc_pen = gap_penalty * 0.5;
      }

      double edge_cost = dxy + 0.35 * std::abs(hz[u] - hz[v]) + dens_pen + allom_pen + cc_pen;
      double new_cost = cur_cost + edge_cost;

      // The threshold controls local crown continuity, not full seed-to-point
      // path length. Accumulated cost is still used to break ties among seeds.
      if (edge_cost <= geodesic_threshold && new_cost < dist[static_cast<size_t>(v)]) {
        int previous_label = label[static_cast<size_t>(v)];
        dist[static_cast<size_t>(v)] = new_cost;
        label[static_cast<size_t>(v)] = lid;
        pq.push(std::make_tuple(new_cost, v, lid));

        if (previous_label != lid) {
          int cnt = ++assigned_count[static_cast<size_t>(sid)];
          double alpha = 1.0 / static_cast<double>(std::max(2, cnt));
          double nx = center_x[static_cast<size_t>(sid)] + alpha * (hx[v] - center_x[static_cast<size_t>(sid)]);
          double ny = center_y[static_cast<size_t>(sid)] + alpha * (hy[v] - center_y[static_cast<size_t>(sid)]);

          double shift = std::sqrt(sqr(nx - center_x[static_cast<size_t>(sid)]) +
                                   sqr(ny - center_y[static_cast<size_t>(sid)]));
          double max_shift = center_drift;
          if (shift > max_shift && shift > 0) {
            double r = max_shift / shift;
            nx = center_x[static_cast<size_t>(sid)] + (nx - center_x[static_cast<size_t>(sid)]) * r;
            ny = center_y[static_cast<size_t>(sid)] + (ny - center_y[static_cast<size_t>(sid)]) * r;
          }

          center_x[static_cast<size_t>(sid)] = nx;
          center_y[static_cast<size_t>(sid)] = ny;
        }
      }
    }
  }

  // Remove tiny disconnected fragments within each tree label (subset space).
  std::vector<int> visited(static_cast<size_t>(nh), 0);
  for (int i = 0; i < nh; ++i) {
    if (label[static_cast<size_t>(i)] <= 0 || visited[static_cast<size_t>(i)] == 1) continue;

    int lid = label[static_cast<size_t>(i)];
    std::vector<int> stack(1, i);
    std::vector<int> comp_points;
    visited[static_cast<size_t>(i)] = 1;

    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();
      comp_points.push_back(u);

      const std::vector<int>& nu = neigh[static_cast<size_t>(u)];
      for (size_t t = 0; t < nu.size(); ++t) {
        int v = nu[t];
        if (visited[static_cast<size_t>(v)] == 1) continue;
        if (label[static_cast<size_t>(v)] != lid) continue;
        visited[static_cast<size_t>(v)] = 1;
        stack.push_back(v);
      }
    }

    if (static_cast<int>(comp_points.size()) < min_component_size) {
      for (size_t t = 0; t < comp_points.size(); ++t) {
        label[static_cast<size_t>(comp_points[t])] = 0;
      }
    }
  }

  // Fast local fill: one lightweight pass over kNN neighbors (subset space).
  for (int i = 0; i < nh; ++i) {
    if (label[static_cast<size_t>(i)] > 0) continue;

    int best_label = 0;
    double best_cost = std::numeric_limits<double>::infinity();
    const std::vector<int>& ni = neigh[static_cast<size_t>(i)];

    for (size_t t = 0; t < ni.size(); ++t) {
      int j = ni[t];
      int lj = label[static_cast<size_t>(j)];
      if (lj <= 0) continue;

      double dxy = std::sqrt(sqr(hx[i] - hx[j]) + sqr(hy[i] - hy[j]));
      double dz = std::abs(hz[i] - hz[j]);
      double cost = dxy + 0.30 * dz;
      if (cost < best_cost) {
        best_cost = cost;
        best_label = lj;
      }
    }

    if (best_label > 0 && best_cost <= geodesic_threshold * 1.5) {
      label[static_cast<size_t>(i)] = best_label;
    }
  }

  // Map subset labels back to the full point cloud before the broader fill.
  std::vector<int> full_label(static_cast<size_t>(n), 0);
  for (int i = 0; i < nh; ++i) {
    full_label[static_cast<size_t>(hi_idx[static_cast<size_t>(i)])] = label[static_cast<size_t>(i)];
  }

  // Broader fill for remaining canopy holes, but only against currently
  // labeled canopy points to keep runtime bounded.
  std::vector<double> lx;
  std::vector<double> ly;
  std::vector<int> llabel;
  lx.reserve(static_cast<size_t>(nh));
  ly.reserve(static_cast<size_t>(nh));
  llabel.reserve(static_cast<size_t>(nh));

  for (int i = 0; i < nh; ++i) {
    if (label[static_cast<size_t>(i)] > 0) {
      lx.push_back(hx[i]);
      ly.push_back(hy[i]);
      llabel.push_back(label[static_cast<size_t>(i)]);
    }
  }

  if (!lx.empty()) {
    NumericVector lxn(lx.begin(), lx.end());
    NumericVector lyn(ly.begin(), ly.end());
    Grid2D labeled_grid(lxn, lyn, std::max(1.0, max_jump * 1.5));

    const double r0 = std::max(2.0, max_jump * 2.0);
    const double r1 = std::max(3.5, max_jump * 3.5);
    const double fill_limit = std::max(geodesic_threshold * 2.0, max_jump * 3.5);

    for (int i = 0; i < n; ++i) {
      if (z[i] < hmin || full_label[static_cast<size_t>(i)] > 0) continue;

      int best_label = 0;
      double best_cost = std::numeric_limits<double>::infinity();

      std::vector<int> cand = labeled_grid.query_radius(lxn, lyn, x[i], y[i], r0);
      if (cand.empty()) {
        cand = labeled_grid.query_radius(lxn, lyn, x[i], y[i], r1);
      }

      for (size_t t = 0; t < cand.size(); ++t) {
        int j = cand[t];
        double dxy = std::sqrt(sqr(lxn[j] - x[i]) + sqr(lyn[j] - y[i]));
        double cost = dxy;
        if (cost < best_cost) {
          best_cost = cost;
          best_label = llabel[static_cast<size_t>(j)];
        }
      }

      if (best_label > 0 && best_cost <= fill_limit) {
        full_label[static_cast<size_t>(i)] = best_label;
      }
    }
  }

  return labels_to_integer_ids(full_label, z, hmin);
}

extern "C" SEXP cpp_allometric_random_walker(SEXP xSEXP,
                                               SEXP ySEXP,
                                               SEXP zSEXP,
                                               SEXP seedXSEXP,
                                               SEXP seedYSEXP,
                                               SEXP seedZSEXP,
                                               SEXP hminSEXP,
                                               SEXP kSEXP,
                                               SEXP alphaSEXP,
                                               SEXP betaSEXP,
                                               SEXP gammaSEXP,
                                               SEXP deltaSEXP,
                                               SEXP etaSEXP,
                                               SEXP crownASEXP,
                                               SEXP crownBSEXP,
                                               SEXP eigenRadiusSEXP,
                                               SEXP densityRadiusSEXP,
                                               SEXP probThresholdSEXP,
                                               SEXP maxIterationsSEXP,
                                               SEXP toleranceSEXP,
                                               SEXP profileIdSEXP) {
  NumericVector x(xSEXP), y(ySEXP), z(zSEXP);
  NumericVector seed_x(seedXSEXP), seed_y(seedYSEXP), seed_z(seedZSEXP);

  const int n = x.size();
  if (y.size() != n || z.size() != n) {
    stop("X, Y, Z must share length.");
  }
  if (n == 0) return IntegerVector();

  const double hmin = as<double>(hminSEXP);
  const int k = std::max(2, as<int>(kSEXP));
  const double alpha = std::max(0.0, as<double>(alphaSEXP));
  const double beta = std::max(0.0, as<double>(betaSEXP));
  const double gamma = std::max(0.0, as<double>(gammaSEXP));
  const double delta = std::max(0.0, as<double>(deltaSEXP));
  const double eta = std::max(0.0, as<double>(etaSEXP));
  const double crown_a = as<double>(crownASEXP);
  const double crown_b = as<double>(crownBSEXP);
  const double eigen_radius = std::max(1e-6, as<double>(eigenRadiusSEXP));
  const double density_radius = std::max(1e-6, as<double>(densityRadiusSEXP));
  const double prob_threshold = std::min(1.0, std::max(0.0, as<double>(probThresholdSEXP)));
  const int max_iterations = std::max(1, as<int>(maxIterationsSEXP));
  const double tolerance = std::max(0.0, as<double>(toleranceSEXP));
  const int profile_id = as<int>(profileIdSEXP);

  const int s = seed_x.size();
  if (seed_y.size() != s || seed_z.size() != s || s == 0) {
    return IntegerVector(n, NA_INTEGER);
  }

  double base_radius = std::max(eigen_radius, density_radius);
  std::vector< std::vector<int> > neigh = build_knn(x, y, k, base_radius, base_radius * 2.5, 5);
  std::vector<int> density = local_density(x, y, density_radius, density_radius);

  std::vector<double> anisotropy, verticality;
  compute_local_eigenfeatures(x, y, z, eigen_radius, anisotropy, verticality);

  std::vector<int> seed_idx = map_seeds_to_points(x, y, z, seed_x, seed_y, hmin);

  std::vector<int> labels(static_cast<size_t>(n), 0);
  std::vector<double> confidence(static_cast<size_t>(n), 0.0);
  std::vector<int> is_seed(static_cast<size_t>(n), 0);
  std::vector<double> seed_height(static_cast<size_t>(s), hmin + 0.1);

  for (int sid = 0; sid < s; ++sid) {
    double H = seed_z[static_cast<size_t>(sid)];
    if (!std::isfinite(H) || H <= hmin) H = hmin + 0.1;
    seed_height[static_cast<size_t>(sid)] = H;
  }

  for (int i = 0; i < s; ++i) {
    int idx = seed_idx[static_cast<size_t>(i)];
    if (idx < 0 || idx >= n) continue;
    labels[static_cast<size_t>(idx)] = i + 1;
    confidence[static_cast<size_t>(idx)] = 1.0;
    is_seed[static_cast<size_t>(idx)] = 1;
  }

  std::vector<int> next_labels = labels;
  std::vector<double> next_conf = confidence;

  for (int iter = 0; iter < max_iterations; ++iter) {
    int changed = 0;

    for (int i = 0; i < n; ++i) {
      if (z[i] < hmin) {
        next_labels[static_cast<size_t>(i)] = 0;
        next_conf[static_cast<size_t>(i)] = 0.0;
        continue;
      }
      if (is_seed[static_cast<size_t>(i)] == 1) continue;

      // Only labels present in local neighbors can contribute, so keep this sparse.
      std::vector<int> cand_labels;
      std::vector<double> cand_scores;
      cand_labels.reserve(16);
      cand_scores.reserve(16);
      double total = 0.0;

      const std::vector<int>& ni = neigh[static_cast<size_t>(i)];
      for (size_t t = 0; t < ni.size(); ++t) {
        int j = ni[t];
        int lj = labels[static_cast<size_t>(j)];
        if (lj <= 0) continue;

        double dij = std::sqrt(sqr(x[i] - x[j]) + sqr(y[i] - y[j]));
        double dz = std::abs(z[i] - z[j]);
        double deigen = std::abs(anisotropy[static_cast<size_t>(i)] - anisotropy[static_cast<size_t>(j)]) +
                        std::abs(verticality[static_cast<size_t>(i)] - verticality[static_cast<size_t>(j)]);

        int dmin = std::min(density[static_cast<size_t>(i)], density[static_cast<size_t>(j)]);
        double dens_gap = std::max(0.0, 2.0 - static_cast<double>(dmin));

        double expo = alpha * dij + beta * dz + gamma * deigen + eta * dens_gap;
        double w = std::exp(-expo);

        bool found = false;
        for (size_t k2 = 0; k2 < cand_labels.size(); ++k2) {
          if (cand_labels[k2] == lj) {
            cand_scores[k2] += w;
            found = true;
            break;
          }
        }
        if (!found) {
          cand_labels.push_back(lj);
          cand_scores.push_back(w);
        }

        total += w;
      }

      int best_label = 0;
      double best_score = 0.0;

      if (total > 0.0 && !cand_labels.empty()) {
        // Apply allometry penalty once per candidate label instead of once per neighbor.
        for (size_t c = 0; c < cand_labels.size(); ++c) {
          int lj = cand_labels[c];
          int sid = lj - 1;

          double dseed = std::sqrt(sqr(x[i] - seed_x[static_cast<size_t>(sid)]) +
                                   sqr(y[i] - seed_y[static_cast<size_t>(sid)]));
          double rall = allometry_radius(z[i], seed_height[static_cast<size_t>(sid)],
                                         crown_a, crown_b,
                                         std::numeric_limits<double>::infinity(),
                                         profile_id);
          double allom_pen = dseed > rall ? (dseed - rall) : 0.0;

          double sc = cand_scores[c] * std::exp(-delta * allom_pen);
          if (sc > best_score) {
            best_score = sc;
            best_label = lj;
          }
        }
      }

      double p = total > 0.0 ? (best_score / total) : 0.0;
      if (p < prob_threshold) best_label = 0;

      if (best_label != labels[static_cast<size_t>(i)] ||
          std::abs(p - confidence[static_cast<size_t>(i)]) > tolerance) {
        ++changed;
      }

      next_labels[static_cast<size_t>(i)] = best_label;
      next_conf[static_cast<size_t>(i)] = p;
    }

    labels.swap(next_labels);
    confidence.swap(next_conf);

    if (changed <= static_cast<int>(tolerance * std::max(1, n))) {
      break;
    }
  }

  return labels_to_integer_ids(labels, z, hmin);
}

extern "C" SEXP cpp_allometric_supervoxel_segment(SEXP xSEXP,
                                                    SEXP ySEXP,
                                                    SEXP zSEXP,
                                                    SEXP seedXSEXP,
                                                    SEXP seedYSEXP,
                                                    SEXP seedZSEXP,
                                                    SEXP hminSEXP,
                                                    SEXP voxelSizeSEXP,
                                                    SEXP kVoxelSEXP,
                                                    SEXP minVoxelPointsSEXP,
                                                    SEXP crownASEXP,
                                                    SEXP crownBSEXP,
                                                    SEXP componentGapSEXP,
                                                    SEXP mergeThresholdSEXP,
                                                    SEXP anisotropyWeightSEXP,
                                                    SEXP verticalityWeightSEXP,
                                                    SEXP densityWeightSEXP,
                                                    SEXP allometryWeightSEXP,
                                                    SEXP profileIdSEXP) {
  NumericVector x(xSEXP), y(ySEXP), z(zSEXP);
  NumericVector seed_x(seedXSEXP), seed_y(seedYSEXP), seed_z(seedZSEXP);

  const int n = x.size();
  if (y.size() != n || z.size() != n) {
    stop("X, Y, Z must share length.");
  }
  if (n == 0) return IntegerVector();

  const double hmin = as<double>(hminSEXP);
  const double voxel_size = std::max(1e-6, as<double>(voxelSizeSEXP));
  const int k_voxel = std::max(1, as<int>(kVoxelSEXP));
  const int min_voxel_points = std::max(1, as<int>(minVoxelPointsSEXP));
  const double crown_a = as<double>(crownASEXP);
  const double crown_b = as<double>(crownBSEXP);
  const double component_gap = std::max(1e-6, as<double>(componentGapSEXP));
  const double merge_threshold = std::max(1e-6, as<double>(mergeThresholdSEXP));
  const double anis_w = std::max(0.0, as<double>(anisotropyWeightSEXP));
  const double vert_w = std::max(0.0, as<double>(verticalityWeightSEXP));
  const double dens_w = std::max(0.0, as<double>(densityWeightSEXP));
  const double allom_w = std::max(0.0, as<double>(allometryWeightSEXP));
  const int profile_id = as<int>(profileIdSEXP);

  const int s = seed_x.size();
  if (seed_y.size() != s || seed_z.size() != s || s == 0) {
    return IntegerVector(n, NA_INTEGER);
  }

  struct VoxelAgg {
    int ix;
    int iy;
    int iz;
    std::vector<int> points;
    double sumx;
    double sumy;
    double sumz;
    double minz;
    double maxz;
    double anis;
    double vert;
    double dens;
    bool active;
  };

  std::unordered_map<long long, int> voxel_lookup;
  std::vector<VoxelAgg> voxels;
  voxels.reserve(static_cast<size_t>(n / 3 + 1));

  std::vector<int> point_voxel(static_cast<size_t>(n), -1);

  for (int i = 0; i < n; ++i) {
    if (z[i] < hmin) continue;

    int ix = static_cast<int>(std::floor(x[i] / voxel_size));
    int iy = static_cast<int>(std::floor(y[i] / voxel_size));
    int iz = static_cast<int>(std::floor(z[i] / voxel_size));

    long long key = key3d(ix, iy, iz);
    auto it = voxel_lookup.find(key);

    int vid;
    if (it == voxel_lookup.end()) {
      vid = static_cast<int>(voxels.size());
      voxel_lookup[key] = vid;

      VoxelAgg v;
      v.ix = ix;
      v.iy = iy;
      v.iz = iz;
      v.sumx = x[i];
      v.sumy = y[i];
      v.sumz = z[i];
      v.minz = z[i];
      v.maxz = z[i];
      v.anis = 0.5;
      v.vert = 0.7;
      v.dens = 0.0;
      v.active = true;
      v.points.push_back(i);
      voxels.push_back(v);
    } else {
      vid = it->second;
      VoxelAgg& v = voxels[static_cast<size_t>(vid)];
      v.sumx += x[i];
      v.sumy += y[i];
      v.sumz += z[i];
      v.minz = std::min(v.minz, static_cast<double>(z[i]));
      v.maxz = std::max(v.maxz, static_cast<double>(z[i]));
      v.points.push_back(i);
    }

    point_voxel[static_cast<size_t>(i)] = vid;
  }

  const int m = static_cast<int>(voxels.size());
  if (m == 0) {
    return IntegerVector(n, NA_INTEGER);
  }

  NumericVector vx(m), vy(m), vz(m);
  std::vector<int> count(static_cast<size_t>(m), 0);

  for (int i = 0; i < m; ++i) {
    VoxelAgg& v = voxels[static_cast<size_t>(i)];
    int c = static_cast<int>(v.points.size());
    count[static_cast<size_t>(i)] = c;
    vx[i] = v.sumx / static_cast<double>(c);
    vy[i] = v.sumy / static_cast<double>(c);
    vz[i] = v.sumz / static_cast<double>(c);
    v.dens = static_cast<double>(c) / std::max(1e-6, voxel_size * voxel_size * voxel_size);
    v.active = false;
  }

  int effective_min_voxel_points = min_voxel_points;
  while (effective_min_voxel_points > 1) {
    int active_count = 0;
    for (int i = 0; i < m; ++i) {
      if (count[static_cast<size_t>(i)] >= effective_min_voxel_points) {
        ++active_count;
      }
    }

    double active_fraction = static_cast<double>(active_count) / static_cast<double>(m);
    if (active_fraction >= 0.50) {
      break;
    }

    --effective_min_voxel_points;
  }

  for (int i = 0; i < m; ++i) {
    voxels[static_cast<size_t>(i)].active = count[static_cast<size_t>(i)] >= effective_min_voxel_points;
  }

  // Estimate voxel eigenfeatures from points in each voxel.
  for (int i = 0; i < m; ++i) {
    VoxelAgg& v = voxels[static_cast<size_t>(i)];
    if (v.points.size() < 4) continue;

    arma::mat cov(3, 3, arma::fill::zeros);
    for (size_t t = 0; t < v.points.size(); ++t) {
      int p = v.points[t];
      arma::vec d(3);
      d(0) = x[p] - vx[i];
      d(1) = y[p] - vy[i];
      d(2) = z[p] - vz[i];
      cov += d * d.t();
    }
    cov /= std::max(1.0, static_cast<double>(v.points.size()) - 1.0);

    arma::vec eigval;
    arma::mat eigvec;
    if (!arma::eig_sym(eigval, eigvec, cov) || eigval.n_elem != 3) continue;

    double l1 = std::max(1e-12, static_cast<double>(eigval(2)));
    double l3 = std::max(0.0, static_cast<double>(eigval(0)));
    v.anis = clamp01((l1 - l3) / l1);
    arma::vec principal = eigvec.col(2);
    v.vert = clamp01(std::abs(principal(2)));
  }

  // Build a centroid kNN graph among voxels. This is more robust than strict
  // adjacency in sparsely occupied voxel grids produced by ALS subsampling.
  std::vector< std::vector<int> > neigh = build_knn(
    vx,
    vy,
    std::max(2, k_voxel),
    voxel_size * 2.0,
    std::max(component_gap, voxel_size * 3.0),
    6
  );

  // Connected components among active voxels.
  std::vector<int> comp(static_cast<size_t>(m), -1);
  int cid = 0;
  for (int i = 0; i < m; ++i) {
    if (!voxels[static_cast<size_t>(i)].active || comp[static_cast<size_t>(i)] >= 0) continue;

    std::vector<int> stack(1, i);
    comp[static_cast<size_t>(i)] = cid;

    while (!stack.empty()) {
      int u = stack.back();
      stack.pop_back();

      const std::vector<int>& nu = neigh[static_cast<size_t>(u)];
      for (size_t t = 0; t < nu.size(); ++t) {
        int v = nu[t];
        if (!voxels[static_cast<size_t>(v)].active || comp[static_cast<size_t>(v)] >= 0) continue;

        double d = std::sqrt(sqr(vx[u] - vx[v]) + sqr(vy[u] - vy[v]) + sqr(vz[u] - vz[v]));
        if (d > component_gap) continue;

        comp[static_cast<size_t>(v)] = cid;
        stack.push_back(v);
      }
    }

    ++cid;
  }

  // Map seeds to nearest active voxel using XY + height proximity. XY-only
  // snapping can anchor seeds to neighboring crowns in dense canopies.
  std::vector<int> seed_voxel(static_cast<size_t>(s), -1);
  for (int i = 0; i < s; ++i) {
    double best_local = std::numeric_limits<double>::infinity();
    double best_global = std::numeric_limits<double>::infinity();
    const double local_xy = std::max(component_gap * 1.5, voxel_size * 3.0);
    const double local_xy2 = local_xy * local_xy;
    const double seed_h = std::isfinite(seed_z[i]) ? seed_z[i] : hmin;
    int best_j = -1;

    for (int j = 0; j < m; ++j) {
      if (!voxels[static_cast<size_t>(j)].active) continue;

      double dxy2 = sqr(vx[j] - seed_x[i]) + sqr(vy[j] - seed_y[i]);
      double dxy = std::sqrt(dxy2);
      double dz = std::abs(vz[j] - seed_h);

      // Prefer nearby voxels, but account for canopy height consistency.
      double local_score = dxy + 0.60 * dz;
      double global_score = dxy + 0.85 * dz;

      if (dxy2 <= local_xy2 && local_score < best_local) {
        best_local = local_score;
        best_j = j;
      }

      if (best_j < 0 && global_score < best_global) {
        best_global = global_score;
        best_j = j;
      }
    }

    seed_voxel[static_cast<size_t>(i)] = best_j;
  }

  // Multi-source Dijkstra on voxel graph.
  std::vector<int> label(static_cast<size_t>(m), 0);
  std::vector<double> dist(static_cast<size_t>(m), std::numeric_limits<double>::infinity());

  typedef std::tuple<double, int, int> Node; // cost, voxel, label
  struct Cmp {
    bool operator()(const Node& a, const Node& b) const {
      if (std::get<0>(a) == std::get<0>(b)) {
        if (std::get<2>(a) == std::get<2>(b)) return std::get<1>(a) > std::get<1>(b);
        return std::get<2>(a) > std::get<2>(b);
      }
      return std::get<0>(a) > std::get<0>(b);
    }
  };

  std::priority_queue<Node, std::vector<Node>, Cmp> pq;

  for (int i = 0; i < s; ++i) {
    int v = seed_voxel[static_cast<size_t>(i)];
    if (v < 0 || v >= m) continue;
    if (!voxels[static_cast<size_t>(v)].active) continue;

    int lid = i + 1;
    if (dist[static_cast<size_t>(v)] > 0.0) {
      dist[static_cast<size_t>(v)] = 0.0;
      label[static_cast<size_t>(v)] = lid;
      pq.push(std::make_tuple(0.0, v, lid));
    }
  }

  double mean_density = 0.0;
  int active_n = 0;
  for (int i = 0; i < m; ++i) {
    if (!voxels[static_cast<size_t>(i)].active) continue;
    mean_density += voxels[static_cast<size_t>(i)].dens;
    ++active_n;
  }
  mean_density = active_n > 0 ? mean_density / static_cast<double>(active_n) : 1.0;

  while (!pq.empty()) {
    Node node = pq.top();
    pq.pop();

    double cur = std::get<0>(node);
    int u = std::get<1>(node);
    int lid = std::get<2>(node);
    int sid = lid - 1;

    if (sid < 0 || sid >= s) continue;
    if (label[static_cast<size_t>(u)] != lid) continue;
    if (cur > dist[static_cast<size_t>(u)]) continue;

    const std::vector<int>& nu = neigh[static_cast<size_t>(u)];
    for (size_t t = 0; t < nu.size(); ++t) {
      int v = nu[t];
      if (!voxels[static_cast<size_t>(v)].active) continue;

      double d = std::sqrt(sqr(vx[u] - vx[v]) + sqr(vy[u] - vy[v]) + sqr(vz[u] - vz[v]));
      double anis_pen = anis_w * std::abs(voxels[static_cast<size_t>(u)].anis - voxels[static_cast<size_t>(v)].anis);
      double vert_pen = vert_w * std::abs(voxels[static_cast<size_t>(u)].vert - voxels[static_cast<size_t>(v)].vert);

      double dens_ratio = std::min(voxels[static_cast<size_t>(u)].dens, voxels[static_cast<size_t>(v)].dens) /
                          std::max(1e-9, mean_density);
      double dens_pen = dens_w * std::max(0.0, 1.0 - dens_ratio);

      double H = seed_z[static_cast<size_t>(sid)];
      if (!std::isfinite(H) || H <= hmin) H = std::max(hmin + 0.1, vz[v]);

      double dseed = std::sqrt(sqr(vx[v] - seed_x[static_cast<size_t>(sid)]) +
               sqr(vy[v] - seed_y[static_cast<size_t>(sid)]));
      double rall = allometry_radius(vz[v], H, crown_a, crown_b,
                                     std::numeric_limits<double>::infinity(),
                                     profile_id);
      double allom_pen = allom_w * std::max(0.0, dseed - rall);

      double comp_pen = 0.0;
      if (comp[static_cast<size_t>(u)] >= 0 &&
          comp[static_cast<size_t>(v)] >= 0 &&
          comp[static_cast<size_t>(u)] != comp[static_cast<size_t>(v)]) {
        comp_pen = component_gap;
      }

      double local_edge = d + anis_pen + vert_pen + dens_pen + comp_pen;
      double edge = local_edge + allom_pen;
      double nd = cur + edge;

      // Gate only on local structural coherence. Keep allometry as a soft
      // cumulative penalty so crowns can still expand in sparse canopies.
      if (local_edge <= merge_threshold && nd < dist[static_cast<size_t>(v)]) {
        dist[static_cast<size_t>(v)] = nd;
        label[static_cast<size_t>(v)] = lid;
        pq.push(std::make_tuple(nd, v, lid));
      }
    }
  }

  // Fill unlabeled voxels by borrowing the nearest plausible labeled voxel.
  // This reduces isolated NA speckles inside crowns caused by hard graph gating.
  std::vector<double> lvx;
  std::vector<double> lvy;
  std::vector<double> lvz;
  std::vector<int> llab;
  lvx.reserve(static_cast<size_t>(m));
  lvy.reserve(static_cast<size_t>(m));
  lvz.reserve(static_cast<size_t>(m));
  llab.reserve(static_cast<size_t>(m));

  for (int i = 0; i < m; ++i) {
    if (label[static_cast<size_t>(i)] > 0) {
      lvx.push_back(vx[i]);
      lvy.push_back(vy[i]);
      lvz.push_back(vz[i]);
      llab.push_back(label[static_cast<size_t>(i)]);
    }
  }

  if (!lvx.empty()) {
    NumericVector lvx_n(lvx.begin(), lvx.end());
    NumericVector lvy_n(lvy.begin(), lvy.end());
    Grid2D labeled_grid(lvx_n, lvy_n, std::max(voxel_size * 1.5, 1e-6));

    const double r0 = std::max(component_gap * 1.5, voxel_size * 3.0);
    const double r1 = std::max(component_gap * 2.5, voxel_size * 5.0);

    for (int i = 0; i < m; ++i) {
      if (label[static_cast<size_t>(i)] > 0) continue;

      std::vector<int> cand = labeled_grid.query_radius(lvx_n, lvy_n, vx[i], vy[i], r0);
      if (cand.empty()) {
        cand = labeled_grid.query_radius(lvx_n, lvy_n, vx[i], vy[i], r1);
      }
      if (cand.empty()) continue;

      int best_label = 0;
      double best_cost = std::numeric_limits<double>::infinity();

      for (size_t t = 0; t < cand.size(); ++t) {
        int j = cand[t];
        int lid = llab[static_cast<size_t>(j)];
        int sid = lid - 1;
        if (sid < 0 || sid >= s) continue;

        double dxy = std::sqrt(sqr(vx[i] - lvx_n[j]) + sqr(vy[i] - lvy_n[j]));
        double dz = std::abs(vz[i] - lvz[static_cast<size_t>(j)]);
        double cost = dxy + 0.30 * dz;

        double H = seed_z[static_cast<size_t>(sid)];
        if (!std::isfinite(H) || H <= hmin) H = std::max(hmin + 0.1, vz[i]);

        double dseed = std::sqrt(sqr(vx[i] - seed_x[static_cast<size_t>(sid)]) +
                                 sqr(vy[i] - seed_y[static_cast<size_t>(sid)]));
        double rall = allometry_radius(vz[i], H, crown_a, crown_b,
                                       std::numeric_limits<double>::infinity(),
                                       profile_id);
        double relaxed_limit = std::max(component_gap * 1.5, 2.3 * rall);
        if (dseed > relaxed_limit) continue;

        if (cost < best_cost) {
          best_cost = cost;
          best_label = lid;
        }
      }

      if (best_label > 0) {
        label[static_cast<size_t>(i)] = best_label;
      }
    }
  }

  // Remove tiny disconnected label islands that appear as implausible speckles,
  // then reattach them to nearby dominant labels.
  {
    std::vector<int> visited(static_cast<size_t>(m), 0);
    std::unordered_map<int, int> label_max_points;

    struct CompInfo {
      int lbl;
      std::vector<int> vox;
      int points;
    };
    std::vector<CompInfo> comps;
    comps.reserve(static_cast<size_t>(m / 4 + 1));

    for (int i = 0; i < m; ++i) {
      int lbl = label[static_cast<size_t>(i)];
      if (lbl <= 0 || visited[static_cast<size_t>(i)] == 1) continue;

      std::vector<int> stack(1, i);
      visited[static_cast<size_t>(i)] = 1;

      CompInfo ci;
      ci.lbl = lbl;
      ci.points = 0;

      while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        ci.vox.push_back(u);
        ci.points += count[static_cast<size_t>(u)];

        const std::vector<int>& nu = neigh[static_cast<size_t>(u)];
        for (size_t t = 0; t < nu.size(); ++t) {
          int v = nu[t];
          if (visited[static_cast<size_t>(v)] == 1) continue;
          if (label[static_cast<size_t>(v)] != lbl) continue;
          visited[static_cast<size_t>(v)] = 1;
          stack.push_back(v);
        }
      }

      auto it = label_max_points.find(lbl);
      if (it == label_max_points.end()) {
        label_max_points[lbl] = ci.points;
      } else {
        it->second = std::max(it->second, ci.points);
      }

      comps.push_back(ci);
    }

    const int min_pts_keep = std::max(6, effective_min_voxel_points * 3);
    for (size_t cidx = 0; cidx < comps.size(); ++cidx) {
      const CompInfo& ci = comps[cidx];
      int max_pts = label_max_points[ci.lbl];
      bool keep = (ci.points >= min_pts_keep) ||
                  (max_pts > 0 && ci.points >= static_cast<int>(0.20 * static_cast<double>(max_pts)));

      if (!keep) {
        for (size_t t = 0; t < ci.vox.size(); ++t) {
          label[static_cast<size_t>(ci.vox[t])] = 0;
        }
      }
    }

    // Reattach removed islands to nearby labels in the local voxel graph.
    for (int pass = 0; pass < 2; ++pass) {
      bool changed = false;

      for (int i = 0; i < m; ++i) {
        if (label[static_cast<size_t>(i)] > 0) continue;

        int best_lbl = 0;
        double best_cost = std::numeric_limits<double>::infinity();

        const std::vector<int>& ni = neigh[static_cast<size_t>(i)];
        for (size_t t = 0; t < ni.size(); ++t) {
          int j = ni[t];
          int lj = label[static_cast<size_t>(j)];
          if (lj <= 0) continue;

          double dxy = std::sqrt(sqr(vx[i] - vx[j]) + sqr(vy[i] - vy[j]));
          double dz = std::abs(vz[i] - vz[j]);
          double cost = dxy + 0.25 * dz;

          if (cost < best_cost) {
            best_cost = cost;
            best_lbl = lj;
          }
        }

        if (best_lbl > 0 && best_cost <= std::max(component_gap * 1.5, voxel_size * 3.5)) {
          label[static_cast<size_t>(i)] = best_lbl;
          changed = true;
        }
      }

      if (!changed) break;
    }
  }

  IntegerVector out(n, NA_INTEGER);
  for (int i = 0; i < n; ++i) {
    if (z[i] < hmin) {
      out[i] = NA_INTEGER;
      continue;
    }

    int vid = point_voxel[static_cast<size_t>(i)];
    if (vid < 0 || vid >= m) {
      out[i] = NA_INTEGER;
      continue;
    }

    int lid = label[static_cast<size_t>(vid)];
    out[i] = lid > 0 ? lid : NA_INTEGER;
  }

  return out;
}
