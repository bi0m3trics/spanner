//  ============================================================================
//  spanner – bole_classify.cpp
//
//  C_ransac_bole_slices():
//    For every tree in a segmented point cloud, iterate over height slices and
//    fit a RANSAC circle to bole-candidate points within each slice.  Returns a
//    data.frame with one row per (TreeID × slice) suitable for downstream
//    volume estimation and quality assessment.
//
//  Author: Andrew Sánchez Meador (bi0m3trics), 2026
//  License: GPL-3
//  ============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include "bole_classify.h"
#include "omp.h"

#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <limits>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// C_ransac_bole_slices
//
// Parameters
//   X, Y, Z        – point coordinates (same length, from las@data)
//   treeIds        – integer tree ID per point (0 = unassigned / ground)
//   tree_X, tree_Y – XY centroids of each tree (from get_raster_eigen_treelocs)
//   tree_ID_vals   – tree IDs corresponding to tree_X / tree_Y rows
//   z_min, z_max   – height window to search for bole points (above ground)
//   dz             – height of each slice
//   overlap        – fractional overlap between adjacent slices (0–<1)
//   search_radius  – max XY distance from tree axis to include a point
//   inlier_tol     – distance from fitted circle surface counting as inlier
//   n_samples      – RANSAC points sampled per iteration
//   confidence     – RANSAC confidence level (0–1)
//   inlier_frac    – expected inlier fraction (0–1)
//   n_best         – number of best RANSAC solutions to median-aggregate
//   max_radius     – reject any cylinder with radius > this value (m)
//   min_pts        – minimum points in a slice to attempt fitting
//
// Returns  DataFrame with columns:
//   TreeID, Segment, Z_low, Z_high, Z_mid,
//   X_ctr, Y_ctr, Radius, RANSAC_err,
//   N_pts, N_inliers, Inlier_frac, Valid
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
DataFrame C_ransac_bole_slices(
    NumericVector X,
    NumericVector Y,
    NumericVector Z,
    IntegerVector treeIds,
    NumericVector tree_X,
    NumericVector tree_Y,
    IntegerVector tree_ID_vals,
    double z_min,
    double z_max,
    double dz,
    double overlap,
    double search_radius,
    double inlier_tol,
    int    n_samples,
    double confidence,
    double inlier_frac,
    int    n_best,
    double max_radius,
    int    min_pts,
    int    ncpu)
{
    int n_pts   = X.size();
    int n_trees = tree_ID_vals.size();

    // -----------------------------------------------------------------------
    // Build lookup: tree_id → (centroid_X, centroid_Y)
    // -----------------------------------------------------------------------
    std::unordered_map<int, std::pair<double, double>> tree_xy_map;
    tree_xy_map.reserve(n_trees);
    for (int t = 0; t < n_trees; ++t) {
        tree_xy_map[tree_ID_vals[t]] = std::make_pair(
            (double)tree_X[t], (double)tree_Y[t]);
    }

    // -----------------------------------------------------------------------
    // Build lookup: tree_id → indices of points belonging to that tree
    // -----------------------------------------------------------------------
    std::unordered_map<int, std::vector<int>> tree_pt_idx;
    tree_pt_idx.reserve(n_trees);
    for (int i = 0; i < n_pts; ++i) {
        int tid = treeIds[i];
        if (tid > 0 && tree_xy_map.count(tid)) {
            tree_pt_idx[tid].push_back(i);
        }
    }

    // -----------------------------------------------------------------------
    // Slice geometry
    // -----------------------------------------------------------------------
    double step = dz * (1.0 - std::max(0.0, std::min(overlap, 0.99)));
    if (step <= 0.0) step = dz;
    int n_slices = std::max(1, (int)std::ceil((z_max - z_min) / step));

    double sr2 = search_radius * search_radius;

    // -----------------------------------------------------------------------
    // Per-tree result storage (pre-allocated so each OMP thread writes to its
    // own slot — no shared-state races, no critical sections needed).
    // -----------------------------------------------------------------------
    struct SliceRow {
        int    treeid, segment, npts, ninliers;
        double zlow, zhigh, zmid;
        double xctr, yctr, radius, err, infrac;
        bool   valid;
    };
    // tree_rows[t] holds all slice rows for tree index t
    std::vector<std::vector<SliceRow>> tree_rows(n_trees);

    // -----------------------------------------------------------------------
    // OpenMP setup
    //   • Clamp thread count to [1, n_trees] so we never spin idle threads.
    //   • Disable nested parallelism: prevents oversubscription if the caller
    //     (R's parallel backend, OpenBLAS, etc.) already holds a thread team.
    // -----------------------------------------------------------------------
    int ncpu_prev = omp_get_max_threads();
    int ncpu_use  = std::max(1, std::min(ncpu, n_trees));
    omp_set_max_active_levels(1);   // no nested OMP teams
    omp_set_num_threads(ncpu_use);

    // -----------------------------------------------------------------------
    // Main loop: tree → slice → RANSAC  (parallel over trees)
    // -----------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int t = 0; t < n_trees; ++t) {
        int tid = tree_ID_vals[t];

        auto it = tree_pt_idx.find(tid);
        if (it == tree_pt_idx.end()) continue;

        double tx = tree_X[t];
        double ty = tree_Y[t];
        const std::vector<int>& pidx = it->second;

        // Per-thread RNG — seeded from hardware entropy xor'd with tree index
        // so different trees get different streams even with a fixed global seed.
        std::random_device rd;
        std::mt19937 rng(rd() ^ (std::mt19937::result_type)(t * 2654435761u));

        unsigned int ns = (unsigned int)std::max(3, n_samples);
        unsigned int nb = (unsigned int)std::max(0, n_best);

        int seg = 0;
        std::vector<SliceRow>& rows = tree_rows[t];  // thread-private slot

        for (int s = 0; s < n_slices; ++s) {
            double zlo  = z_min + s * step;
            double zhi  = zlo + dz;
            double zmid = 0.5 * (zlo + zhi);

            // Collect bole-candidate points in this slice
            std::vector<std::vector<double>> slice(3);
            for (int i : pidx) {
                double zi = Z[i];
                if (zi < zlo || zi > zhi) continue;
                double dx = X[i] - tx;
                double dy = Y[i] - ty;
                if (dx * dx + dy * dy > sr2) continue;
                slice[0].push_back(X[i]);
                slice[1].push_back(Y[i]);
                slice[2].push_back(zi);
            }

            int npts = (int)slice[0].size();
            ++seg;

            SliceRow row;
            row.treeid  = tid;
            row.segment = seg;
            row.zlow    = zlo;
            row.zhigh   = zhi;
            row.zmid    = zmid;
            row.npts    = npts;

            if (npts < min_pts) {
                row.xctr = row.yctr = row.radius = row.err = row.infrac =
                    std::numeric_limits<double>::quiet_NaN();
                row.ninliers = 0;
                row.valid    = false;
                rows.push_back(row);
                continue;
            }

            // Thread-safe RANSAC using per-thread RNG
            std::vector<double> fit =
                ransacCircleTs(slice, ns, confidence, inlier_frac, nb, rng);

            double cx  = fit[0];
            double cy  = fit[1];
            double r   = fit[2];
            double err = fit[3];

            if (r > max_radius || r <= 0.0 || std::isnan(r) || std::isnan(cx)) {
                row.xctr     = cx;
                row.yctr     = cy;
                row.radius   = r;
                row.err      = err;
                row.ninliers = 0;
                row.infrac   = std::numeric_limits<double>::quiet_NaN();
                row.valid    = false;
                rows.push_back(row);
                continue;
            }

            // Count inliers
            int n_inliers = 0;
            for (int i = 0; i < npts; ++i) {
                double dx = slice[0][i] - cx;
                double dy = slice[1][i] - cy;
                if (std::fabs(std::sqrt(dx*dx + dy*dy) - r) <= inlier_tol)
                    ++n_inliers;
            }

            row.xctr     = cx;
            row.yctr     = cy;
            row.radius   = r;
            row.err      = err;
            row.ninliers = n_inliers;
            row.infrac   = (double)n_inliers / (double)npts;
            row.valid    = true;
            rows.push_back(row);

        } // end slice loop
    } // end tree loop

    // Restore caller's thread count
    omp_set_num_threads(ncpu_prev);

    // -----------------------------------------------------------------------
    // Flatten per-tree results into output vectors (sequential, no races)
    // -----------------------------------------------------------------------
    std::vector<int>    out_treeid, out_segment, out_npts, out_ninliers;
    std::vector<double> out_zlow, out_zhigh, out_zmid;
    std::vector<double> out_xctr, out_yctr, out_radius, out_err, out_infrac;
    std::vector<bool>   out_valid;

    for (int t = 0; t < n_trees; ++t) {
        for (const SliceRow& row : tree_rows[t]) {
            out_treeid.push_back(row.treeid);
            out_segment.push_back(row.segment);
            out_zlow.push_back(row.zlow);
            out_zhigh.push_back(row.zhigh);
            out_zmid.push_back(row.zmid);
            out_npts.push_back(row.npts);
            out_xctr.push_back(std::isnan(row.xctr) ? NA_REAL : row.xctr);
            out_yctr.push_back(std::isnan(row.yctr) ? NA_REAL : row.yctr);
            out_radius.push_back(std::isnan(row.radius) ? NA_REAL : row.radius);
            out_err.push_back(std::isnan(row.err) ? NA_REAL : row.err);
            out_ninliers.push_back(row.ninliers);
            out_infrac.push_back(std::isnan(row.infrac) ? NA_REAL : row.infrac);
            out_valid.push_back(row.valid);
        }
    }

    // -----------------------------------------------------------------------
    // Assemble output DataFrame
    // -----------------------------------------------------------------------
    return DataFrame::create(
        Named("TreeID")      = out_treeid,
        Named("Segment")     = out_segment,
        Named("Z_low")       = out_zlow,
        Named("Z_high")      = out_zhigh,
        Named("Z_mid")       = out_zmid,
        Named("X_ctr")       = out_xctr,
        Named("Y_ctr")       = out_yctr,
        Named("Radius")      = out_radius,
        Named("RANSAC_err")  = out_err,
        Named("N_pts")       = out_npts,
        Named("N_inliers")   = out_ninliers,
        Named("Inlier_frac") = out_infrac,
        Named("Valid")       = out_valid,
        Named("stringsAsFactors") = false
    );
}
