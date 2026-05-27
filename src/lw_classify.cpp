// ============================================================================
// spanner – lw_classify.cpp
//
// C_knn_majority_vote()
//   Iterates LeWoS-style majority-vote label propagation over a pre-built kNN
//   index (as returned by FNN::knn.index) entirely in C++, updating labels
//   in-place across iterations.  Returns the converged labels and per-point
//   vote fractions (BoleProb) from the final round — all in a single pass
//   without allocating any intermediate R objects.
//
// Why: the equivalent R loop creates an n × k logical matrix on every
// iteration.  For 1 M points and k = 15 neighbours that is 15 M allocations
// per iteration repeated n_iter + 1 times, putting heavy pressure on the
// R garbage collector.  The C++ version uses a fixed-size std::vector<int>
// that is reused across all iterations.
//
// Author: Andrew Sánchez Meador, 2026
// License: GPL-3
// ============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "omp.h"   // project-local OpenMP shim (safe when OMP unavailable)
#include <vector>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// C_knn_majority_vote
//
// Parameters
//   nn_idx    – n × k integer matrix of 1-based kNN indices (FNN::knn.index).
//               Row i holds the k nearest neighbour indices for point i.
//   labels_in – length-n logical vector: initial bole label per point.
//   n_iter    – number of majority-vote propagation iterations.
//   ncpu      – OpenMP thread count for the vote-counting inner loop.
//
// Returns a named list:
//   $labels    – length-n logical: converged bole labels after n_iter rounds.
//   $vote_frac – length-n double in [0,1]: fraction of k neighbours that
//                voted TRUE in the final round (BoleProb).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List C_knn_majority_vote(IntegerMatrix nn_idx,
                         LogicalVector labels_in,
                         int n_iter,
                         int ncpu = 1)
{
    const int n = nn_idx.nrow();
    const int k = nn_idx.ncol();
    const int half = k / 2;          // majority threshold: vote_count > half

    // Working copy as plain ints (faster than accessing LogicalVector in hot loop)
    std::vector<int> lbl(n);
    for (int i = 0; i < n; ++i) lbl[i] = labels_in[i] ? 1 : 0;

    std::vector<int> counts(n, 0);   // reused every iteration

    // ------------------------------------------------------------------
    // Propagation iterations — no heap allocation inside the loop
    // ------------------------------------------------------------------
    for (int iter = 0; iter < n_iter; ++iter) {

        // Count bole-labelled neighbours for every point
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(ncpu)
#endif
        for (int i = 0; i < n; ++i) {
            int c = 0;
            for (int j = 0; j < k; ++j) {
                c += lbl[nn_idx(i, j) - 1];  // nn_idx is 1-based (from R/FNN)
            }
            counts[i] = c;
        }

        // Update labels: majority wins
        for (int i = 0; i < n; ++i) {
            lbl[i] = (counts[i] > half) ? 1 : 0;
        }
    }

    // ------------------------------------------------------------------
    // Final vote fractions (BoleProb) — one more count pass, no new labels
    // ------------------------------------------------------------------
    NumericVector vote_frac(n);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(ncpu)
#endif
    for (int i = 0; i < n; ++i) {
        int c = 0;
        for (int j = 0; j < k; ++j) {
            c += lbl[nn_idx(i, j) - 1];
        }
        vote_frac[i] = static_cast<double>(c) / static_cast<double>(k);
    }

    // Convert back to a logical R vector
    LogicalVector labels_out(n);
    for (int i = 0; i < n; ++i) labels_out[i] = (lbl[i] != 0);

    return List::create(Named("labels")    = labels_out,
                        Named("vote_frac") = vote_frac);
}
