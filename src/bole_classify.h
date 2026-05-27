//  ============================================================================
//  spanner – bole_classify.h
//
//  Batch RANSAC circle fitting across height slices for individual trees.
//  Used by the bole-segmentation pipeline:
//    classify_stem_points() → segment_bole() → compute_bole_volume()
//    → assess_tree_quality() → refit_trees()
//
//  Author: Andrew Sánchez Meador (bi0m3trics), 2026
//  License: GPL-3
//  ============================================================================

#ifndef BOLE_CLASSIFY_HPP
#define BOLE_CLASSIFY_HPP

#include "algorithms.h"   // ransacCircle(), rmatrix2cpp(), etc.

// Declared here; implementation + Rcpp::export tag in bole_classify.cpp
// C_ransac_bole_slices  – see bole_classify.cpp

#endif // BOLE_CLASSIFY_HPP
