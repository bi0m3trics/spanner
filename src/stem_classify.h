//  ============================================================================
//  spanner – stem_classify.h
//
//  Batch RANSAC circle fitting across height slices for individual trees.
//  Used by the stem-segmentation pipeline:
//    classify_stem_points() → segment_stem() → compute_stem_volume()
//    → assess_tree_quality() → refit_trees()
//
//  Author: Andrew Sánchez Meador (bi0m3trics), 2026
//  License: GPL-3
//  ============================================================================

#ifndef STEM_CLASSIFY_HPP
#define STEM_CLASSIFY_HPP

#include "algorithms.h"   // ransacCircle(), rmatrix2cpp(), etc.

// Declared here; implementation + Rcpp::export tag in stem_classify.cpp
// C_ransac_stem_slices  – see stem_classify.cpp

#endif // STEM_CLASSIFY_HPP
