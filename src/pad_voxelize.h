//  ============================================================================
//  spanner – pad_voxelize.h
//
//  C_voxelize_pad():
//    Memory-efficient, optionally OpenMP-parallel voxelization of a
//    height-normalized lidar point cloud for Beer-Lambert PAD estimation.
//    Reads X/Y/Z by pointer (no intermediate copy), builds compact per-voxel
//    pulse counts via hash maps, and returns P_Intercepted / P_Directed /
//    P_Transmitted vectors that compute_pad_voxels() uses for PAD arithmetic.
//
//  Author: Andrew Sánchez Meador (bi0m3trics), 2026
//  License: GPL-3
//  ============================================================================

#ifndef PAD_VOXELIZE_H
#define PAD_VOXELIZE_H

// Declaration only – implementation + [[Rcpp::export]] in pad_voxelize.cpp

#endif // PAD_VOXELIZE_H
