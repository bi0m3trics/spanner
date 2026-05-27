//  ============================================================================
//  spanner – pad_voxelize.cpp
//
//  C_voxelize_pad():
//    Voxelize a height-normalized lidar point cloud and compute per-voxel
//    pulse-tracking counts (P_Intercepted, P_Directed, P_Transmitted) for the
//    Beer-Lambert plant area density (PAD) pipeline in compute_pad_voxels().
//
//  Algorithm
//  ---------
//  1. Parallel voxelization  (per-thread hash maps, then serial merge)
//     Each point independently maps to (ix, iy, iz) = floor(coord/vox_size).
//     Per-thread std::unordered_maps accumulate P_Intercepted (3-D) and
//     column totals (2-D).  Thread maps are merged serially.  Memory peak:
//     O(unique_voxels) rather than O(n_points) as in the pure-R version.
//
//  2. Flatten hash map -> compact VoxEntry vector
//     Only occupied voxels are stored; size << n_points for typical TLS/MLS.
//
//  3. Sort by (ix, iy, iz DESC) — top-down column order for cumsum.
//     Only the compact voxel vector is sorted, not the point table.
//
//  4. Parallel per-column cumulative sum  (OpenMP over independent columns)
//     P_Directed[i]    = col_total - cumint_up_to_and_including_i + P_Int[i]
//     P_Transmitted[i] = col_total - cumint_up_to_and_including_i
//
//  5. Build output List (vx, vy, vz, P_Intercepted, P_Directed, P_Transmitted)
//     Voxel-centre coordinates: ix * vox_size + vox_size/2, etc.
//
//  Author: Andrew Sánchez Meador (bi0m3trics), 2026
//  License: GPL-3
//  ============================================================================

#include <Rcpp.h>
#include "omp.h"

#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// 3-D voxel key: three int32_t indices packed into a struct.
// The hash uses boost-style hash_combine with a Wang-scrambled int32 mixer so
// that adjacent voxels produce well-distributed bucket assignments.
// ---------------------------------------------------------------------------
struct VoxKey3 {
    int32_t ix, iy, iz;
    bool operator==(const VoxKey3& o) const noexcept {
        return ix == o.ix && iy == o.iy && iz == o.iz;
    }
};

static inline uint32_t scramble32(int32_t v) noexcept {
    // Wang hash on reinterpreted bit pattern
    uint32_t u = static_cast<uint32_t>(v);
    u ^= u >> 16;
    u *= 0x45d9f3bu;
    u ^= u >> 16;
    return u;
}

struct VoxKey3Hash {
    size_t operator()(const VoxKey3& k) const noexcept {
        // boost hash_combine pattern
        size_t h = static_cast<size_t>(scramble32(k.ix));
        h ^= static_cast<size_t>(scramble32(k.iy)) + 0x9e3779b9u + (h << 6) + (h >> 2);
        h ^= static_cast<size_t>(scramble32(k.iz)) + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
    }
};

// 2-D column key: pack uint32 bit-patterns of (ix, iy) into uint64_t.
// Casting via uint32_t first makes the shift well-defined for negative indices.
static inline uint64_t col_key(int32_t ix, int32_t iy) noexcept {
    return (static_cast<uint64_t>(static_cast<uint32_t>(ix)) << 32) |
            static_cast<uint64_t>(static_cast<uint32_t>(iy));
}

using Map3D = std::unordered_map<VoxKey3,  int32_t, VoxKey3Hash>;
using Map2D = std::unordered_map<uint64_t, int32_t>;

// ---------------------------------------------------------------------------
// C_voxelize_pad
//
// Parameters
//   x, y, z    – first-return point coordinates (height-normalized, filtered
//                to min_height already on the R side)
//   vox_size   – cubic voxel edge length (m)
//   ncpu       – number of OpenMP threads (1 = serial)
//
// Returns a named List with six equally-sized vectors:
//   vx, vy, vz          – voxel-centre coordinates (double)
//   P_Intercepted        – first returns whose Z falls in this voxel (int)
//   P_Directed           – cumulative first returns at height >= voxel bottom (int)
//   P_Transmitted        – P_Directed - P_Intercepted (int)
//
// Rows are ordered by (vx ASC, vy ASC, vz DESC) — top-down within column.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List C_voxelize_pad(NumericVector x,
                    NumericVector y,
                    NumericVector z,
                    double vox_size,
                    int    ncpu = 1)
{
    if (x.size() != y.size() || x.size() != z.size())
        stop("x, y, z must have the same length.");
    if (vox_size <= 0.0)
        stop("vox_size must be positive.");

    const int n = static_cast<int>(x.size());
    if (n == 0) {
        return List::create(Named("vx")            = NumericVector(0),
                            Named("vy")            = NumericVector(0),
                            Named("vz")            = NumericVector(0),
                            Named("P_Intercepted") = IntegerVector(0),
                            Named("P_Directed")    = IntegerVector(0),
                            Named("P_Transmitted") = IntegerVector(0));
    }

    const double* px = x.begin();
    const double* py = y.begin();
    const double* pz = z.begin();

    // -----------------------------------------------------------------------
    // Resolve thread count
    // -----------------------------------------------------------------------
    int nthreads = 1;
    int original_threads = 1;
#ifdef _OPENMP
    original_threads = omp_get_max_threads();
    nthreads = (ncpu > 0) ? std::min(ncpu, original_threads) : original_threads;
    omp_set_max_active_levels(1);
#endif

    // -----------------------------------------------------------------------
    // Step 1: Voxelization – build 3-D count map and 2-D column-total map.
    //
    // With OpenMP: each thread accumulates into private maps; serial merge.
    // Serial fallback: single pass directly into the shared maps.
    // -----------------------------------------------------------------------
    Map3D count3d;
    Map2D col_total;

#ifdef _OPENMP
    {
        std::vector<Map3D> t3(nthreads);
        std::vector<Map2D> t2(nthreads);

        // Rough pre-allocation: assume ~4 points per 3-D voxel on average
        {
            size_t hint3 = static_cast<size_t>(n) / (4 * nthreads);
            size_t hint2 = static_cast<size_t>(n) / (8 * nthreads);
            if (hint3 < 16) hint3 = 16;
            if (hint2 < 16) hint2 = 16;
            for (int t = 0; t < nthreads; ++t) {
                t3[t].reserve(hint3);
                t2[t].reserve(hint2);
            }
        }

        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (int i = 0; i < n; ++i) {
            const int tid = omp_get_thread_num();
            const int32_t ix = static_cast<int32_t>(std::floor(px[i] / vox_size));
            const int32_t iy = static_cast<int32_t>(std::floor(py[i] / vox_size));
            const int32_t iz = static_cast<int32_t>(std::floor(pz[i] / vox_size));
            t3[tid][{ix, iy, iz}]++;
            t2[tid][col_key(ix, iy)]++;
        }

        // Serial merge: move thread-0 maps, fold in remaining threads
        count3d   = std::move(t3[0]);
        col_total = std::move(t2[0]);
        for (int t = 1; t < nthreads; ++t) {
            for (auto& kv : t3[t]) count3d[kv.first]   += kv.second;
            for (auto& kv : t2[t]) col_total[kv.first] += kv.second;
        }
    }
#else
    {
        const size_t hint3 = static_cast<size_t>(n) / 4;
        const size_t hint2 = static_cast<size_t>(n) / 8;
        count3d.reserve(hint3 > 16 ? hint3 : 16);
        col_total.reserve(hint2 > 16 ? hint2 : 16);

        for (int i = 0; i < n; ++i) {
            const int32_t ix = static_cast<int32_t>(std::floor(px[i] / vox_size));
            const int32_t iy = static_cast<int32_t>(std::floor(py[i] / vox_size));
            const int32_t iz = static_cast<int32_t>(std::floor(pz[i] / vox_size));
            count3d[{ix, iy, iz}]++;
            col_total[col_key(ix, iy)]++;
        }
    }
#endif

    // -----------------------------------------------------------------------
    // Step 2: Flatten 3-D map into a compact vector of VoxEntry structs.
    //         Store col_total alongside each entry so the cumsum step needs
    //         no further map lookups.
    // -----------------------------------------------------------------------
    struct VoxEntry {
        int32_t ix, iy, iz;
        int32_t P_Intercepted;
        int32_t col_total_val;
    };

    const size_t nv = count3d.size();
    std::vector<VoxEntry> voxels;
    voxels.reserve(nv);

    for (auto& kv : count3d) {
        const int32_t ct = col_total.at(col_key(kv.first.ix, kv.first.iy));
        voxels.push_back({kv.first.ix, kv.first.iy, kv.first.iz,
                          kv.second, ct});
    }

    // Free hash maps immediately to recover peak memory before output allocation
    { Map3D tmp; count3d.swap(tmp); }
    { Map2D tmp; col_total.swap(tmp); }

    // -----------------------------------------------------------------------
    // Step 3: Sort voxels by (ix ASC, iy ASC, iz DESC) — top-down within
    //         each column.  Only the compact voxel vector is sorted, not the
    //         original point table.
    // -----------------------------------------------------------------------
    std::sort(voxels.begin(), voxels.end(),
              [](const VoxEntry& a, const VoxEntry& b) {
                  if (a.ix != b.ix) return a.ix < b.ix;
                  if (a.iy != b.iy) return a.iy < b.iy;
                  return a.iz > b.iz;  // descending Z → top to bottom
              });

    // -----------------------------------------------------------------------
    // Step 4: Per-column cumulative sum to derive P_Directed / P_Transmitted.
    //
    //   Sorted top-to-bottom within each (ix, iy) column:
    //     cumint  — running sum of P_Intercepted down to and including row i
    //
    //   P_Directed[i]    = col_total - cumint + P_Int[i]
    //                    (= pulses directed at this voxel from above)
    //   P_Transmitted[i] = col_total - cumint
    //                    (= pulses that pass through without interception)
    //
    //   Verification (simple 3-voxel column: P_Int = 5, 3, 2; total = 10):
    //     vz=top:  cumint=5, PDir=10-5+5=10✓  PTrans=10-5=5✓
    //     vz=mid:  cumint=8, PDir=10-8+3=5✓   PTrans=10-8=2✓
    //     vz=bot:  cumint=10,PDir=10-10+2=2✓  PTrans=10-10=0✓
    //
    //   Columns are independent → parallelise over block indices with OpenMP.
    // -----------------------------------------------------------------------

    // Find contiguous-block boundaries (each block is one vertical column)
    const int inv = static_cast<int>(nv);
    std::vector<int> bstart;
    bstart.reserve(nv / 4 + 1);
    bstart.push_back(0);
    for (int i = 1; i < inv; ++i) {
        if (voxels[i].ix != voxels[i - 1].ix ||
            voxels[i].iy != voxels[i - 1].iy)
            bstart.push_back(i);
    }
    bstart.push_back(inv);

    const int n_blocks = static_cast<int>(bstart.size()) - 1;

    std::vector<int32_t> P_Directed(nv);
    std::vector<int32_t> P_Transmitted(nv);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 32) num_threads(nthreads)
#endif
    for (int b = 0; b < n_blocks; ++b) {
        const int start = bstart[b];
        const int end   = bstart[b + 1];
        int32_t cumint  = 0;
        for (int i = start; i < end; ++i) {
            cumint            += voxels[i].P_Intercepted;
            P_Directed[i]    = voxels[i].col_total_val - cumint + voxels[i].P_Intercepted;
            P_Transmitted[i] = voxels[i].col_total_val - cumint;
        }
    }

    // Restore original thread count
#ifdef _OPENMP
    omp_set_num_threads(original_threads);
#endif

    // -----------------------------------------------------------------------
    // Step 5: Build output vectors (voxel-centre coordinates + pulse counts).
    // -----------------------------------------------------------------------
    const double hs = vox_size / 2.0;

    NumericVector out_vx(nv), out_vy(nv), out_vz(nv);
    IntegerVector out_pi(nv), out_pd(nv), out_pt(nv);

    for (int i = 0; i < inv; ++i) {
        out_vx[i] = voxels[i].ix * vox_size + hs;
        out_vy[i] = voxels[i].iy * vox_size + hs;
        out_vz[i] = voxels[i].iz * vox_size + hs;
        out_pi[i] = voxels[i].P_Intercepted;
        out_pd[i] = P_Directed[i];
        out_pt[i] = P_Transmitted[i];
    }

    return List::create(Named("vx")            = out_vx,
                        Named("vy")            = out_vy,
                        Named("vz")            = out_vz,
                        Named("P_Intercepted") = out_pi,
                        Named("P_Directed")    = out_pd,
                        Named("P_Transmitted") = out_pt);
}

// ============================================================================
// C_voxelize_pad_rays
//
// 3-D ray-tracing voxelizer for TLS or trajectory-ALS data.
//
// For each first-return point (x[i], y[i], z[i]), a ray is cast from the
// corresponding scanner position (sx[i], sy[i], sz[i]) and walked through
// the voxel grid using the Amanatides & Woo (1987) DDA algorithm.
//
// Per-voxel counters:
//   P_Directed      – incremented for every voxel the ray passes through or
//                     terminates in.
//   P_Intercepted   – incremented only for the terminal voxel (return point).
//   P_Transmitted   – incremented for every non-terminal voxel on the ray
//                     (P_Directed - P_Intercepted).
//
// Note: For present voxels, P_Directed = P_Intercepted + P_Transmitted by
// construction, so the occluded fraction (1 - (I+T)/D) is always 0.
// Truly occluded voxels receive zero directed pulses and are absent from the
// output table — the same behaviour as the column model.
//
// Parameters
//   x, y, z        – first-return point coordinates (n elements)
//   sx, sy, sz     – per-point scanner positions (n elements).
//                    For TLS: repeat a single fixed position for all points.
//                    For ALS: supply per-pulse positions interpolated from the
//                    flight trajectory (R-side) and Z-clipped to the canopy
//                    top to avoid traversing thousands of empty-air voxels.
//   vox_size       – cubic voxel edge length (m)
//   ncpu           – number of OpenMP threads
//
// Returns the same named List as C_voxelize_pad.
// ============================================================================

// ---------------------------------------------------------------------------
// Amanatides & Woo (1987) 3-D DDA helper.
// Walks from scanner (sx,sy,sz) to return point (ex,ey,ez) through a voxel
// grid of cell size vs, invoking fn(ix, iy, iz, is_terminal) for each voxel.
// Both endpoints are visited.  The guard of 100 000 steps prevents infinite
// loops on degenerate inputs.
// ---------------------------------------------------------------------------
template<typename Fn>
static inline void dda_walk(double sx, double sy, double sz,
                              double ex, double ey, double ez,
                              double vs, Fn fn)
{
    int32_t ix = static_cast<int32_t>(std::floor(sx / vs));
    int32_t iy = static_cast<int32_t>(std::floor(sy / vs));
    int32_t iz = static_cast<int32_t>(std::floor(sz / vs));
    const int32_t tx = static_cast<int32_t>(std::floor(ex / vs));
    const int32_t ty = static_cast<int32_t>(std::floor(ey / vs));
    const int32_t tz = static_cast<int32_t>(std::floor(ez / vs));

    if (ix == tx && iy == ty && iz == tz) { fn(ix, iy, iz, true);  return; }

    const double dx = ex - sx, dy = ey - sy, dz = ez - sz;
    const double len = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (len < 1e-12) { fn(ix, iy, iz, true); return; }

    const double ndx = dx / len, ndy = dy / len, ndz = dz / len;
    const int stepx = (ndx >= 0) ? 1 : -1;
    const int stepy = (ndy >= 0) ? 1 : -1;
    const int stepz = (ndz >= 0) ? 1 : -1;

    // Distance (world units) to cross one full voxel along each axis.
    const double tDX = (std::fabs(ndx) > 1e-15) ? vs / std::fabs(ndx) : 1e30;
    const double tDY = (std::fabs(ndy) > 1e-15) ? vs / std::fabs(ndy) : 1e30;
    const double tDZ = (std::fabs(ndz) > 1e-15) ? vs / std::fabs(ndz) : 1e30;

    // Distance to first boundary crossing in each axis.
    auto first_t = [&](int32_t idx, int step, double coord, double nd) -> double {
        if (std::fabs(nd) <= 1e-15) return 1e30;
        const double next = (step > 0) ? static_cast<double>(idx + 1) * vs
                                       : static_cast<double>(idx)     * vs;
        return (next - coord) / nd;   // sign matches by construction
    };
    double tMaxX = first_t(ix, stepx, sx, ndx);
    double tMaxY = first_t(iy, stepy, sy, ndy);
    double tMaxZ = first_t(iz, stepz, sz, ndz);

    for (int g = 0; g < 100000; ++g) {
        const bool terminal = (ix == tx && iy == ty && iz == tz);
        fn(ix, iy, iz, terminal);
        if (terminal) break;
        if (tMaxX <= tMaxY && tMaxX <= tMaxZ) { ix += stepx; tMaxX += tDX; }
        else if (tMaxY  <= tMaxZ)              { iy += stepy; tMaxY += tDY; }
        else                                   { iz += stepz; tMaxZ += tDZ; }
    }
}

// ---------------------------------------------------------------------------
// Per-voxel accumulator for the ray-tracing path.
// ---------------------------------------------------------------------------
struct RayVoxData {
    int32_t directed    = 0;
    int32_t intercepted = 0;
    int32_t transmitted = 0;
};
using RayMap3D = std::unordered_map<VoxKey3, RayVoxData, VoxKey3Hash>;

// [[Rcpp::export]]
List C_voxelize_pad_rays(NumericVector x,  NumericVector y,  NumericVector z,
                          NumericVector sx, NumericVector sy, NumericVector sz,
                          double vox_size,
                          int    ncpu = 1)
{
    if (x.size() != y.size() || x.size() != z.size() ||
        x.size() != sx.size() || x.size() != sy.size() || x.size() != sz.size())
        stop("Point and scanner-position vectors must all have the same length.");
    if (vox_size <= 0.0)
        stop("vox_size must be positive.");

    const int n = static_cast<int>(x.size());
    if (n == 0) {
        return List::create(Named("vx")            = NumericVector(0),
                            Named("vy")            = NumericVector(0),
                            Named("vz")            = NumericVector(0),
                            Named("P_Intercepted") = IntegerVector(0),
                            Named("P_Directed")    = IntegerVector(0),
                            Named("P_Transmitted") = IntegerVector(0));
    }

    int nthreads = 1;
    int original_threads = 1;
#ifdef _OPENMP
    original_threads = omp_get_max_threads();
    nthreads = (ncpu > 0) ? std::min(ncpu, original_threads) : original_threads;
    omp_set_max_active_levels(1);
#endif

    const double* px  = x.begin();
    const double* py  = y.begin();
    const double* pz  = z.begin();
    const double* psx = sx.begin();
    const double* psy = sy.begin();
    const double* psz = sz.begin();

    // Per-thread accumulation maps
    std::vector<RayMap3D> tmaps(nthreads);
    {
        const size_t hint = static_cast<size_t>(n) / std::max(1, nthreads) / 2;
        for (auto& m : tmaps) m.reserve(hint > 16 ? hint : 16);
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64) num_threads(nthreads)
#endif
    for (int i = 0; i < n; ++i) {
#ifdef _OPENMP
        RayMap3D& m = tmaps[omp_get_thread_num()];
#else
        RayMap3D& m = tmaps[0];
#endif
        dda_walk(psx[i], psy[i], psz[i],
                 px[i],  py[i],  pz[i],
                 vox_size,
                 [&](int32_t ix, int32_t iy, int32_t iz, bool terminal) {
                     RayVoxData& v = m[{ix, iy, iz}];
                     v.directed++;
                     if (terminal) v.intercepted++;
                     else          v.transmitted++;
                 });
    }

    // Serial merge: fold threads 1..N-1 into thread 0's map
    for (int t = 1; t < nthreads; ++t) {
        for (auto& kv : tmaps[t]) {
            RayVoxData& dst = tmaps[0][kv.first];
            dst.directed    += kv.second.directed;
            dst.intercepted += kv.second.intercepted;
            dst.transmitted += kv.second.transmitted;
        }
        { RayMap3D tmp; tmaps[t].swap(tmp); }   // free thread map memory
    }

#ifdef _OPENMP
    omp_set_num_threads(original_threads);
#endif

    const RayMap3D& merged = tmaps[0];
    const size_t nv = merged.size();
    const double hs = vox_size / 2.0;

    NumericVector out_vx(nv), out_vy(nv), out_vz(nv);
    IntegerVector out_pi(nv), out_pd(nv), out_pt(nv);

    size_t idx = 0;
    for (const auto& kv : merged) {
        out_vx[idx] = kv.first.ix * vox_size + hs;
        out_vy[idx] = kv.first.iy * vox_size + hs;
        out_vz[idx] = kv.first.iz * vox_size + hs;
        out_pd[idx] = kv.second.directed;
        out_pi[idx] = kv.second.intercepted;
        out_pt[idx] = kv.second.transmitted;
        ++idx;
    }

    return List::create(Named("vx")            = out_vx,
                        Named("vy")            = out_vy,
                        Named("vz")            = out_vz,
                        Named("P_Intercepted") = out_pi,
                        Named("P_Directed")    = out_pd,
                        Named("P_Transmitted") = out_pt);
}
