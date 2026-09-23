//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file corecel/data/CollectionBuilder.hh
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>
#include <algorithm>

#include "corecel/data/Collection.hh"
#include "corecel/math/ArrayUtils.hh"
#include "corecel/math/Constant.hh"
#include "orange/surf/detail/MeshData.hh"
#include "orange/surf/detail/VoxelLocatorData.hh"

namespace celeritas
{

constexpr int default_tracked_triangles_per_voxel = 50;

template <int MaxIntersections,
         int TrackedTrianglesPerVoxel = default_tracked_triangles_per_voxel>
class VoxelAcceleratedMeshIntersector
{
  public:
    using DistancesT = Array<int, MaxIntersections>;
    using NormalsT = Array<Real3, MaxIntersections>;
    using TriangleVertsT = Array<Real3, 3>;
    using TrianglesT = Array<TriangleVertsT, MaxIntersections>;
    using IndexTrackerT = Array<int, TrackedTrianglesPerVoxel>;

    explicit CELER_FUNCTION VoxelAcceleratedMeshIntersector(
            DeviceCRef<MeshData> mesh,
            DeviceCRef<VoxelLocatorData> locator)
        : mesh_{mesh}
        , locator_{locator}
    {
    }

    //! Find triangles that intersect with ray
    CELER_FUNCTION int operator()(
            Real3 start, const Real3& direction,
            DistancesT& distances, NormalsT& normals)
    {
        IndexTrackerT prev_tested_tris;

        // direction should be a unit vector
        CELER_ENSURE(std::fabs(direction[0]*direction[0]
                    + direction[1]*direction[1]
                    + direction[2]*direction[2] - 1) < eps);

        if (!in_bounds(start)) {
            const real_type dist = distance_to_cell_edge(start, direction);
            if (!std::isfinite(dist))
                return 0;

            start[0] += dist * direction[0];
            start[1] += dist * direction[1];
            start[2] += dist * direction[2];
        }

        // get corners of first voxel
        const Array<int, 3> first_voxel_idx = get_voxel_by_point(start);

        Real3 d_edge, t_max;
        Array<int, 3> idx_step;
        for (int dim = 0; dim < 3; dim++)
        {
            real_type first_lo = first_voxel_idx[dim]
                * locator_.step_dimensions[dim];
            real_type first_hi = (first_voxel_idx[dim]+1)
                * locator_.step_dimensions[dim];

            if (std::fabs(direction[dim]) < eps)
                t_max[dim] = inf;

            d_edge[dim] = std::fabs(
                    locator_.step_dimensions[dim] / direction[dim]);
            t_max[dim] = direction[dim] > 0 ?
                (first_hi - start[dim]) / direction[dim] :
                (first_lo - start[dim]) / direction[dim];
            idx_step[dim] = direction[dim] > 0 ? 1 : 0;
        }

        Array<int, 3> voxel_coord = first_voxel_idx;
        int num_intersections = 0;
        while (voxel_in_bounds(voxel_coord))
        {
            IndexTrackerT tested_tris;
            num_intersections += intersect_within_voxel(
                    voxel_coord, start, direction,
                    normals, distances, num_intersections,
                    tested_tris, prev_tested_tris);

            prev_tested_tris = std::move(tested_tris);

            if (num_intersections == MaxIntersections)
                break;

            auto idx = std::min_element(d_edge.begin(), d_edge.end())
                - d_edge.begin();
            voxel_coord[idx] += idx_step[idx];
            t_max[idx] += d_edge[idx];
        }

        return num_intersections;
    }

  private:

    DeviceCRef<MeshData> mesh_;
    DeviceCRef<VoxelLocatorData> locator_;

    const real_type inf = std::numeric_limits<real_type>::infinity();
    const real_type eps = std::numeric_limits<real_type>::epsilon();

    /*!
     * For a particular voxel, loop through all child triangles and see which
     * ones our ray intersects with. Avoid repeating work by keeping a list of
     * already tested triangles and testing against the list of triangles
     * tested within the previous voxel. Since all tested voxels are adjacent,
     * this should never test the same triangle twice.
     */
    CELER_FORCEINLINE_FUNCTION int intersect_within_voxel(
            const Array<int, 3>& voxel,
            const Real3& start,
            const Real3& direction,
            NormalsT& norm_output,
            DistancesT& dist_output,
            int found_already,
            IndexTrackerT& tested_tris,
            IndexTrackerT& prev_tested_tris)
    {
        tested_tris.fill(-1);
        int found = 0;

        int flat_idx = voxel_to_index(voxel);
        const int ntris = locator_.voxels_ntris[ItemId<int>(flat_idx)];
        const int start_idx = locator_.voxels_start[ItemId<int>(flat_idx)];

        for (int i = 0; i < ntris; i++)
        {
            if (found_already == MaxIntersections)
                return found;

            const int map_idx = start_idx + i;
            const int tri_idx = locator_.triangle_map[ItemId<int>(map_idx)];

            if (i < TrackedTrianglesPerVoxel)
                tested_tris[i] = tri_idx;

            for (int j = 0; j < TrackedTrianglesPerVoxel; j++)
                if (prev_tested_tris[j] == tri_idx)
                    continue;

            const Triangle& tri = mesh_.triangles[ItemId<Triangle>(tri_idx)];
            Array<Real3, 3> tri_verts = {
                mesh_.global_verts[ItemId<Real3>(tri.verts[0])],
                mesh_.global_verts[ItemId<Real3>(tri.verts[1])],
                mesh_.global_verts[ItemId<Real3>(tri.verts[2])]
            };

            Real3 norm;
            real_type dist = run_moller_trumbore(tri_verts, start, direction, &norm);

            if (std::isfinite(dist))
            {
                norm_output[found_already] = norm;
                dist_output[found_already++] = dist;
                found++;
            }
        }

        return found;
    }

    /*!
     * Run the Moller-Trumbore triangle intersection algorithm on a triangle.
     * Probably could be moved into the Algorithms file or something...
     * Adapted from https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
     */
    CELER_FORCEINLINE_FUNCTION static real_type run_moller_trumbore(
            const Array<Real3, 3>& vertices,
            const Real3& start, const Real3& direction,
            Real3* norm)
    {
        Real3 edge1, edge2, s;

        constexpr real_type inf = std::numeric_limits<real_type>::infinity();
        constexpr real_type eps = std::numeric_limits<real_type>::epsilon();

        for (int dim = 0; dim < 3; dim++)
        {
            edge1[dim] = vertices[1][dim] - vertices[0][dim];
            edge2[dim] = vertices[2][dim] - vertices[0][dim];
            s[dim] = start[dim] - vertices[0][dim];
        }

        Real3 ray_cross_e2 = cross_product(direction, edge2);
        real_type det = dot_product(edge1, ray_cross_e2);

        if (std::fabs(det) < eps) return inf;

        real_type inv_det = 1.0 / det;
        real_type u = inv_det * dot_product(s, ray_cross_e2);

        if (u < -eps || u - 1 > eps) return inf;

        Real3 s_cross_e1 = cross_product(s, edge1);
        real_type v = inv_det * dot_product(direction, s_cross_e1);

        if (v < -eps || u + v - 1 > eps) return inf;

        real_type t = inv_det * dot_product(edge2, s_cross_e1);

        if (t > eps) {
            Real3 normal = make_unit_vector(cross_product(edge1, edge2));
            *norm = normal;
            return t;
        }
        return inf;
    }

    /*!
     * Given a ray (start and direction), get the minimum distance before it
     * hits a boundary of the voxel grid.
     */
    CELER_FORCEINLINE_FUNCTION real_type distance_to_cell_edge(Real3 start,
            Real3 direction) const
    {
        //const real_type inf = std::numeric_limits<real_type>::infinity();
        real_type distance = inf;

        for (int dim = 0; dim < 3; dim++)
        {
            if (direction[dim] != 0)
            {
                const real_type dmin = locator_.lower_corner[dim];
                const real_type dmax = dmin + locator_.dimensions[dim];
                const real_type d_dmin = (dmin - start[dim]) / direction[dim];
                const real_type d_dmax = (dmax - start[dim]) / direction[dim];

                if (d_dmin >= 0)
                    distance = std::min(distance, d_dmin);
                if (d_dmax >= 0)
                    distance = std::min(distance, d_dmax);
            }
        }

        return distance;
    }

    /*!
     * Determine whether a point is within the bounds of the voxel grid.
     */
    CELER_FORCEINLINE_FUNCTION bool in_bounds(Real3 point)
    {
        const Real3& lower_corner = locator_.lower_corner;
        const Real3& dimensions = locator_.dimensions;

        return (point[0] >= lower_corner[0] &&
                point[0] <= lower_corner[0] + dimensions[0] &&
                point[1] >= lower_corner[1] &&
                point[1] <= lower_corner[1] + dimensions[0] &&
                point[2] >= lower_corner[2] &&
                point[2] <= lower_corner[2] + dimensions[0]);
    }

    /*!
     * Given a point within our voxel universe, get the position of the voxel
     * that contains that point.
     */
    CELER_FORCEINLINE_FUNCTION Array<int, 3> get_voxel_by_point(Real3 point)
    {
        return {
            (int) std::floor(point[0] / locator_.step_dimensions[0]),
            (int) std::floor(point[1] / locator_.step_dimensions[1]),
            (int) std::floor(point[2] / locator_.step_dimensions[2])
        };
    }

    /*!
     * Given 3-dimensional voxel coordinates, get an index into the flat array
     * of data within this class's VoxelLocatorData.
     */
    CELER_FORCEINLINE_FUNCTION int voxel_to_index(Array<int, 3> voxel)
    {
        return voxel[0] + voxel[1]*locator_.voxels_x
            + voxel[2]*locator_.voxels_x*locator_.voxels_y;
    }

    /*!
     * Given 3-dimensional voxel coordinates, tell whether they are within the
     * voxel grid.
     */
    CELER_FORCEINLINE_FUNCTION bool voxel_in_bounds(Array<int, 3> voxel)
    {
        return (voxel[0] >= 0 && voxel[0] < locator_.voxels_x
                && voxel[1] >= 0 && voxel[1] < locator_.voxels_y
                && voxel[2] >= 0 && voxel[2] < locator_.voxels_z);
    }
};

}
