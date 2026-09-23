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
#include "orange/surf/detail/MeshData.hh"
#include "orange/surf/detail/VoxelLocatorData.hh"

namespace celeritas
{

struct Voxel
{
    int start;
    int ntris;

    explicit Voxel(int start, int ntris)
        : start{start}, ntris{ntris}
    {
    }
};

class VoxelLocatorBuilder
{
  public:
    using Storage = VoxelLocatorData<Ownership::value, MemSpace::host>;

    // TODO: using a raw pointer is too dirty for modern C++, but not sure what
    // celeritas pointer conventions are
    explicit VoxelLocatorBuilder(
            Storage* storage,
            int voxels_x, int voxels_y, int voxels_z,
            real_type length_x, real_type length_y, real_type length_z,
            Real3 lower_point, Real3 upper_point)
        : triangle_map_{&storage->triangle_map}
        , voxels_start_{&storage->voxels_start}
        , voxels_ntris_{&storage->voxels_ntris}
        , max_tris_{&storage->max_tris}
        , voxels_x_{storage->voxels_x}
        , voxels_y_{storage->voxels_y}
        , voxels_z_{storage->voxels_z}
        , lower_corner_{storage->lower_corner}
        , dimensions_{storage->dimensions}
        , step_dimensions_{storage->dimensions}
    {
    }

    // NOTE: the way this works currently, it assumes that no triangle will
    // span more than 8 voxels; for a more generalizable VoxelLocator, this
    // algorithm must be updated
    void assign_triangles(HostVal<MeshData>& data)
    {
        std::vector<std::vector<int>> local_map;
        int num_voxels = voxels_x_ * voxels_y_ * voxels_z_;
        local_map.reserve(num_voxels);

        for (int i = 0; i < data.triangles.size(); i++)
        {
            std::vector<int> voxels_added;
            Array<Real3, 3> tri_verts = data.triangle_at(i);
            Real3 bbox_lo = {
                std::min({ tri_verts[0][0], tri_verts[1][0], tri_verts[2][0] }),
                std::min({ tri_verts[0][1], tri_verts[1][1], tri_verts[2][1] }),
                std::min({ tri_verts[0][2], tri_verts[1][2], tri_verts[2][2] })
            };
            Real3 bbox_hi = {
                std::max({ tri_verts[0][0], tri_verts[1][0], tri_verts[2][0] }),
                std::max({ tri_verts[0][1], tri_verts[1][1], tri_verts[2][1] }),
                std::max({ tri_verts[0][2], tri_verts[1][2], tri_verts[2][2] })
            };

            for (int p = 0; p < 8; p++)
            {
                Real3 bbox_p = {
                    p&1 ? bbox_hi[0] : bbox_lo[0],
                    (p>>1)&1 ? bbox_hi[1] : bbox_lo[1],
                    (p>>2)&1 ? bbox_hi[2] : bbox_lo[2]
                };
                int voxel_p = voxel_to_index(get_point_voxel(bbox_p));

                if (std::find(voxels_added.begin(), voxels_added.end(), voxel_p)
                        == voxels_added.end())
                {
                    voxels_added.push_back(voxel_p);
                    local_map[voxel_p].push_back(i);
                }
            }
        }

        int max_voxel_tris = 0;
        std::vector<int> global_map;
        std::vector<int> global_start(num_voxels);
        std::vector<int> global_ntris(num_voxels);
        int start = 0;
        for (int i = 0; i < num_voxels; i++)
        {
            const int voxel_tris = local_map[i].size();
            for (int j = 0; j < voxel_tris; j++)
                global_map.push_back(local_map[i][j]);
            global_start[i] = start;
            global_ntris[i] = voxel_tris;
            start += voxel_tris;

            if (voxel_tris > max_voxel_tris)
                max_voxel_tris = voxel_tris;
        }

        triangle_map_.insert_back(global_map.begin(), global_map.end());
        voxels_start_.insert_back(global_start.begin(), global_start.end());
        voxels_ntris_.insert_back(global_ntris.begin(), global_ntris.end());

        *max_tris_ = max_voxel_tris;
    }

    CELER_FUNCTION int voxels_along_ray(Real3 start, Real3 direction)
    {

    }

    CELER_FUNCTION bool in_bounds(Real3 point)
    {
        return (point[0] >= lower_corner_[0] &&
                point[0] <= lower_corner_[0] + step_dimensions_[0] &&
                point[1] >= lower_corner_[1] &&
                point[1] <= lower_corner_[1] + step_dimensions_[1] &&
                point[2] >= lower_corner_[2] &&
                point[2] <= lower_corner_[2] + step_dimensions_[2]);
    }

private:
    MeshData<MemSpace::host>* data_;

    const int voxels_x_, voxels_y_, voxels_z_;
    const Real3 lower_corner_;
    const Real3 dimensions_;
    const Real3 step_dimensions_;

    // Accessors to VoxelLocatorData members
    CollectionBuilder<int> triangle_map_;
    CollectionBuilder<int> voxels_start_;
    CollectionBuilder<int> voxels_ntris_;
    int* max_tris_;

    inline int voxel_to_index(Array<int, 3> voxel)
    {
        return voxel[0] + voxel[1]*voxels_x_ + voxel[2]*voxels_x_*voxels_y_;
    }

    inline Array<int, 3> get_point_voxel(Real3 point)
    {
        return {
            (int) std::floor(point[0] / step_dimensions_[0]),
            (int) std::floor(point[1] / step_dimensions_[1]),
            (int) std::floor(point[2] / step_dimensions_[2])
        };
    }
};

}
