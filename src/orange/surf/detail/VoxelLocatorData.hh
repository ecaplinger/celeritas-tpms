
//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/detail/BvhData.hh
//! \todo move to orange/BvhTreeData
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Types.hh"
#include "corecel/cont/EnumArray.hh"
#include "corecel/data/Collection.hh"
#include "geocel/BoundingBox.hh"  // IWYU pragma: keep
#include "orange/surf/detail/MeshData.hh"

namespace celeritas
{

struct SingleVoxel
{
    ItemId<ItemId<Triangle>> start;
    int ntris;
};

//---------------------------------------------------------------------------//
/*!
 * Persistent data used by all BVH trees.
 */
template<Ownership W, MemSpace M>
struct VoxelLocatorData
{
    template<class T>
    using Items = Collection<T, W, M>;

    // Low-level storage
    Items<int> voxels_start;
    Items<int> voxels_ntris;
    Items<int> triangle_map;

    int voxels_x, voxels_y, voxels_z;
    Real3 lower_corner;
    Real3 dimensions;
    Real3 step_dimensions;

    int max_tris;

    //! True if assigned
    explicit CELER_FUNCTION operator bool() const
    {
        return !voxels_start.empty() && !voxels_ntris.empty()
            && !triangle_map.empty();
    }

    //! Assign from another set of data
    template<Ownership W2, MemSpace M2>
    VoxelLocatorData& operator=(VoxelLocatorData<W2, M2> const& other)
    {
        voxels_start = other.voxels_start;
        voxels_ntris = other.voxels_ntris;
        triangle_map = other.triangle_map;

        max_tris = other.max_tris;

        CELER_ENSURE(static_cast<bool>(*this) == static_cast<bool>(other));
        return *this;
    }
};

}
