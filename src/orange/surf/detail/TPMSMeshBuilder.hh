//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/surf/detail/TPMSMeshBuilder.hh
//---------------------------------------------------------------------------//
#pragma once

#include "orange/surf/detail/MarchingCubes.hh"
#include "corecel/data/Collection.hh"
#include "corecel/data/CollectionBuilder.hh"
#include "corecel/Constants.hh"
#include "corecel/Types.hh"
#include "corecel/cont/Array.hh"

namespace celeritas
{

class TPMSMeshBuilder
{
  public:
    explicit TPMSMeshBuilder(real_type(* func)(const Real3&),
            Real3(* grad)(const Real3&),
            int side_steps, real_type offset)
        : iso_mesh_generator_(func, grad,
                side_steps, side_steps, side_steps,
                1.f, 1.f, 1.f)
        , offset_{offset}
    {
    }

    //-----------------------------------------------------------------------//
    /*!
     * Do all the work of generating an isosurface, and offsetting it in both
     * positive and negative directions to 
     *
     * Later, we use the cross product of \f[ (v_2-v_1) \times (v_3-v_1) \f]
     * to determine the normal of the triangle. By swapping \f[ v_2 \f] with
     * \f[ v_3 \f], we make the normals of the positive-offset surface point
     * the opposite direction (back toward the isosurface).
     */
    void operator()(HostVal<MeshData>& offset_data)
    {
        CELER_ENSURE(offset_data.global_verts.empty());
        CELER_ENSURE(offset_data.triangles.empty());

        const size_type ntris = iso_mesh_generator_(zero_data_);

        std::vector<Real3> verts_pos;
        std::vector<Real3> verts_neg;

        std::vector<Triangle> tris_pos;

        const int nverts = zero_data_.global_verts.size();

        verts_neg.reserve(nverts);
        verts_pos.reserve(nverts);

        for (size_type i = 0; i < nverts; i++)
        {
            verts_neg[i] = move_along_normal(
                    zero_data_.global_verts[ItemId<Real3>(i)],
                    zero_data_.global_normals[ItemId<Real3>(i)], -offset_);
            verts_pos[i] = move_along_normal(
                    zero_data_.global_verts[ItemId<Real3>(i)],
                    zero_data_.global_normals[ItemId<Real3>(i)], offset_);
        }

        tris_pos.reserve(ntris);

        // rewind positive verts so they reference positive vertices and point
        // back toward the zero-surface when we take cross-product
        for (size_type i = 0; i < ntris; i++)
        {
            const Triangle& tri_neg = zero_data_.triangles[ItemId<Triangle>(i)];
            tris_pos[i] = Triangle({
                    tri_neg[0] + nverts,
                    tri_neg[2] + nverts,
                    tri_neg[1] + nverts
            });
        }

        CollectionBuilder{&offset_data.global_verts}.insert_back(
                verts_neg.begin(), verts_neg.end());
        CollectionBuilder{&offset_data.global_verts}.insert_back(
                verts_pos.begin(), verts_pos.end());

        offset_data.triangles = zero_data_.triangles;
        CollectionBuilder{&offset_data.triangles}.insert_back(
                tris_pos.begin(), tris_pos.end());
    }
    //-----------------------------------------------------------------------//

  private:
    MarchingCubes<MemSpace::host> iso_mesh_generator_;
    DirectionalMeshData<MemSpace::host> zero_data_;

    /*!
     * Given a point and a normalized direction vector, move the point along
     * the direction vector by a scalar distance and return the resulting
     * point.
     */
    static inline Real3 move_along_normal(const Real3& point,
            const Real3& normal,
            real_type distance)
    {
        return {
            point[0] + normal[0] * distance,
            point[1] + normal[1] * distance,
            point[2] + normal[2] * distance
        };
    }

    /*!
     * Later, we use the cross product of \f[ (v_2-v_1) \times (v_3-v_1) \f]
     * to determine the normal of the triangle. By swapping \f[ v_2 \f] with
     * \f[ v_3 \f], we make the normal point the opposite direction (back
     * toward the isosurface.
     */
    static inline Triangle rewind_triangle(const Triangle& tri)
    {
        return Triangle({ tri.verts[0], tri.verts[2], tri.verts[1] });
    }

    real_type offset_;
};

}
