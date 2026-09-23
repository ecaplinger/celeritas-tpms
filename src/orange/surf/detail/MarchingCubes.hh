//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/surf/detail/TPMSDefinition.hh
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>

#include "corecel/Constants.hh"
#include "corecel/Types.hh"
#include "corecel/cont/Array.hh"
#include "corecel/math/NewtonSolver.hh"
#include "corecel/math/Interval.hh"
#include "corecel/data/CollectionBuilder.hh"

#include "orange/surf/detail/MeshData.hh"
#include "orange/surf/detail/MarchingCubesLUT.hh"

namespace celeritas {

template<MemSpace M>
class MarchingCubes
{
  public:

    //! Marching cubes with gradient
    MarchingCubes(real_type(* func)(const Real3&),
            Real3(* grad)(const Real3&),
            const int, const int, const int,
            const real_type, const real_type, const real_type);

    size_type operator()(DirectionalMeshData<M>&);

  private:

    std::vector<real_type> vertex_vals_;
    std::vector<size_type> x_verts_, y_verts_, z_verts_;

    const int voxels_x_, voxels_y_, voxels_z_;
    const int pts_x_, pts_y_, pts_z_;
    const real_type length_x_, length_y_, length_z_;
    const real_type step_x_, step_y_, step_z_;

    size_type num_triangles_ = 0;

    real_type(* func_)(const Real3&);
    Real3(* grad_)(const Real3&);

    CollectionBuilder<Triangle>* triangle_builder_;
    CollectionBuilder<Real3>* vert_builder_;
    CollectionBuilder<Real3>* normal_builder_;

    std::vector<Triangle> triangles_;
    std::vector<Real3> verts_;
    std::vector<Real3> normals_;

    void compute_vertex_vals();
    void compute_intersection_points();

    real_type get_data(size_t i, size_t j, size_t k) { return vertex_vals_.at(i + j*pts_x_ + k*pts_x_*pts_y_); }

    void process_cube(Array<real_type, 8> cube, int lut_entry, int i, int j, int k);

    inline Real3 index_to_coord(int i, int j, int k) { return {i * step_x_, j * step_y_, k * step_z_}; }

    inline void set_x_vert(int idx, int i, int j, int k) { x_verts_.at(i + j*pts_x_ + k*pts_x_*pts_y_) = idx; }
    inline void set_y_vert(int idx, int i, int j, int k) { y_verts_.at(i + j*pts_x_ + k*pts_x_*pts_y_) = idx; }
    inline void set_z_vert(int idx, int i, int j, int k) { z_verts_.at(i + j*pts_x_ + k*pts_x_*pts_y_) = idx; }

    inline int get_x_vert(int i, int j, int k) { return x_verts_.at(i + j*pts_x_ + k*pts_x_*pts_y_); }
    inline int get_y_vert(int i, int j, int k) { return y_verts_.at(i + j*pts_x_ + k*pts_x_*pts_y_); }
    inline int get_z_vert(int i, int j, int k) { return z_verts_.at(i + j*pts_x_ + k*pts_x_*pts_y_); }

    inline size_type add_x_vertex(Array<real_type, 8> cube, int i, int j, int k);
    inline size_type add_y_vertex(Array<real_type, 8> cube, int i, int j, int k);
    inline size_type add_z_vertex(Array<real_type, 8> cube, int i, int j, int k);
    inline size_type add_c_vertex(Array<real_type, 8>, int, int, int);

    void add_triangle(const char* trig, char n, int i, int j, int k, int v12 = -1);

    bool test_face(signed char face, Array<real_type, 8> cube);
    bool test_interior(signed char interior, signed char cube_case, signed char config, signed char subconfig, Array<real_type, 8> cube);

};

}
