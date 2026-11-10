//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/surf/detail/MarchingCubesMeshGenerator.hh
//---------------------------------------------------------------------------//
#pragma once

#include <vector>

#include "corecel/Types.hh"
#include "corecel/cont/Array.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Generates vertices, normals, and connectivity for the isosurface of a
 * scalar field sampled on a regular voxel grid.
 *
 * Adapted from Thomas Lewiner's topologically accurate algorithm:
 * https://thomas.lewiner.org/pub/marching_cubes_jgt.html
 */
class MarchingCubesMeshGenerator
{
  public:
    //!@{
    //! \name Type aliases
    using GeneratingFunctionT = real_type (*)(Real3);
    using GradientFunctionT = Real3 (*)(Real3);
    using VerticesT = std::vector<Real3>;
    using NormalsT = std::vector<Real3>;
    using TriangleT = Array<int, 3>;
    using TrianglesT = std::vector<TriangleT>;
    //!@}

    //! Generator configuration
    struct Input
    {
      Array<int, 3> voxels;
      Real3 length;
      Real3 offset;
      real_type isovalue = 0;
    };

    //! Algorithm output
    struct Output
    {
      VerticesT global_verts;
      NormalsT global_normals;
      TrianglesT triangles;
    };

  public:
    //! Construct with the given configuration
    explicit MarchingCubesMeshGenerator(GeneratingFunctionT, GradientFunctionT, Input const&);

    //! Append the generated isosurface to the given mesh data
    Output operator()();

  private:
    //!@{
    //! \name Mesh generation
    void compute_vertex_vals();
    void compute_intersection_points();
    void process_cube(Array<real_type, 8> const& cube, int lut_entry,
                      int i, int j, int k);
    void add_triangle(const char* trig, char n, int i, int j, int k,
                      int v12 = -1);
    bool test_face(signed char face, Array<real_type, 8> const& cube) const;
    bool test_interior(signed char interior, signed char cube_case,
                       signed char config, signed char subconfig,
                       Array<real_type, 8> const& cube) const;
    //!@}

    //!@{
    //! \name Vertex insertion
    int add_x_vertex(Array<real_type, 8> const& cube, int i, int j, int k);
    int add_y_vertex(Array<real_type, 8> const& cube, int i, int j, int k);
    int add_z_vertex(Array<real_type, 8> const& cube, int i, int j, int k);
    int add_c_vertex(Array<real_type, 8> const& cube, int i, int j, int k);
    //!@}

    //!@{
    //! \name Grid access
    real_type get_data(int i, int j, int k) const;
    Real3 index_to_coord(int i, int j, int k) const;
    void set_x_vert(int idx, int i, int j, int k);
    void set_y_vert(int idx, int i, int j, int k);
    void set_z_vert(int idx, int i, int j, int k);
    int get_x_vert(int i, int j, int k) const;
    int get_y_vert(int i, int j, int k) const;
    int get_z_vert(int i, int j, int k) const;
    //!@}

  private:
    const int voxels_x_, voxels_y_, voxels_z_;
    const int pts_x_, pts_y_, pts_z_;
    const real_type length_x_, length_y_, length_z_;
    const real_type offset_x_, offset_y_, offset_z_;
    const real_type step_x_, step_y_, step_z_;
    const real_type iso_;

    std::vector<real_type> vertex_vals_;
    std::vector<int> x_verts_, y_verts_, z_verts_;

    GeneratingFunctionT func_;
    GradientFunctionT grad_;

    std::vector<TriangleT> triangles_;
    std::vector<Real3> verts_;
    std::vector<Real3> normals_;

    //! Per-cube debug records
    std::vector<std::string> dbg_lines_;
    int dbg_cube_ = -1;
};

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace celeritas
