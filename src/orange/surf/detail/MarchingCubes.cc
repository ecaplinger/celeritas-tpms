//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/detail/BvhBuilder.cc
//---------------------------------------------------------------------------//
#include "MarchingCubes.hh"

#include "corecel/data/Collection.hh"
#include "corecel/math/ArrayUtils.hh"

#include "orange/surf/detail/MarchingCubesLUT.hh"

namespace celeritas
{

template<MemSpace M>
size_type MarchingCubes<M>::add_x_vertex(Array<real_type, 8> cube, int i, int j, int k)
{
    Real3 coord = index_to_coord(i, j, k);

    real_type u = cube[0] / (cube[0] - cube[1]);
    coord[0] += u;

    size_type r = verts_.size();

    verts_.push_back(coord);

    Real3 normal = grad_(coord);
    normal = make_unit_vector(normal);
    normals_.push_back(normal);

    return r;
}

template<MemSpace M>
size_type MarchingCubes<M>::add_y_vertex(Array<real_type, 8> cube, int i, int j, int k)
{
    Real3 coord = index_to_coord(i, j, k);

    real_type u = cube[0] / (cube[0] - cube[3]);
    coord[1] += u;

    size_type r = verts_.size();

    verts_.push_back(coord);

    Real3 normal = grad_(coord);
    normal = make_unit_vector(normal);
    normals_.push_back(normal);

    return r;
}

template<MemSpace M>
size_type MarchingCubes<M>::add_z_vertex(Array<real_type, 8> cube, int i, int j, int k)
{
    Real3 coord = index_to_coord(i, j, k);

    real_type u = cube[0] / (cube[0] - cube[4]);
    coord[2] += u;

    size_type r = verts_.size();

    verts_.push_back(coord);

    Real3 normal = grad_(coord);
    normal = make_unit_vector(normal);
    normals_.push_back(normal);

    return r;
}

template<MemSpace M>
size_type MarchingCubes<M>::add_c_vertex(Array<real_type, 8> cube, int i, int j, int k)
{
    real_type u = 0;
    size_type r = verts_.size();

    Real3 vertex{0, 0, 0};

    auto add_vert = [&] (int vid)
    {
        if (vid == -1) return;

        u += 1;
        const Real3& vx = verts_.at(vid);
        vertex[0] += vx[0];
        vertex[1] += vx[1];
        vertex[2] += vx[2];
    };

    add_vert(get_x_vert(i, j, k));
    add_vert(get_y_vert(i+1, j, k));
    add_vert(get_x_vert(i, j+1, k));
    add_vert(get_y_vert(i, j, k));
    add_vert(get_x_vert(i, j, k+1));
    add_vert(get_y_vert(i+1, j, k+1));
    add_vert(get_x_vert(i, j+1, k+1));
    add_vert(get_y_vert(i, j, k+1));
    add_vert(get_z_vert(i, j, k));
    add_vert(get_z_vert(i+1, j, k));
    add_vert(get_z_vert(i+1, j+1, k));
    add_vert(get_z_vert(i, j+1, k));

    vertex[0] /= u;
    vertex[1] /= u;
    vertex[2] /= u;

    verts_.push_back(vertex);

    Real3 normal = grad_(vertex);
    normal = make_unit_vector(normal);
    normals_.push_back(normal);

    return r;
}

template<MemSpace M>
MarchingCubes<M>::MarchingCubes(real_type(* func)(const Real3&),
        Real3(* grad)(const Real3&),
        const int voxels_x, const int voxels_y, const int voxels_z,
        const real_type length_x,
        const real_type length_y,
        const real_type length_z)
    : func_{func}
    , grad_{grad}
    , voxels_x_{voxels_x}
    , voxels_y_{voxels_y}
    , voxels_z_{voxels_z}
    , length_x_{length_x}
    , length_y_{length_y}
    , length_z_{length_z}
    , step_x_{length_x / voxels_x}
    , step_y_{length_y / voxels_y}
    , step_z_{length_z / voxels_z}
    , pts_x_{voxels_x + 1}
    , pts_y_{voxels_y + 1}
    , pts_z_{voxels_z + 1}
    , x_verts_(length_x * length_y * length_z, -1)
    , y_verts_(length_x * length_y * length_z, -1)
    , z_verts_(length_x * length_y * length_z, -1)
{
}

template<MemSpace M>
void MarchingCubes<M>::compute_vertex_vals()
{
    vertex_vals_.reserve(pts_x_ * pts_y_ * pts_z_);

    for (int k = 0; k < pts_z_; k++)
        for (int j = 0; j < pts_y_; j++)
            for (int i = 0; i < pts_x_; i++)
                vertex_vals_[i + j*pts_x_ + k*pts_x_*pts_y_] =
                    func_(index_to_coord(i, j, k));
}

template<MemSpace M>
void MarchingCubes<M>::compute_intersection_points()
{
    for (int i = 0; i < pts_x_; i++)
        for (int j = 0; j < pts_y_; j++)
            for (int k = 0; k < pts_z_; k++)
            {
                Array<real_type, 8> cube;
                cube[0] = get_data(i, j, k);

                if (i < voxels_x_)
                    cube[1] = get_data(i+1, j, k);
                else
                    cube[1] = cube[0];

                if (j < voxels_y_)
                    cube[3] = get_data(i, j+1, k);
                else
                    cube[3] = cube[0];

                if (k < voxels_z_)
                    cube[4] = get_data(i, j, k+1);
                else
                    cube[4] = cube[0];

                if (cube[0] < 0)
                {
                    if (cube[1] > 0) set_x_vert(add_x_vertex(cube, i, j, k), i, j, k);
                    if (cube[3] > 0) set_y_vert(add_y_vertex(cube, i, j, k), i, j, k);
                    if (cube[4] > 0) set_z_vert(add_z_vertex(cube, i, j, k), i, j, k);
                }
                else
                {
                    if (cube[1] < 0) set_x_vert(add_x_vertex(cube, i, j, k), i, j, k);
                    if (cube[3] < 0) set_y_vert(add_y_vertex(cube, i, j, k), i, j, k);
                    if (cube[4] < 0) set_z_vert(add_z_vertex(cube, i, j, k), i, j, k);
                }
            }
}

template<MemSpace M>
size_type MarchingCubes<M>::operator()(DirectionalMeshData<M>& data)
{
    compute_vertex_vals();

    for (int i = 0; i < voxels_x_; i++)
        for (int j = 0; j < voxels_y_; j++)
            for (int k = 0; k < voxels_z_; k++)
            {
                Array<real_type, 8> cube;
                int lut_entry = 0;
                for (int p = 0; p < 8; p++)
                {
                    cube[p] = vertex_vals_[i + (p&1)
                        + (j + ((p>>1)&1))*(voxels_x_+1)
                        + (k + ((p>>2)&1))*(voxels_y_+1)*(voxels_x_+1)];
                    if (cube[p] > 0) lut_entry += 1 << p;
                }

                process_cube(cube, lut_entry, i, j, k);
            }

    CollectionBuilder{&data.global_verts}.insert_back(
            verts_.begin(), verts_.end());
    CollectionBuilder{&data.global_normals}.insert_back(
            normals_.begin(), normals_.end());
    CollectionBuilder{&data.triangles}.insert_back(
            triangles_.begin(), triangles_.end());

    return triangles_.size();
}

template<MemSpace M>
void MarchingCubes<M>::process_cube(Array<real_type, 8> cube, int lut_entry, int i, int j, int k)
{
    char cube_case = cases[lut_entry][0];
    char config = cases[lut_entry][1];
    char subconfig = 0;

    int v12;

    switch (cube_case)
    {
        case 0:
            break;

        case 1:
            add_triangle(tiling1[config], 1, i, j, k);
            break;

        case 2:
            add_triangle(tiling2[config], 2, i, j, k);

        case 3:
            if (test_face(test3[config], cube))
                add_triangle(tiling3_2[config], 4, i, j, k);
            else
                add_triangle(tiling3_1[config], 2, i, j, k);
            break;

        case 4:
            if (test_interior(test4[config], cube_case, config, subconfig, cube))
                add_triangle(tiling4_1[config], 2, i, j, k);
            else
                add_triangle(tiling4_2[config], 6, i, j, k);
            break;

        case 5:
            add_triangle(tiling5[config], 3, i, j, k);
            break;

        case 6:
            if (test_face(test6[config][0], cube))
                add_triangle(tiling6_2[config], 5, i, j, k);
            else
            {
                if (test_interior(test6[config][1], cube_case, config, subconfig, cube))
                    add_triangle(tiling6_1_1[config], 3, i, j, k);
                else
                {
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling6_1_2[config], 9, i, j, k, v12);
                }
            }

        case 7:
            if (test_face(test7[config][0], cube)) subconfig += 1;
            if (test_face(test7[config][1], cube)) subconfig += 2;
            if (test_face(test7[config][2], cube)) subconfig += 4;
            switch (subconfig)
            {
                case 0:
                    add_triangle(tiling7_1[config], 3, i, j, k); break;
                case 1:
                    add_triangle(tiling7_2[config][0], 5, i, j, k); break;
                case 2:
                    add_triangle(tiling7_2[config][1], 5, i, j, k); break;
                case 3:
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling7_3[config][0], 9, i, j, k, v12); break;
                case 4 :
                    add_triangle(tiling7_2[config][2], 5, i, j, k); break;
                case 5 :
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling7_3[config][1], 9, i, j, k, v12); break;
                case 6 :
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling7_3[config][2], 9, i, j, k, v12); break;
                case 7 :
                    if (test_interior(test7[config][3], cube_case, config, subconfig, cube))
                        add_triangle(tiling7_4_2[config], 9, i, j, k);
                    else
                        add_triangle(tiling7_4_1[config], 5, i, j, k);
                    break;
            }
            break;

        case 8:
            add_triangle(tiling8[config], 2, i, j, k);
            break;

        case 9:
            add_triangle(tiling9[config], 4, i, j, k);
            break;

        case 10:
            if (test_face(test10[config][0], cube))
            {
                if (test_face(test10[config][1], cube))
                    add_triangle(tiling10_1_1_[config], 4, i, j, k); // 10.1.1
                else
                {
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling10_2[config], 8, i, j, k, v12); // 10.2
                }
            }
            else
            {
                if (test_face(test10[config][1], cube))
                {
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling10_2_[config], 8, i, j, k, v12); // 10.2
                }
                else
                {
                    if (test_interior(test10[config][2], cube_case, config, subconfig, cube))
                        add_triangle(tiling10_1_1[config], 4, i, j, k); // 10.1.1
                    else
                        add_triangle(tiling10_1_2[config], 8, i, j, k); // 10.1.2
                }
            }
            break;

        case 11 :
            add_triangle(tiling11[config], 4, i, j, k);
            break ;

        case 12 :
            if (test_face(test12[config][0], cube))
            {
                if (test_face(test12[config][1], cube))
                    add_triangle(tiling12_1_1_[config], 4, i, j, k); // 12.1.1
                else
                {
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling12_2[config], 8, i, j, k, v12); // 12.2
                }
            }
            else
            {
                if (test_face(test12[config][1], cube))
                {
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling12_2_[config], 8, i, j, k, v12); // 12.2
                }
                else
                {
                    if (test_interior(test12[config][2], cube_case, config, subconfig, cube))
                        add_triangle(tiling12_1_1[config], 4, i, j, k); // 12.1.1
                    else
                        add_triangle(tiling12_1_2[config], 8, i, j, k); // 12.1.2
                }
            }
            break ;

        case 13 :
            if (test_face(test13[config][0], cube)) subconfig += 1;
            if (test_face(test13[config][1], cube)) subconfig += 2;
            if (test_face(test13[config][2], cube)) subconfig += 4;
            if (test_face(test13[config][3], cube)) subconfig += 8;
            if (test_face(test13[config][4], cube)) subconfig += 16;
            if (test_face(test13[config][5], cube)) subconfig += 32;
            switch (subconfig13[subconfig])
            {
                case 0:/* 13.1 */
                    add_triangle(tiling13_1[config], 4, i, j, k); break;

                case 1:/* 13.2 */
                    add_triangle(tiling13_2[config][0], 6, i, j, k); break;
                case 2:/* 13.2 */
                    add_triangle(tiling13_2[config][1], 6, i, j, k); break;
                case 3:/* 13.2 */
                    add_triangle(tiling13_2[config][2], 6, i, j, k); break;
                case 4:/* 13.2 */
                    add_triangle(tiling13_2[config][3], 6, i, j, k); break;
                case 5:/* 13.2 */
                    add_triangle(tiling13_2[config][4], 6, i, j, k); break;
                case 6:/* 13.2 */
                    add_triangle(tiling13_2[config][5], 6, i, j, k); break;

                case 7:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][0], 10, i, j, k, v12); break;
                case 8:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][1], 10, i, j, k, v12); break;
                case 9:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][2], 10, i, j, k, v12); break;
                case 10:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][3], 10, i, j, k, v12); break;
                case 11:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][4], 10, i, j, k, v12); break;
                case 12:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][5], 10, i, j, k, v12); break;
                case 13:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][6], 10, i, j, k, v12); break;
                case 14:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][7], 10, i, j, k, v12); break;
                case 15:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][8], 10, i, j, k, v12); break;
                case 16:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][9], 10, i, j, k, v12); break;
                case 17:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][10], 10, i, j, k, v12); break;
                case 18:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3[config][11], 10, i, j, k, v12); break;

                case 19:/* 13.4 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_4[config][0], 12, i, j, k, v12); break;
                case 20:/* 13.4 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_4[config][1], 12, i, j, k, v12); break;
                case 21:/* 13.4 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_4[config][2], 12, i, j, k, v12); break;
                case 22:/* 13.4 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_4[config][3], 12, i, j, k, v12); break;

                case 23:/* 13.5 */
                    subconfig = 0 ;
                    if (test_interior(test13[config][6], cube_case, config, subconfig, cube))
                        add_triangle(tiling13_5_1[config][0], 6, i, j, k);
                    else
                        add_triangle(tiling13_5_2[config][0], 10, i, j, k);
                    break ;
                case 24:/* 13.5 */
                    subconfig = 1 ;
                    if (test_interior(test13[config][6], cube_case, config, subconfig, cube))
                        add_triangle(tiling13_5_1[config][1], 6, i, j, k);
                    else
                        add_triangle(tiling13_5_2[config][1], 10, i, j, k);
                    break ;
                case 25:/* 13.5 */
                    subconfig = 2 ;
                    if (test_interior(test13[config][6], cube_case, config, subconfig, cube))
                        add_triangle(tiling13_5_1[config][2], 6, i, j, k);
                    else
                        add_triangle(tiling13_5_2[config][2], 10, i, j, k);
                    break ;
                case 26:/* 13.5 */
                    subconfig = 3 ;
                    if (test_interior(test13[config][6], cube_case, config, subconfig, cube))
                        add_triangle(tiling13_5_1[config][3], 6, i, j, k);
                    else
                        add_triangle(tiling13_5_2[config][3], 10, i, j, k);
                    break ;

                case 27:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][0], 10, i, j, k, v12); break;
                case 28:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][1], 10, i, j, k, v12); break;
                case 29:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][2], 10, i, j, k, v12); break;
                case 30:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][3], 10, i, j, k, v12); break;
                case 31:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][4], 10, i, j, k, v12); break;
                case 32:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][5], 10, i, j, k, v12); break;
                case 33:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][6], 10, i, j, k, v12); break;
                case 34:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][7], 10, i, j, k, v12); break;
                case 35:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k);
                    add_triangle(tiling13_3_[config][8], 10, i, j, k, v12); break;
                case 36:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k) ;
                    add_triangle(tiling13_3_[config][9], 10, i, j, k, v12); break;
                case 37:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k) ;
                    add_triangle(tiling13_3_[config][10], 10, i, j, k, v12); break;
                case 38:/* 13.3 */
                    v12 = add_c_vertex(cube, i, j, k) ;
                    add_triangle(tiling13_3_[config][11], 10, i, j, k, v12); break;

                case 39:/* 13.2 */
                    add_triangle(tiling13_2_[config][0], 6, i, j, k); break;
                case 40:/* 13.2 */
                    add_triangle(tiling13_2_[config][1], 6, i, j, k); break;
                case 41:/* 13.2 */
                    add_triangle(tiling13_2_[config][2], 6, i, j, k); break;
                case 42:/* 13.2 */
                    add_triangle(tiling13_2_[config][3], 6, i, j, k); break;
                case 43:/* 13.2 */
                    add_triangle(tiling13_2_[config][4], 6, i, j, k); break;
                case 44:/* 13.2 */
                    add_triangle(tiling13_2_[config][5], 6, i, j, k); break;

                case 45:/* 13.1 */
                    add_triangle(tiling13_1_[config], 4, i, j, k); break;

                default:
                    break;
                }
                break ;

        case 14 :
            add_triangle(tiling14[config], 4, i, j, k);
            break ;

        default:
            break;
    }
}

template<MemSpace M>
void MarchingCubes<M>::add_triangle(const char* trig, char n, int i, int j, int k, int v12)
{
    int tv[3];

    for (int t = 0; t < 3*n; t++)
    {
        switch (trig[t])
        {
            case 0:
                tv[t%3] = get_x_vert(i, j, k); break;
            case 1:
                tv[t%3] = get_y_vert(i+1, j, k); break;
            case 2:
                tv[t%3] = get_x_vert(i, j+1, k); break;
            case 3:
                tv[t%3] = get_y_vert(i, j, k); break;
            case 4:
                tv[t%3] = get_x_vert(i, j, k+1); break;
            case 5:
                tv[t%3] = get_y_vert(i+1, j, k+1); break;
            case 6:
                tv[t%3] = get_x_vert(i, j+1, k+1); break;
            case 7:
                tv[t%3] = get_y_vert(i, j, k+1); break;
            case 8:
                tv[t%3] = get_z_vert(i, j, k); break;
            case 9:
                tv[t%3] = get_z_vert(i+1, j, k); break;
            case 10:
                tv[t%3] = get_z_vert(i+1, j+1, k); break;
            case 11:
                tv[t%3] = get_z_vert(i, j+1, k); break;
            case 12:
                tv[t%3] = v12; break;
            default:
                break;
        }

        if (tv[t%3] == -1)
        {
            printf("Marching Cubes: invalid triangle %d\n", triangles_.size());
        }

        if (t%3 == 2)
        {
            triangles_.push_back(Triangle(tv));
        }
    }
}

template<MemSpace M>
bool MarchingCubes<M>::test_face(signed char face, Array<real_type, 8> cube)
{
    real_type A, B, C, D;

    switch (face)
    {
        case -1: case 1:
            A = cube[0]; B = cube[4]; C = cube[5]; D = cube[1]; break;
        case -2: case 2:
            A = cube[1]; B = cube[5]; C = cube[6]; D = cube[2]; break;
        case -3: case 3:
            A = cube[2]; B = cube[6]; C = cube[7]; D = cube[3]; break;
        case -4: case 4:
            A = cube[3]; B = cube[7]; C = cube[4]; D = cube[0]; break;
        case -5: case 5:
            A = cube[0]; B = cube[3]; C = cube[2]; D = cube[1]; break;
        case -6: case 6:
            A = cube[4]; B = cube[7]; C = cube[6]; D = cube[5]; break;
        default:
            printf( "Invalid face code %d\n", face ); A = B = C = D = 0;
    };

    return face * A * (A*C - B*D) >= 0;
}

template<MemSpace M>
bool MarchingCubes<M>::test_interior(signed char interior,
        signed char cube_case,
        signed char config,
        signed char subconfig,
        Array<real_type, 8> cube)
{
    real_type t, At = 0, Bt = 0, Ct = 0, Dt = 0, a, b;
    char test = 0;
    char edge = -1;
    switch (cube_case)
    {
        case 4:
        case 10:
            a = (cube[4] - cube[0]) * (cube[6] - cube[2]) - (cube[7] - cube[3]) * (cube[5] - cube[1]);
            b = cube[2] * (cube[4] - cube[0]) + cube[0] * (cube[6] - cube[2]) - cube[1] * (cube[7] -cube[3]) - cube[3] * (cube[5] - cube[1]);
            t = -b / (2*a);
            if (t < 0 || t > 1) return interior > 0;

            At = cube[0] + (cube[4] - cube[0]) * t;
            Bt = cube[3] + (cube[7] - cube[3]) * t;
            Ct = cube[2] + (cube[6] - cube[2]) * t;
            Dt = cube[1] + (cube[5] - cube[1]) * t;
            break;

        case 6:
        case 7:
        case 12:
        case 13:
            switch (cube_case)
            {
                case 6: edge = test6[config][2]; break;
                case 7: edge = test7[config][4]; break;
                case 12: edge = test12[config][3]; break;
                case 13: edge = tiling13_5_1[config][subconfig][0]; break;
            };

            switch (edge)
            {
                case 0 :
                    t = cube[0] / (cube[0] - cube[1]);
                    At = 0;
                    Bt = cube[3] + (cube[2] - cube[3]) * t;
                    Ct = cube[7] + (cube[6] - cube[7]) * t;
                    Dt = cube[4] + (cube[5] - cube[4]) * t;
                    break;
                case 1:
                    t = cube[1] / (cube[1] - cube[2]);
                    At = 0;
                    Bt = cube[0] + (cube[3] - cube[0]) * t;
                    Ct = cube[4] + (cube[7] - cube[4]) * t;
                    Dt = cube[5] + (cube[6] - cube[5]) * t;
                    break;
                case 2:
                    t = cube[2] / (cube[2] - cube[3]);
                    At = 0;
                    Bt = cube[1] + (cube[0] - cube[1]) * t;
                    Ct = cube[5] + (cube[4] - cube[5]) * t;
                    Dt = cube[6] + (cube[7] - cube[6]) * t;
                    break;
                case 3:
                    t = cube[3] / (cube[3] - cube[0]);
                    At = 0;
                    Bt = cube[2] + (cube[1] - cube[2]) * t;
                    Ct = cube[6] + (cube[5] - cube[6]) * t;
                    Dt = cube[7] + (cube[4] - cube[7]) * t;
                    break;
                case 4:
                    t = cube[4] / (cube[4] - cube[5]);
                    At = 0;
                    Bt = cube[7] + (cube[6] - cube[7]) * t;
                    Ct = cube[3] + (cube[2] - cube[3]) * t;
                    Dt = cube[0] + (cube[1] - cube[0]) * t;
                    break;
                case 5:
                    t = cube[5] / (cube[5] - cube[6]);
                    At = 0;
                    Bt = cube[4] + (cube[7] - cube[4]) * t;
                    Ct = cube[0] + (cube[3] - cube[0]) * t;
                    Dt = cube[1] + (cube[2] - cube[1]) * t;
                    break;
                case 6:
                    t = cube[6] / (cube[6] - cube[7]);
                    At = 0;
                    Bt = cube[5] + (cube[4] - cube[5]) * t;
                    Ct = cube[1] + (cube[0] - cube[1]) * t;
                    Dt = cube[2] + (cube[3] - cube[2]) * t;
                    break;
                case 7:
                    t = cube[7] / (cube[7] - cube[4]);
                    At = 0;
                    Bt = cube[6] + (cube[5] - cube[6]) * t;
                    Ct = cube[2] + (cube[1] - cube[2]) * t;
                    Dt = cube[3] + (cube[0] - cube[3]) * t;
                    break;
                case 8:
                    t = cube[0] / (cube[0] - cube[4]);
                    At = 0;
                    Bt = cube[3] + (cube[7] - cube[3]) * t;
                    Ct = cube[2] + (cube[6] - cube[2]) * t;
                    Dt = cube[1] + (cube[5] - cube[1]) * t;
                    break;
                case 9:
                    t = cube[1] / (cube[1] - cube[5]);
                    At = 0;
                    Bt = cube[0] + (cube[4] - cube[0]) * t;
                    Ct = cube[3] + (cube[7] - cube[3]) * t;
                    Dt = cube[2] + (cube[6] - cube[2]) * t;
                    break;
                case 10:
                    t = cube[2] / (cube[2] - cube[6]);
                    At = 0;
                    Bt = cube[1] + (cube[5] - cube[1]) * t;
                    Ct = cube[0] + (cube[4] - cube[0]) * t;
                    Dt = cube[3] + (cube[7] - cube[3]) * t;
                    break;
                case 11:
                    t = cube[3] / (cube[3] - cube[7]);
                    At = 0;
                    Bt = cube[2] + (cube[6] - cube[2]) * t;
                    Ct = cube[1] + (cube[5] - cube[1]) * t;
                    Dt = cube[0] + (cube[4] - cube[0]) * t;
                    break;
                default : printf("Invalid edge %d\n", edge); break;
            }
    }

    if (At >= 0) test += 1;
    if (Bt >= 0) test += 2;
    if (Ct >= 0) test += 4;
    if (Dt >= 0) test += 8;

    switch (test)
    {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 6:
        case 8:
        case 9:
        case 12:
            return interior > 0;
        case 5:
            if (At * Ct - Bt * Dt <  FLT_EPSILON)
                return interior > 0;
        case 10:
            if (At * Ct - Bt * Dt >= FLT_EPSILON)
                return interior > 0;
        case 7:
        case 11:
        case 13:
        case 14:
        case 15:
            return interior < 0 ;
    }

    return interior < 0;
}

}
