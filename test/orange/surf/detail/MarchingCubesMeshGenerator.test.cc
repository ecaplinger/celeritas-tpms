//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/surf/detail/MarchingCubesMeshGenerator.test.cc
//---------------------------------------------------------------------------//
#include "orange/surf/detail/MarchingCubesMeshGenerator.hh"
#include "orange/surf/detail/TPMSDefinition.hh"

#include "TestMacros.hh"

#include "celeritas_test.hh"
#include "corecel/Types.hh"
#include "corecel/Constants.hh"

#include <algorithm>
#include <array>
#include <vector>
#include <fstream>
#include <nlohmann/json.hpp>
// #include "MarchingCubesMeshGenerator.test.hh"

namespace celeritas
{
namespace test
{
namespace detail
{

using constants::pi;

namespace {
// algorithm description:
// for each triangle in expected, find a topologically matching triangle in
// result.
void expect_triangles_eq(
        const ::celeritas::detail::MarchingCubesMeshGenerator::Output& expected,
        const ::celeritas::detail::MarchingCubesMeshGenerator::Output& result)
{
    auto is_near = [] (real_type a, real_type b) { return std::fabs(a - b) < 1e-5; };

    for (int i = 0; i < expected.triangles.size(); i++)
    {
        bool found = false;
        Real3 exp_v1 = expected.global_verts[expected.triangles[i][0]];
        Real3 exp_v2 = expected.global_verts[expected.triangles[i][1]];
        Real3 exp_v3 = expected.global_verts[expected.triangles[i][2]];

        Real3 exp_n1 = expected.global_normals[expected.triangles[i][0]];
        Real3 exp_n2 = expected.global_normals[expected.triangles[i][1]];
        Real3 exp_n3 = expected.global_normals[expected.triangles[i][2]];


        for (int j = 0; j < result.triangles.size(); j++)
        {
            // test all winding orders
            for (int w = 0; w < 3; w++)
            {
                bool this_one = true;

                Real3 res_v1 = result.global_verts[result.triangles[j][w]];
                Real3 res_v2 = result.global_verts[result.triangles[j][(w+1)%3]];
                Real3 res_v3 = result.global_verts[result.triangles[j][(w+2)%3]];

                Real3 res_n1 = result.global_normals[result.triangles[j][w]];
                Real3 res_n2 = result.global_normals[result.triangles[j][(w+1)%3]];
                Real3 res_n3 = result.global_normals[result.triangles[j][(w+2)%3]];

                for (int d = 0; d < 3; d++)
                {
                    this_one = this_one && is_near(exp_v1[d], res_v1[d]);
                    this_one = this_one && is_near(exp_v2[d], res_v2[d]);
                    this_one = this_one && is_near(exp_v3[d], res_v3[d]);

                    this_one = is_near(exp_n1[d], res_n1[d]);
                    this_one = this_one && is_near(exp_n2[d], res_n2[d]);
                    this_one = this_one && is_near(exp_n3[d], res_n3[d]);

                    if (this_one) break;
                }

                found = found || this_one;

                if (found) break;
            }

            if (found) break;
        }

        EXPECT_TRUE(found) << "expected triangle at index " << i <<
            " not found in results";

        if (!found) break;
    }
}
}  // namespace

//---------------------------------------------------------------------------//
class SphereIsosurface
{
  public:
    static inline CELER_FUNCTION real_type func(Real3 point, real_type radius)
    {
        return point[0]*point[0] + point[1]*point[1] + point[2]*point[2] - radius*radius;
    }

    static inline CELER_FUNCTION Real3 grad(Real3 point, real_type radius)
    {
        (void) radius;

        return {
            2 * point[0],
            2 * point[1],
            2 * point[2]
        };
    }
};

class MarchingCubesMeshGeneratorTest : public ::celeritas::test::Test
{
  protected:
    // load in reference results from JSON to compare against
    void load_expected(std::string filename)
    {
        std::string filepath = test_data_path("orange", filename);

        std::ifstream infile(filepath);
        CELER_VALIDATE(
            infile, << "failed to open geometry at '" << filename << '\'');

        auto j = nlohmann::json::parse(infile);

        std::vector<std::vector<real_type>> mid_verts =
            j["verts"].get<std::vector<std::vector<real_type>>>();
        std::vector<std::vector<real_type>> mid_normals =
            j["normals"].get<std::vector<std::vector<real_type>>>();
        std::vector<std::vector<int>> mid_triangles =
            j["triangles"].get<std::vector<std::vector<int>>>();

        CELER_EXPECT(mid_verts.size() == mid_normals.size());

        for (int i = 0; i < mid_verts.size(); i++)
        {
            CELER_EXPECT(mid_verts[i].size() == 3);
            CELER_EXPECT(mid_normals[i].size() == 3);
            expected_.global_verts.push_back(
                    {mid_verts[i][0], mid_verts[i][1], mid_verts[i][2]});
            expected_.global_normals.push_back(
                    {mid_normals[i][0], mid_normals[i][1], mid_normals[i][2]});
        }

        for (size_type i = 0; i < mid_triangles.size(); i++)
        {
            CELER_EXPECT(mid_triangles[i].size() == 3);
            expected_.triangles.push_back(
                    {mid_triangles[i][0], mid_triangles[i][1], mid_triangles[i][2]});
        }
    }

    ::celeritas::detail::MarchingCubesMeshGenerator::Output expected_;
};

TEST_F(MarchingCubesMeshGeneratorTest, schwarz_p)
{
    load_expected("mc-schwarz-p.org.json");

    const ::celeritas::detail::MarchingCubesMeshGenerator::Input inp{
        {10, 10, 10},
        {1.0, 1.0, 1.0},
        {0.0, 0.0, 0.0},
        0.0
    };

    auto test_func = [] (Real3 point) { return SchwarzIsosurface::func(point, 1.0); };
    auto test_grad = [] (Real3 point) { return SchwarzIsosurface::grad(point, 1.0); };

    ::celeritas::detail::MarchingCubesMeshGenerator mc(
            test_func, test_grad, inp);
    ::celeritas::detail::MarchingCubesMeshGenerator::Output result = mc();

    EXPECT_EQ(expected_.global_verts.size(), result.global_verts.size());
    EXPECT_EQ(expected_.global_normals.size(), result.global_normals.size());
    EXPECT_EQ(expected_.triangles.size(), result.triangles.size());

    expect_triangles_eq(expected_, result);
    expect_triangles_eq(result, expected_);
}

TEST_F(MarchingCubesMeshGeneratorTest, sphere)
{
    load_expected("mc-sphere.org.json");

    const ::celeritas::detail::MarchingCubesMeshGenerator::Input inp{
        {10, 10, 10},
        {2.0, 2.0, 2.0},
        {-1.0, -1.0, -1.0},
        0.0
    };

    auto test_func = [] (Real3 point) { return SphereIsosurface::func(point, 1.0); };
    auto test_grad = [] (Real3 point) { return SphereIsosurface::grad(point, 1.0); };

    ::celeritas::detail::MarchingCubesMeshGenerator mc(
            test_func, test_grad, inp);
    ::celeritas::detail::MarchingCubesMeshGenerator::Output result = mc();

    EXPECT_EQ(expected_.global_verts.size(), result.global_verts.size());
    EXPECT_EQ(expected_.global_normals.size(), result.global_normals.size());
    EXPECT_EQ(expected_.triangles.size(), result.triangles.size());

    expect_triangles_eq(expected_, result);
    expect_triangles_eq(result, expected_);
}

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace test
}  // namespace celeritas
