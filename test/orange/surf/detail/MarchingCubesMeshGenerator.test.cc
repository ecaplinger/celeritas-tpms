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
        0.0, 1.0
    };

    ::celeritas::detail::MarchingCubesMeshGenerator mc(
            &SchwarzIsosurface::func, &SchwarzIsosurface::grad, inp);
    ::celeritas::detail::MarchingCubesMeshGenerator::Output result = mc();

    EXPECT_VEC_NEAR(expected_.global_verts, result.global_verts, 1e-7);
    EXPECT_VEC_NEAR(expected_.global_normals, result.global_normals, 1e-7);
    EXPECT_VEC_EQ(expected_.triangles, result.triangles);
}

TEST_F(MarchingCubesMeshGeneratorTest, sphere)
{
    load_expected("mc-sphere.org.json");

    const ::celeritas::detail::MarchingCubesMeshGenerator::Input inp{
        {10, 10, 10},
        {2.0, 2.0, 2.0},
        {-1.0, -1.0, -1.0},
        0.0, 1.0
    };

    ::celeritas::detail::MarchingCubesMeshGenerator mc(
            &SphereIsosurface::func, &SphereIsosurface::grad, inp);
    ::celeritas::detail::MarchingCubesMeshGenerator::Output result = mc();

    EXPECT_VEC_NEAR(expected_.global_verts, result.global_verts, 1e-7);
    EXPECT_VEC_NEAR(expected_.global_normals, result.global_normals, 1e-7);
    EXPECT_EQ(expected_.triangles.size(), result.triangles.size());

    int matches = 0;
    for (const Array<int, 3>& expected_tri : expected_.triangles)
        for (const Array<int, 3>& result_tri : result.triangles)
            if (expected_tri == result_tri)
            {
                matches++;
                break;
            }

    EXPECT_EQ(matches, result.triangles.size());
}

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace test
}  // namespace celeritas
