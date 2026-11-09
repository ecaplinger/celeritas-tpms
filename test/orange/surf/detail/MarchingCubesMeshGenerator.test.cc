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
    void SetUp() override {}

    void load_expected(std::string filename)
    {
        std::string filepath = test_data_path("orange", filename);

        std::ifstream infile(filename);
        CELER_VALIDATE(
            infile, << "failed to open geometry at '" << filename << '\'');

        auto j = nlohmann::json::parse(infile);
        expected.global_verts = j["verts"].get<std::vector<Real3>>();
        expected.global_normals = j["normals"].get<std::vector<Real3>>();
        expected.triangles = j["triangles"].get<std::vector<Array<int, 3>>>();
    }

    ::celeritas::detail::MarchingCubesMeshGenerator::Output expected;
};

TEST_F(MarchingCubesMeshGeneratorTest, schwarz_p)
{
    load_expected("mc_schwarz_p.org.json");

    const ::celeritas::detail::MarchingCubesMeshGenerator::Input inp{
        {100, 100, 100},
        {1.0, 1.0, 1.0},
        0.0, 1.0
    };

    ::celeritas::detail::MarchingCubesMeshGenerator mc(
            &SchwarzSurface::func, &SchwarzSurface::grad, inp);
    ::celeritas::detail::MarchingCubesMeshGenerator::Output result = mc();

    EXPECT_VEC_SOFT_EQ(expected.global_normals, result.global_normals);
    EXPECT_VEC_SOFT_EQ(expected.global_verts, result.global_normals);
    EXPECT_VEC_SOFT_EQ(expected.triangles, result.triangles);
}

TEST_F(MarchingCubesMeshGeneratorTest, sphere)
{
    load_expected("mc_sphere.org.json");

    const ::celeritas::detail::MarchingCubesMeshGenerator::Input inp{
        {100, 100, 100},
        {1.0, 1.0, 1.0},
        0.0, 1.0
    };

    ::celeritas::detail::MarchingCubesMeshGenerator mc(
            &SphereIsosurface::func, &SphereIsosurface::grad, inp);
    ::celeritas::detail::MarchingCubesMeshGenerator::Output result = mc();

    EXPECT_VEC_SOFT_EQ(expected.global_normals, result.global_normals);
    EXPECT_VEC_SOFT_EQ(expected.global_verts, result.global_normals);
    EXPECT_VEC_SOFT_EQ(expected.triangles, result.triangles);
}

// TEST_F(MarchingCubesMeshGeneratorTest, TEST_IF_CELER_DEVICE(device))
// {
//     MCMGTestInput input;
//     mcmg_test(input);
// }

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace test
}  // namespace celeritas
