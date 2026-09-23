//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file corecel/math/FerrariSolver.hh
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>

#include "corecel/Constants.hh"
#include "corecel/Types.hh"
#include "corecel/cont/Array.hh"
#include "corecel/math/NewtonSolver.hh"
#include "corecel/math/Interval.hh"
#include "orange/surf/detail/TPMSDefinition.hh"
#include "orange/surf/detail/VoxelAcceleratedMeshIntersector.hh"
#include "corecel/data/Collection.hh"

namespace celeritas {

template<class SurfDefinition>
class TPMSUnitCellIntersector {
  public:
    explicit CELER_FUNCTION TPMSUnitCellIntersector(
            DeviceCRef<MeshData> mesh,
            DeviceCRef<VoxelLocatorData> locator,
            const real_type cell_size, const real_type t)
        : mesh_intersector_{mesh, locator}
        , t_{t}
        , p_{cell_size}
        , a_{2*pi / cell_size}
    {
    }

    CELER_FUNCTION real_type operator()(
            const Real3& start, const Real3& direction)
    {
        Array<Real3, MaxIntersections> mesh_normals;
        Array<real_type, MaxIntersections> mesh_distances;

        real_type dmax = distance_to_cell_edge(start, direction);

        // mesh intersector operator() takes mesh_normals and mesh_distances
        // as references, and populates them. but, it may be better to return
        // by value?
        int num_intersections = mesh_intersector_(
                start, direction, mesh_distances, mesh_normals);

        for (int i = 0; i < num_intersections; i++)
        {
            Real3 tri_normal = mesh_normals[i];
            Real3 guess = {
                mesh_distances[i],
                std::acos(mesh_normals[i][2]),
                std::atan2(mesh_normals[i][1], mesh_normals[i][0])
            };

            real_type newton_guess = run_newton(start, direction, guess);
            if (std::isfinite(newton_guess)
                    && is_in_cell(start, direction, newton_guess))
                return newton_guess;
        }

        return dmax;
    }

  private:

    constexpr static int MaxIntersections = SurfDefinition::MaxIntersections;

    CELER_FUNCTION real_type run_newton(const Real3& start,
            const Real3& direction, Real3 guess)
    {
        auto eval_here = [&] (const Real3& point, Real3& p_eval,
                Real9& p_jacobian)
        {
            SurfDefinition::evaluate_point(
                    start, direction, t_, a_,
                    point, p_eval, p_jacobian);
        };

        NewtonSolver3 nsolver{eval_here};

        Real3 result;
        bool success = nsolver(guess, result);

        if (success)
            return result[0];
        return std::numeric_limits<real_type>::infinity();
    }

    CELER_FUNCTION bool is_in_cell(const Real3& start, const Real3& direction,
            const real_type distance)
    {
        const Real3 point = {
            start[0] + distance * direction[0],
            start[1] + distance * direction[1],
            start[2] + distance * direction[2]
        };

        return point[0] >= 0 && point[0] <= p_
            && point[1] >= 0 && point[1] <= p_
            && point[2] >= 0 && point[2] <= p_;
    }

    CELER_FUNCTION real_type distance_to_cell_edge(const Real3& start,
            const Real3& direction)
    {
        real_type distance = std::numeric_limits<real_type>::infinity();
        for (int dim = 0; dim < 3; dim++)
        {
            real_type ddim;
            if (direction[dim] < 0)
                ddim = (0 - start[dim]) / direction[dim];
            else if (direction[dim] > 0)
                ddim = (p_ - start[dim]) / direction[dim];

            distance = std::min(distance, ddim);
        }

        return distance;
    }

    VoxelAcceleratedMeshIntersector<MaxIntersections> mesh_intersector_;

    constexpr static int intersections_expected = 10;

    real_type t_, p_, a_;
};

}
