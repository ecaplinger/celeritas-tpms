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

namespace celeritas {

using constants::pi;
using constants::sqrt_two;
using constants::sqrt_three;

class SchwarzIsosurface
{
  public:
    static inline real_type func(Real3 point, real_type period)
    {
        real_type a = 2 * pi / period;
        return std::cos(a * point[0])
            + std::cos(a * point[1])
            + std::cos(a * point[2]);
    }

    static inline Real3 grad(Real3 point, real_type period)
    {
        real_type a = 2 * pi / period;
        return {
            -a * std::sin(a * point[0]),
            -a * std::sin(a * point[1]),
            -a * std::sin(a * point[2])
        };
    }

  private:
    real_type a_;
};

}

