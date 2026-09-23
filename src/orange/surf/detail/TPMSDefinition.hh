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

namespace celeritas {

using Real9 = Array<real_type, 9>;
using Interval3 = Array<Interval, 3>;
using Interval9 = Array<Interval, 9>;

using constants::pi;
using constants::sqrt_two;
using constants::sqrt_three;

struct TPMSPoint {
    Real3 point_eval;
    Real9 point_jacobian;
};

struct TPMSInterval {
    Real3 point_eval;
    Real9 point_jacobian;
    Real9 lischitz_limits;
    Interval9 interval_jacobian;
};

template<class T>
struct TPMSParams
{
    using T3 = Array<T, 3>;

    CELER_FUNCTION TPMSParams(
            const Real3& start,
            const Real3& direction,
            real_type t,
            real_type a)
        : x0(start[0]), y0(start[1]), z0(start[2]),
        ux(direction[0]), uy(direction[1]), uz(direction[2]),
        t(t), a(a)
    {
        set_params({0, 0, 0});
        static_assert(
                std::is_same_v<T, real_type> || std::is_same_v<T, Interval>,
                "TPMSParams only supports real_type and Interval.");
    }

    inline CELER_FUNCTION void set_params(const T3& params)
    {
        d = params[0];
        th = params[1];
        ph = params[2];

        if constexpr (std::is_same_v<T, real_type>)
        {
            sth = std::sin(th);
            cth = std::cos(th);
            sph = std::sin(ph);
            cph = std::cos(ph);
        }
        else if constexpr (std::is_same_v<T, Interval>)
        {
            sth = Interval::isin(th);
            cth = Interval::icos(th);
            sph = Interval::isin(ph);
            cph = Interval::icos(ph);
        }

        X = a * (x0 + ux*d + t*sth*cph);
        Y = a * (y0 + uy*d + t*sth*sph);
        Z = a * (z0 + uz*d + t*cth);

        if constexpr (std::is_same_v<T, real_type>)
        {
            cx = std::cos(X);
            cy = std::cos(Y);
            cz = std::cos(Z);
            sx = std::sin(X);
            sy = std::sin(Y);
            sz = std::sin(Z);
        }
        else if constexpr (std::is_same_v<T, Interval>)
        {
            cx = Interval::icos(X);
            cy = Interval::icos(Y);
            cz = Interval::icos(Z);
            sx = Interval::isin(X);
            sy = Interval::isin(Y);
            sz = Interval::isin(Z);
        }

        const T cc = a*t*cth*cph;
        const T sc = a*t*sth*cph;
        const T cs = a*t*cth*sph;
        const T ss = a*t*sth*sph;

        dX = {a*ux, cc, -ss};
        dY = {a*uy, cs, sc};
        dZ = {a*uz, -t*sth, 0};

        dXth = {0, -sc, -cs};
        dYth = {0, -ss, cc};
        dZth = {0, -t*cth, 0};

        dXph = {0, -cs, -sc};
        dYph = {0, cc, -ss};
        dZph = {0, 0, 0};
    }

    real_type x0, y0, z0;
    real_type ux, uy, uz;
    real_type t, a;

    T d, th, ph;
    T sth, cth, sph, cph;

    T X, Y, Z;
    T cx, cy, cz, sx, sy, sz;

    T3 dX, dY, dZ;
    T3 dXth, dYth, dZth;
    T3 dXph, dYph, dZph;
};

class SchwarzDefinition
{
  public:

    constexpr static int MaxIntersections = 4;

    CELER_FUNCTION SchwarzDefinition(
            Real3 start, Real3 direction, real_type t, real_type p)
        : start_{start}
        , direction_{direction}
        , t_{t}
        , a_{2 * pi / p}
    {
        generate_lipschitz();
    }

    CELER_FUNCTION inline real_type evaluate_cartesian(
            const Real3& point) const
    {
        return std::cos(a_ * point[0])
            + std::cos(a_ * point[1])
            + std::cos(a_ + point[2]);
    }

    CELER_FUNCTION inline Real3 grad_cartesian(const Real3& point) const
    {
        return {
            -a_ * std::sin(a_ * point[0]),
            -a_ * std::sin(a_ * point[1]),
            -a_ * std::sin(a_ * point[2])
        };
    }

    CELER_FUNCTION const Real9& lipschitz_first() const { return lipschitz_first_; }

    CELER_FUNCTION const Real9& lipschitz_second() const { return lipschitz_second_; }

    CELER_FUNCTION inline void evaluate_point(
            const Real3& point, Real3& p_eval, Real9& p_jacobian) const
    {
        TPMSParams<real_type> params{start_, direction_, t_, a_};
        params.set_params(point);

        const real_type d1 = -params.sx;
        const real_type d2 = -params.sy;
        const real_type d3 = -params.sz;
        const real_type dd1 = -params.cx;
        const real_type dd2 = -params.cy;
        const real_type dd3 = -params.cz;

        p_jacobian[0] = params.dX[0] * d1
            + params.dY[0] * d2
            + params.dZ[0] * d3;

        p_jacobian[1] = params.dX[1] * d1
            + params.dY[1] * d2
            + params.dZ[1] * d3;

        p_jacobian[2] = params.dX[2] * d1
            + params.dY[2] * d2;

        p_jacobian[3] = params.dX[1] * params.dX[0] * dd1
            + params.dY[1] * params.dY[0] * dd2
            + params.dZ[1] * params.dZ[0] * dd3;

        p_jacobian[4] = params.dX[1] * params.dX[1] * dd1
            + params.dY[1] * params.dY[1] * dd2
            + params.dZ[1] * params.dZ[1] * dd3
            + params.dXth[1] * d1
            + params.dYth[1] * d2
            + params.dZth[1] * d3;

        p_jacobian[5] = params.dX[1] * params.dX[2] * dd1
            + params.dY[1] * params.dY[2] * dd2
            + params.dXth[2] * d1
            + params.dYth[2] * d2;

        p_jacobian[6] = params.dX[2] * params.dX[0] * dd1
            + params.dY[2] * params.dY[0] * dd2;

        p_jacobian[7] = p_jacobian[5];

        p_jacobian[8] = params.dX[2] * params.dX[2] * dd1
            + params.dY[2] * params.dY[2] * dd2
            + params.dXph[2] * d1
            + params.dYph[2] * d2;

        p_eval[0] = params.cx + params.cy + params.cz;
        p_eval[1] = p_jacobian[1];
        p_eval[2] = p_jacobian[2];
    }

  private:

    CELER_FUNCTION void generate_lipschitz()
    {
        const int ux = std::fabs(direction_[0]);
        const int uy = std::fabs(direction_[1]);
        const int uz = std::fabs(direction_[2]);

        // Max of f_d
        lipschitz_first_[0] = a_ * (ux + uy + uz);

        // Max of f_{\theta}
        lipschitz_first_[1] = a_*t_ * sqrt_three;

        // Max of f_{\phi}
        lipschitz_first_[2] = a_*t_ * sqrt_two;

        // Max of f_{\theta d}
        lipschitz_first_[3] = a_*t_ * lipschitz_first_[0];

        // Max of f_{\theta \theta}
        lipschitz_first_[4] = sqrt_three * a_*t_ * (a_*t_ + 1);

        // Max of f_{\theta \phi}
        lipschitz_first_[5] = a_*t_ * (0.5 * a_*t_ + sqrt_two);

        // Max of f_{\phi d}
        lipschitz_first_[6] = a_*a_*t_ * (std::fabs(ux) + std::fabs(uy));

        // Max of f_{\phi \theta}
        lipschitz_first_[7] = lipschitz_first_[5];

        // Max of f_{\phi \phi}
        lipschitz_first_[8] = sqrt_two * a_*t_ * (a_*t_ + 1);

        // Max of f_{dd}
        lipschitz_second_[0] = a_*a_;

        // Max of f_{\theta \theta}
        lipschitz_second_[1] = lipschitz_first_[4];

        // Max of f_{\phi \phi}
        lipschitz_second_[2] = lipschitz_first_[8];

        // Max of f_{\theta dd}
        lipschitz_second_[3] = a_*a_*a_ * t_ * (ux * ux
                + uy * uy
                + uz * uz);

        // Max of f_{\theta \theta \theta}
        lipschitz_second_[4] = a_*t_ * (sqrt_three*a_*a_*t_*t_
                + 1.5*a_*t_*(1 + sqrt_two) + sqrt_three);

        // Max of f_{\theta \phi \phi}
        lipschitz_second_[5] = a_*t_ * (0.5*a_*a_*t_*t_ + 1.5*sqrt_two*t_ + sqrt_two);

        // Max of f_{\phi dd}
        lipschitz_second_[6] = a_*a_*a_ * t_ * (ux * ux
                + uy * uy);

        // Max of f_{\phi \theta \theta}
        lipschitz_second_[7] = a_*t_ * (0.5*a_*a_*t_*t_ + 3*a_*t_ + sqrt_two);

        // Max of f_{\phi \phi \phi}
        lipschitz_second_[8] = a_*t_ * (sqrt_two*a_*a_*t_*t_ + 3*a_*t_ + sqrt_two);
    }

private:
    const Real3 start_, direction_;
    const real_type a_, t_;

    Real9 lipschitz_first_;
    Real9 lipschitz_second_;
};

class GyroidDefinition
{
  public:

    static constexpr int MaxIntersections = 8;

    CELER_FUNCTION GyroidDefinition(
            Real3 start, Real3 direction, real_type t, real_type p)
        : start_{start}
        , direction_{direction}
        , t_{t}
        , a_{2*pi / p}
    {
        generate_lipschitz();
    }

    static inline CELER_FUNCTION void evaluate_point(
            const Real3& start,
            const Real3& direction,
            const real_type t, const real_type a,
            const Real3& point,
            Real3& p_eval,
            Real9& p_jacobian)
    {
        TPMSParams<real_type> params{start, direction, t, a};
        params.set_params(point);

        const real_type d1 = params.cx*params.cy - params.sz*params.sx;
        const real_type d2 = params.cy*params.cz - params.sx*params.sy;
        const real_type d3 = params.cz*params.cx - params.sy*params.sz;
        const Real3 dd1c = {a*(params.sx*params.cy + params.sz*params.cx),
            a*params.cx*params.sy,
            a*params.cz*params.sx};
        const Real3 dd2c = {a*params.cx*params.sy,
            a*(params.sy*params.cz + params.sx*params.cy),
            a*params.cy*params.sz};
        const Real3 dd3c = {a*params.cz*params.sx,
            a*params.cy*params.sz,
            a*(params.sz*params.cx + params.sy*params.cz)};

        auto dd = [&] (int i, const Real3& ddc) { return -params.dX[i] * ddc[0]
                - params.dY[i]*ddc[1]
                - params.dZ[i]*ddc[2];
        };

        p_jacobian[0] = params.dX[0] * d1
            + params.dY[0] * d2
            + params.dZ[0] * d3;

        p_jacobian[1] = params.dX[1] * d1
            + params.dY[1] * d2
            + params.dZ[1] * d3;

        p_jacobian[2] = params.dX[2] * d1
            + params.dY[2] * d2;

        real_type dd1 = dd(0, dd1c);
        real_type dd2 = dd(0, dd2c);
        real_type dd3 = dd(0, dd3c);

        p_jacobian[3] = params.dX[1] * dd1
            + params.dY[1] * dd2
            + params.dZ[1] * dd3;

        p_jacobian[6] = params.dX[2] * dd1
            + params.dY[2] * dd2;

        p_jacobian[4] = params.dX[1] * dd(1, dd1c)
            + params.dY[1] * dd(1, dd2c)
            + params.dZ[1] * dd(1, dd3c)
            + params.dXth[1] * d1
            + params.dYth[1] * d2
            + params.dZth[1] * d3;

        dd1 = dd(2, dd1c);
        dd2 = dd(2, dd2c);
        dd3 = dd(2, dd3c);

        p_jacobian[5] = params.dX[1] * dd1
            + params.dY[1] * dd2
            + params.dZ[1] * dd3
            + params.dXth[2] * d1
            + params.dYth[2] * d2;

        p_jacobian[7] = p_jacobian[5];

        p_jacobian[8] = params.dX[2] * dd1
            + params.dY[2] * dd2
            + params.dXph[2] * d1
            + params.dYph[2] * d2;

        p_eval[0] = params.sx*params.cy + params.sy*params.cz
            + params.sz*params.cx;
        p_eval[1] = p_jacobian[1];
        p_eval[2] = p_jacobian[2];
    }

  private:

    inline CELER_FUNCTION void generate_lipschitz()
    {
        const real_type ux = std::fabs(direction_[0]);
        const real_type uy = std::fabs(direction_[1]);
        const real_type uz = std::fabs(direction_[2]);

        const real_type sqrt_six = std::sqrt(6);

        // Max of f_d
        lipschitz_first_[0] = a_ * sqrt_two * (ux + uy + uz);

        // Max of f_{\theta}
        lipschitz_first_[1] = sqrt_six * a_ * t_;

        // Max of f_{\phi}
        lipschitz_first_[2] = 2 * a_ * t_;

        // Max of f_{\theta d}
        lipschitz_first_[3] = 2 * t_ * (ux + uy + uz);

        // Max of f_{\theta \theta}
        lipschitz_first_[4] = a_*a_*t_*t_ * (1 + sqrt_six + sqrt_two) + a_*t_*sqrt_six;

        // Max of f_{\theta \phi}
        lipschitz_first_[5] = 1.75*sqrt_two*t_*t_ + 2*t_;

        // Max of f_{\phi d}
        lipschitz_first_[6] = t_ * (sqrt_three * (ux + uy) + sqrt_two * uz);

        // Max of f_{\phi \theta}
        lipschitz_first_[7] = lipschitz_first_[5];

        // Max of f_{\phi \phi}
        lipschitz_first_[8] = 3*a_*a_*t_*t_ + 2*a_*t_;

        // Max of f_{dd}
        lipschitz_second_[0] = a_*a_ * (2 * (ux*uy + uy*uz + uz*ux) + sqrt_two);

        // Max of f_{\theta \theta}
        lipschitz_second_[1] = lipschitz_first_[4];

        // Max of f_{\phi \phi}
        lipschitz_second_[2] = lipschitz_first_[8];

        // Max of f_{\theta dd}
        lipschitz_second_[3] = 2 * a_*a_*a_ * t_ * (1 + sqrt_two * (ux*uy + uy*uz + uz*ux));

        // Max of f_{\theta \theta \theta}
        lipschitz_second_[4] = a_*t_ * (a_*a_*t_*t_*(sqrt_six + 1.5*sqrt_two + 3)
                + a_*t_*(1.5*(1 + sqrt_two) + 9) + sqrt_six);

        // Max of f_{\theta \phi \phi}
        lipschitz_second_[5] = a_*t_ * (a_*a_*t_*t_*(0.25*(1 + std::sqrt(5)) + 1.5*sqrt_two)
                + a_*t_*(2*sqrt_two + 2.5) + 2);

        // Max of f_{\phi dd}
        lipschitz_second_[6] = a_*a_*a_ * t_ * (sqrt_three*(ux*ux + uy*uy) + sqrt_two*uz*uz
                + 2*(sqrt_two*ux*uy + uy*uz + uz*ux));

        // Max of f_{\phi \theta \theta}
        lipschitz_second_[7] = a_*t_ * (a_*a_*t_*t_*(2*sqrt_two + 0.25*(1 + std::sqrt(5)))
                + 6*sqrt_two*a_*t_ + 2);

        // Max of f_{\phi \phi \phi}
        lipschitz_second_[8] = a_*t_ * (a_*a_*t_*t_*(2 + 1.5*sqrt_two) + 4.5*sqrt_two*t_ + 2);
    }

    const Real3 start_, direction_;
    const real_type a_, t_;

    Real9 lipschitz_first_;
    Real9 lipschitz_second_;
};

}

