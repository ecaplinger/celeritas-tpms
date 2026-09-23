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

namespace celeritas {

using Real9 = Array<real_type, 9>;

/*
 *!
 * Find a root using Newton's method for a system of 3 equations and 3 unknowns.
 */
class NewtonSolver3
{
public:
    struct NewtonCalcs {
        Real3 eval;
        Real9 jacobian;
    };

    CELER_FUNCTION NewtonSolver3(
            void(* func)(const Real3&, Real3&, Real9&))
        : func_{func}, tol_{default_tol_}
    {
    }

    CELER_FUNCTION NewtonSolver3(
            void(* func)(const Real3&, Real3&, Real9&), real_type tol)
        : func_{func}, tol_{tol}
    {
    }

    /*!
     * Solve for a root given an initial guess.
     * Perform Newton iterations \f[
     * x_{k+1} = x_{k} - J_{f} (x_{k}) f(x_{k})
     * \f] until all elements of the Newton step \f[
     * \left| J_{f} (x) f(x) \right| < \epsilon
     * \f].
     *
     * \return True if we are able to find a root, False if we burn through
     * all allowed iterations without finding anything, or we encounter a
     * singular Jacobian evaluation.
     */
    inline CELER_FUNCTION bool operator()(const Real3& x0, Real3& r)
    {
        Array<real_type, 3> x{x0};

        for (int i = 0; i < max_iters_; i++)
        {
            Real3 eval;
            Real9 jacobian;
            func_(x, eval, jacobian);
            Real3 step;

            bool result = cramer_3x3(jacobian, eval, step);
            if (!result)
                return false;

            for (int j = 0; j < 3; j++)
                x[j] -= step[j];

            if (std::fabs(step[0]) < tol_
                    && std::fabs(step[1]) < tol_
                    && std::fabs(step[2]) < tol_)
            {
                r = x;
                return true;
            }
        }

        return false;
    }

private:

    //! Default tolerance for Newton solve, taken from Orange `Tolerance`.
    static constexpr real_type default_tol_
        = (std::is_same_v<real_type, double> ? 1e-5 : 5e-2f);

    //! Maximum amount of iterations.
    static constexpr int max_iters_ = 5;

    real_type tol_;
    void(* func_)(const Real3&, Real3&, Real9&);

    /*
     *!
     * Perform a matrix solve
     * \f[
     * Ax = b
     * \f]
     * where A is a 3x3 matrix and b is a 1x3 column vector.
     * We hardcode this because it is expected to be faster than whatever
     * general matrix solve would be implemented in a 3rd party library.
     *
     * \return True if operation is successful, False if A is singular.
     */
    static inline CELER_FUNCTION bool cramer_3x3(
            const Real9& a, const Real3& b, Real3& r)
    {
        // a(ei - fh) - b(di - fg) + c(dh - eg)
        auto det_cols = [] (
                const real_type* c1,
                const real_type* c2,
                const real_type* c3) {
            return c1[0] * (c2[1]*c3[2] - c3[1]*c2[2])
                + c2[0] * (c1[1]*c3[2] - c3[1]*c1[2])
                + c3[0] * (c1[1]*c2[2] - c2[1]*c1[2]);
        };

        const real_type* ac1 = &a.data()[0];
        const real_type* ac2 = &a.data()[3];
        const real_type* ac3 = &a.data()[6];
        const real_type* bp = b.data();

        real_type det = det_cols(ac1, ac2, ac3);
        if (det == 0)
            return false;

        r[0] = det_cols(bp, ac2, ac3) / det;
        r[1] = det_cols(ac1, bp, ac3) / det;
        r[2] = det_cols(ac1, ac2, bp) / det;

        return true;
    }
};

}
