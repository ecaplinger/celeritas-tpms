//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file corecel/math/Turn.hh
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>
#include <cstdlib>
#include <type_traits>
#include <algorithm>

#include "corecel/Constants.hh"
#include "corecel/Types.hh"
#include "corecel/math/Algorithms.hh"
#include "corecel/math/NumericLimits.hh"

namespace celeritas {
using constants::pi;
auto nan = numeric_limits<real_type>::quiet_NaN();

//---------------------------------------------------------------------------//
/*!
 * A data type that stores an interval and provides interval-arithmetic
 * operators.
 *
 * Arithmetic operators take sets as operands instead of scalar values, and
 * provide the set of possible results as their output.
 */
class Interval {
public:
    CELER_FUNCTION Interval(real_type a, real_type b)
    : a_(a), b_(b)
    {
    }

    CELER_FUNCTION Interval(real_type a)
        : a_(a), b_(a)
    {
    }

    CELER_FUNCTION Interval()
        : a_(nan), b_(nan)
    {
    }

    //// ACCESSORS ////

    //! Lesser of the two interval values
    CELER_FUNCTION real_type a() const { return a_; }

    //! Greater of the two interval values
    CELER_FUNCTION real_type b() const { return b_; }

    CELER_FUNCTION bool is_null() const { return a_ == nan; }

    //! Distance between the two values
    CELER_FORCEINLINE_FUNCTION real_type width() const { return b_ - a_; }

    CELER_FORCEINLINE_FUNCTION real_type midpoint() const { return a_ + width() * 0.5; }

    //// OPERATOR OVERLOADS ////

    //! In-place addition with another interval
    CELER_FUNCTION Interval& operator+=(const Interval& rhs)
    {
        a_ += rhs.a_;
        b_ += rhs.b_;
        return *this;
    }

    //! In-place addition with a scalar value
    CELER_FUNCTION Interval& operator+=(const real_type rhs)
    {
        a_ += rhs;
        b_ += rhs;
        return *this;
    }

    //! Unary negative operator
    CELER_FUNCTION friend Interval operator-(Interval rhs) {
        return Interval(-rhs.b_, -rhs.a_);
    }

    //! In-place subtraction with another interval
    CELER_FUNCTION Interval& operator-=(const Interval& rhs)
    {
        a_ -= rhs.b_;
        b_ -= rhs.a_;
        return *this;
    }

    //! In-place subtraction with a scalar value
    CELER_FUNCTION Interval& operator-=(const real_type rhs)
    {
        a_ -= rhs;
        b_ -= rhs;
        return *this;
    }

    //! In-place multiplication with another interval
    CELER_FUNCTION Interval& operator*=(const Interval& rhs)
    {
        const auto vals = {
            a_ * rhs.a_,
            a_ * rhs.b_,
            b_ * rhs.a_,
            b_ * rhs.b_
        };

        a_ = *std::min_element(vals.begin(), vals.end());
        b_ = *std::max_element(vals.begin(), vals.end());
        return *this;
    }

    //! In-place multiplication with a scalar value
    CELER_FUNCTION Interval& operator*=(const real_type rhs)
    {
        real_type prod_a = a_ * rhs;
        real_type prod_b = b_ * rhs;

        a_ = std::min(prod_a, prod_b);
        b_ = std::max(prod_a, prod_b);
        return *this;
    }

    //! Binary addition with interval
    CELER_FUNCTION friend Interval operator+(Interval lhs, const Interval& rhs)
    {
        lhs += rhs;
        return lhs;
    }

    //! Binary addition with scalar
    CELER_FUNCTION friend Interval operator+(Interval lhs, const real_type rhs)
    {
        lhs += rhs;
        return lhs;
    }

    //! Binary addition with scalar
    CELER_FUNCTION friend Interval operator+(const real_type lhs, Interval rhs)
    {
        rhs += lhs;
        return rhs;
    }

    //! Binary subtraction with interval
    CELER_FUNCTION friend Interval operator-(Interval lhs, const Interval& rhs)
    {
        lhs -= rhs;
        return lhs;
    }

    //! Binary subtraction with scalar
    CELER_FUNCTION friend Interval operator-(Interval lhs, const real_type rhs)
    {
        lhs -= rhs;
        return lhs;
    }

    //! Binary subtraction with scalar
    CELER_FUNCTION friend Interval operator-(const real_type lhs, Interval rhs)
    {
        rhs = -rhs + lhs;
        return rhs;
    }

    //! Binary multiplication with interval
    CELER_FUNCTION friend Interval operator*(Interval lhs, const Interval& rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    //! Binary multiplication with scalar
    CELER_FUNCTION friend Interval operator*(Interval lhs, const real_type rhs)
    {
        lhs *= rhs;
        return lhs;
    }

    //! Binary multiplication with scalar
    CELER_FUNCTION friend Interval operator*(const real_type lhs, Interval rhs)
    {
        rhs *= lhs;
        return rhs;
    }

    //! Unary intersection with another interval
    CELER_FUNCTION Interval& intersect(const Interval& other)
    {
        if (other.b_ < a_ || b_ < other.a_) {
            a_ = nan;
            b_ = nan;
            return *this;
        }

        a_ = std::max(a_, other.a_);
        b_ = std::min(b_, other.b_);
        return *this;
    }

    //! Intersection of two intervals
    CELER_FUNCTION static inline Interval intersection(Interval lhs, const Interval& rhs)
    {
        lhs.intersect(rhs);
        return lhs;
    }

    //! Interval sine
    CELER_FUNCTION static inline Interval isin(const Interval& lhs);

    //! Interval cosine
    CELER_FUNCTION static inline Interval icos(const Interval& lhs);

    CELER_FUNCTION static inline Interval null_interval()
    {
        Interval r = Interval();
        r.a_ = nan;
        r.b_ = nan;
        return r;
    }


private:

    real_type a_;
    real_type b_;
};

CELER_FUNCTION Interval Interval::icos(const Interval& rhs) {
    Constant twopi = 2 * pi;

    if (rhs.width() > twopi) {
        return Interval(-1, 1);
    }

    real_type k = std::floor(rhs.a_ / twopi);
    real_type a = rhs.a_ - k * twopi;
    real_type b = rhs.b_ - k * twopi;

    real_type cos_a = std::cos(a);
    real_type cos_b = std::cos(b);

    real_type lower = std::min(cos_a, cos_b);
    real_type upper = std::max(cos_a, cos_b);

    if (a <= pi && b >= pi) {
        lower = -1;
    }

    if (b >= twopi) {
        upper = 1;
    }

    return Interval(lower, upper);
}

CELER_FUNCTION Interval Interval::isin(const Interval& rhs)
{
    Constant twopi = 2 * pi;

    if (rhs.width() > twopi) {
        return Interval(-1, 1);
    }

    real_type k = std::floor(rhs.a_ / twopi);
    real_type a = rhs.a_ - k * twopi;
    real_type b = rhs.b_ - k * twopi;

    real_type sin_a = std::sin(a);
    real_type sin_b = std::sin(b);

    real_type lower = std::min(sin_a, sin_b);
    real_type upper = std::max(sin_a, sin_b);

    Constant limit_pos = pi / 2;
    Constant limit_neg = 3 * pi / 2;

    if (a <= limit_pos && limit_pos <= b) {
        upper = 1;
    }

    if (a <= limit_neg && limit_neg <= b) {
        lower = -1;
    }

    return Interval(lower, upper);
}

}  // namespace celeritas
