//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/surf/Toroid.hh
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>

#include "corecel/cont/Array.hh"
#include "corecel/cont/Span.hh"
#include "corecel/math/Algorithms.hh"
#include "corecel/math/ArrayOperators.hh"
#include "corecel/math/ArrayUtils.hh"
#include "corecel/math/FerrariSolver.hh"
#include "orange/OrangeTypes.hh"
#include "orange/SenseUtils.hh"

#include "orange/surf/detail/TPMSDefinition.hh"
#include "orange/surf/detail/TPMSIntersector.hh"

namespace celeritas
{

/*!
 * Axis-aligned triply periodic minimal surface.
 *
 * This is a template class because multiple kinds of TPMS are supported;
 * currently Schwarz TPMS and gyroids are implemented. See TPMSDefinition.hh
 * for examples of implementations.
 *
 * The currently implemented TPMS types all have cubic unit cells, so this only
 * takes a single dimensional input for unit_size; for 
 */
template<class SurfDefinition>
class TPMSLattice
{
  public:

    explicit CELER_FUNCTION TPMSLattice(
            Real3 lower_corner, real_type unit_size)
        : cell_intersector_(
    {
    }


  private:
    SurfDefinition tpms_def_;
    TPMSUnitCellIntersector<SurfDefinition> cell_intersector_;

    Real3 model_to_

};

}
