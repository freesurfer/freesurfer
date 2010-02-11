/**
 * @file  cudatypeutils.cu
 * @brief Holds utility routines for the native CUDA types
 *
 * Holds things like mathematical operators and stream operators
 * for the native CUDA types
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/11 15:29:56 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iomanip>

#include "cudatypeutils.hpp"


// =======================================================
// Stream operators


std::ostream& operator<<( std::ostream& os, const float3& me ) {
  const unsigned int w = 6;
  const unsigned int p = 2;

  os << "(";
  os << std::setw(w) << std::setprecision(p) << me.x;
  os << std::setw(w) << std::setprecision(p) << me.y;
  os << std::setw(w) << std::setprecision(p) << me.z;
  os << " )";
  return( os );
}

std::ostream& operator<<( std::ostream& os, const dim3& me ) {
  const unsigned int w = 6;
  os << "(";
  os << std::setw(w) << me.x;
  os << std::setw(w) << me.y;
  os << std::setw(w) << me.z;
  os << " )";
  return( os );
}
