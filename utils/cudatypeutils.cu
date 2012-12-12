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
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iomanip>

#include "cudatypeutils.hpp"


// =======================================================
// Stream operators


std::ostream& operator<<( std::ostream& os, const float3& me )
{
  const unsigned int w = 6;
  const unsigned int p = 2;

  os << "(";
  os << std::setw(w) << std::setprecision(p) << me.x;
  os << std::setw(w) << std::setprecision(p) << me.y;
  os << std::setw(w) << std::setprecision(p) << me.z;
  os << " )";
  return( os );
}

std::ostream& operator<<( std::ostream& os, const dim3& me )
{
  const unsigned int w = 6;
  os << "(";
  os << std::setw(w) << me.x;
  os << std::setw(w) << me.y;
  os << std::setw(w) << me.z;
  os << " )";
  return( os );
}
