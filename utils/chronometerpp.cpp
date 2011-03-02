/**
 * @file  chronometerpp.cpp
 * @brief Holds C++ chronometer stream insertion operator
 *
 * Holds C++ chronometer stream insertion operator
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
 *    $Revision: 1.2 $
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

#include "chronometer.hpp"

namespace SciGPU {
  namespace Utilities {

std::ostream& operator<<( std::ostream& os,
			  const SciGPU::Utilities::Chronometer& timer ) {
  
  os << std::setw(9)
     << std::setprecision(6)
     << timer.GetAverageTime() << " ms (avg) ";

  os << std::setw(9)
     << std::setprecision(6)
     << timer.GetTime() << " ms (tot)";

  return( os );
}

  }
}
