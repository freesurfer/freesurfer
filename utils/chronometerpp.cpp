/**
 * @file  chronometerpp.cpp
 * @brief Holds C++ chronometer stream insertion operator
 *
 * Holds C++ chronometer stream insertion operator
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/04 18:45:48 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
