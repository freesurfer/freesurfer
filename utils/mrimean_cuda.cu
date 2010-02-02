/**
 * @file  mrimean_cuda.cu
 * @brief Holds MRI mean function for the GPU
 *
 * Implements MRI mean function on the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/02 17:34:19 $
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

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>


extern "C" {
#include "mri.h"
}

//#include "cuPrintf.cu"

#include "chronometer.hpp"
#include "cudacheck.h"


#include "mriframegpu.hpp"


#include "mrimean_cuda.h"






namespace GPU {
  namespace Algorithms {

    const unsigned int kMRImeanBlockSize = 16;

  }
}


MRI* MRImean_cuda( const MRI* src, MRI* dst,
		   const unsigned int wSize ) {
  /*!
    Wrapper around the GPU routine, to be called from the
    original MRImean routine
  */


  return( dst );
}
