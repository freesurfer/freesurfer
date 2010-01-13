/**
 * @file  mrivol2vol_cuda.cu
 * @brief Holds MRI volume to volume transform for CUDA
 *
 * Implements MRIvol2vol for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/01/13 20:59:52 $
 *    $Revision: 1.2 $
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
using namespace std;


extern "C" {
#include "mri.h"
}


#include "chronometer.hpp"
#include "cudacheck.h"

#include "mri_transfer.h"
#include "mrislices_cuda.h"


#include "mrivol2vol_cuda.h"

// ==============================================================


int MRIvol2vol_cuda( const MRI* src, MRI* targ, 
		     const float* transformMatrx,
		     const int InterpMode,
		     const float param ) {
  /*!
    Reimplementation of MRIvol2vol for the GPU.
    Is meant to be called from within MRIvol2vol,
    once that routine has performed necessary checks
  */

  cout << __FUNCTION__ << ": Begin" << endl;

  // Need to track 'actual' and GPU dimensions
  // Recall that GPU dimensions will be padded
  // to multiples of 16
  uint4 srcDimsAct, srcDimsGPU;
  uint4 targDimsAct, targDimsGPU;

  srcDimsAct = MRI_GetSlicesDims( src );
  srcDimsGPU = MRI_GetSlicesDimsonGPU( src );
  targDimsAct = MRI_GetSlicesDims( targ );
  targDimsGPU = MRI_GetSlicesDimsonGPU( targ );

  // Get the source MRI onto the GPU

  // Promote source MRI data to float on the GPU

  // Release original src data on GPU

  // Allocate destination float array on GPU

  // Allocate CUDA array to hold one frame of data


  // Loop over frames
  for( unsigned int iFrame=0; iFrame < src->nframes; iFrame++ ) {
    // Copy src frame into CUDA array

    // Bind the array to a texture

    // Run the kernel on this frame

    // Unbind texture and array
  }

  // Release CUDA array

  // Release src float data

  // Convert destination data back to desired data type

  // Get destination data back from the GPU

  return( 0 );
}
    
