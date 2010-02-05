/**
 * @file  mriconvolve1d_cuda.h
 * @brief Header file for GPU-based MRI 1D convolutions
 *
 * Implements 1D MRI convolutions.... somewhat messy right now
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/05 16:14:19 $
 *    $Revision: 1.4 $
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

#ifndef MRI_CONVOLVE_1D_CUDA_H
#define MRI_CONVOLVE_1D_CUDA_H




#include "mri.h"


#if defined(__cplusplus)
extern "C" {
#endif



  MRI* MRIconvolve1d_cuda( const MRI* src, MRI* dst,
			   const float *kernel,
			   const unsigned int kernelLength,
			   const int axis,
			   const int srcFrame, const int dstFrame );
    
  MRI* MRIconvolveGaussian_cuda( const MRI* src, MRI* dst,
				 const float *kernel,
				 const unsigned int kernelLength );


#if defined(__cplusplus)
};
#endif



#endif
