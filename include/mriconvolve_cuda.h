/**
 * @file  mriconvolve1d_cuda.h
 * @brief Header file for GPU-based MRI 1D convolutions
 *
 * Implements 1D MRI convolutions.... somewhat messy right now
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.5 $
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
