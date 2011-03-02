/**
 * @file  mrimean_cuda.h
 * @brief Header file for GPU-based MRI mean filter
 *
 * Implements MRI mean
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
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


#ifndef MRI_MEAN_CUDA_H
#define MRI_MEAN_CUDA_H

#include "mri.h"

#if defined(__cplusplus)
extern "C" {
#endif



  MRI* MRImean_cuda( const MRI* src, MRI* dst,
		     const unsigned int wSize );

#if defined(__cplusplus)
};
#endif




#endif
