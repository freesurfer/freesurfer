/**
 * @file  mrivol2vol_cuda.h
 * @brief Header file for GPU-based MRI vol2vol transforms
 *
 * Implements MRI vol2vol transforms.... somewhat messy right now
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

#ifndef MRI_VOL_2_VOL_CUDA_H
#define MRI_VOL_2_VOL_CUDA_H


#include "mri.h"
#include "matrix.h"

#if defined(__cplusplus)
extern "C" {
#endif


  //! Implementation of MRIvol2vol for the GPU
  int MRIvol2vol_cuda( const MRI* src, MRI* targ, 
		       const MATRIX* transformMatrix,
		       const int InterpMode,
		       const float param );

  //! Print timers
  void MRIvol2volShowTimers( void );
    

#if defined(__cplusplus)
};
#endif



#endif
