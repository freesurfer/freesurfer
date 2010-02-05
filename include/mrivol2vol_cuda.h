/**
 * @file  mrivol2vol_cuda.h
 * @brief Header file for GPU-based MRI vol2vol transforms
 *
 * Implements MRI vol2vol transforms.... somewhat messy right now
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
