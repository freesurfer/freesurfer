/**
 * @file  gcamfots_cuda.h
 * @brief C interface to gcamFindOptimalTimeStep for the GPU
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
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

#ifndef GCA_MORPH_FOTS_H
#define GCA_MORPH_FOTS_H

#include "gcamorph.h"
#include "mri.h"

#if defined(__cplusplus)
extern "C" {
#endif

  //! Wrapper, to be called form gcamFindOptimalTimestep
  float gcamFindOptimalTimestepGPU( GCA_MORPH *gcam,
				    GCA_MORPH_PARMS *parms,
				    MRI *mri );

#if defined(__cplusplus)
};
#endif


#endif
