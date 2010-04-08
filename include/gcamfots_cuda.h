/**
 * @file  gcamfots_cuda.h
 * @brief C interface to gcamFindOptimalTimeStep for the GPU
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/04/08 14:09:32 $
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
