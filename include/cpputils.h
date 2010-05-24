/**
 * @file  cpputils.h
 * @brief include files for cpp utils
 *
 */
/*
 * Original Author: Krish Subramaniam
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2010/05/24 15:36:53 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

// include guard
#ifndef cpputils_h
#define cpputils_h

# ifdef __cplusplus
extern "C"
{
#endif


// The following is usable from C
MRI* MRISfillInterior2(MRI_SURFACE *mris, MRI* mri_interior);

#ifdef __cplusplus
}
#endif

#endif 
