/**
 * @file  testgpu.h
 * @brief A place to chuck simple GPU test functions
 *
 * A place to chuck simple GPU test functions, which have to be
 * called from the main freesurfer code. Likely to be buggy
 * and error-prone, so should never appear in production side.
 * If this file is non-empty, it's because something is being tested
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/02 20:36:42 $
 *    $Revision: 1.3 $
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
 *
 */

#ifndef TEST_GPU_H
#define TEST_GPU_H

#include "gcamorph.h"


#if defined(__cplusplus)
extern "C" {
#endif

  void GCAMorphSendBefore( const GCAM* src );
  void GCAMorphSendAfter( const GCAM* src );
  void GCAMorphCompareBeforeAfter( void );
#if defined(__cplusplus)
};
#endif


#endif
