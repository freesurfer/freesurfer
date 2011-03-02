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

#ifndef TEST_GPU_H
#define TEST_GPU_H

#include "gcamorph.h"


#if defined(__cplusplus)
extern "C" {
#endif

  void GCAMorphSendBefore( const GCAM* src );
  void GCAMorphSendAfter( const GCAM* src );
  void GCAMorphCompareBeforeAfter( GCAM* dst );
#if defined(__cplusplus)
};
#endif


#endif
