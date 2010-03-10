/**
 * @file  gcamorphtestutils.h
 * @brief Utilities to help with testing GCAmorph routines
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/10 20:49:08 $
 *    $Revision: 1.2 $
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

#ifndef GCA_MORPH_TEST_UTILS_H
#define GCA_MORPH_TEST_UTILS_H

#include "gcamorph.h"

#if defined(__cplusplus)
extern "C" {
#endif


  //! Writes the parts required by gcamComputeMetricProperties
  void WriteGCAMforMetricProperties( const GCAM* src, const char* fName );
  //! Reads the parts required by gcamComputeMetricProperties
  void ReadGCAMforMetricProperties( GCAM** dst, const char* fName );




#if defined(__cplusplus)
};
#endif




#endif
