/**
 * @file  ctrpoints.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
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


//
// ctrpoints.h
//
// purpose: read and write control points
//
#ifndef ctrpoints_h
#define ctrpoints_h

#include "volume_io/basic.h" /* defines Real */

typedef struct
{
  Real x;
  Real y;
  Real z ;
}
MPoint;

// reading control points
// returning array of MGHPoint
// if useRealRAS = 1, then it is in scanner RAS
// if useRealRAS = 0, then it is in surface RAS
MPoint *MRIreadControlPoints(const char *fname, int *count, int *useRealRAS);

// writing control points
// returning array of MGHPoint
// we will write whether they are in scannerRAS or surfaceRAS
// Note that use must tell which coordinate system it is.
int MRIwriteControlPoints(MPoint *pointArray, int count, int useRealRAS, char *fname);

#endif // inclusion guard
