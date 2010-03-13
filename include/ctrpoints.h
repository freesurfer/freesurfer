/**
 * @file  ctrpoints.h
 * @brief read and write control points
 *
 */
/*
 * Original Author: Y.Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/03/13 01:32:40 $
 *    $Revision: 1.5 $
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

#ifndef ctrpoints_h
#define ctrpoints_h

typedef struct
{
  double x;
  double y;
  double z ;
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
int MRIwriteControlPoints(MPoint *pointArray,
                          int count,
                          int useRealRAS,
                          char *fname);

#endif // inclusion guard
