/**
 * @brief read and write control points
 *
 */
/*
 * Original Author: Y.Tosa
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef ctrpoints_h
#define ctrpoints_h

#include "transform.h"

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
int MRIwriteControlPoints(const MPoint *pointArray,
                          int count,
                          int useRealRAS,
                          const char *fname);
													
// map control points using LTA (contains src and trg geometries):
MPoint *MRImapControlPoints(const MPoint *pointArray, int count, int useRealRAS,
                            MPoint *trgArray, LTA* lta);
MPoint *ControlPoints2Vox(MPoint *ras, int npoints, int UseRealRAS, MRI *vol);	

MATRIX *ControlPoints2TalMatrix(char *subject);
MPoint *ControlPointsApplyMatrix(MPoint *srcctr, int nctrpoints, MATRIX *M, MPoint *outctr);
MPoint *GetTalControlPoints(char **subjectlist, int nsubjects, int *pnctrtot);
MPoint *GetTalControlPointsSFile(const char *subjectlistfile, int *pnctrtot);

#endif // inclusion guard
