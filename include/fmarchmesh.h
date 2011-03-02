/**
 * @file  fmarchmesh.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.4 $
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


/* Header function for computing geodesic distance on a triangular
 * surfmesh
 */

#ifndef FMARCHING_H
#define FMARCHING_H

#include <math.h>
#include "heap.h"
#include "mrisurf.h"

#define ALIVE 1
#define NBAND 2
#define FAWAY 3

#undef INFINITY
#define INFINITY 1000000000

typedef struct
{
  float x;
  float y;
  float z;
}
MyVector;

/* function declarations */
float *FastMarchMesh(MRI_SURFACE *mesh, int *contour, int numinitvert, float thred);


#endif

