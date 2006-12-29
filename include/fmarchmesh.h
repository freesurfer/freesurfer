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
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.3 $
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

