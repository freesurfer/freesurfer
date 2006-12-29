/**
 * @file  face.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:19 $
 *    $Revision: 1.2 $
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


#include "topology/face.h"

Face::Face(void) {
  marked=0;
}

Face::~Face(void) {}

const Face &Face::operator=(const Face &face) {
  for (int n = 0 ; n < 3 ; n++) {
    v[n] = face.v[n];
    f[n] = face.f[n];
  }
  marked=face.marked;
  nx=face.nx;
  ny=face.ny;
  nz=face.nz;
  x=face.x;
  y=face.y;
  z=face.z;

  return face;
}
