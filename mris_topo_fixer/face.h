/**
 * @file  face.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
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


#ifndef TOPOLOGY_FACE_H
#define TOPOLOGY_FACE_H

#include "globals.h"

class Face
{
public:
  int v[3]; // list of vertices
  int f[3]; // list of neighboring faces

  int marked; //for computational purposes

  double nx,ny,nz; // the face normal
  double x,y,z; // the coordinates of the face

  //constructor/destructor
  Face(void);
  ~Face(void);
  const Face& operator=(const Face &face);
};

#endif
