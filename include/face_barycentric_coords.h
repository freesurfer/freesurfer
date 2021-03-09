#pragma once
/**
 * @brief MRI_SURFACE utilities.
 *
 * Utilities, constants and structure definitions for manipulation
 * and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl, extracted by Bevin Brett
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



#include "mrisurf.h"

template <class Surface>
int face_barycentric_coords_template(
  Surface surface,
  int     fno,
  int     which_vertices,
  double  cx,
  double  cy,
  double  cz,
  double *pl1,
  double *pl2,
  double *pl3)
{
  auto face = surface.faces(fno);
  
  double d[3][3];

  for (int i = 0; i < 3; i++) {
    float x,y,z;
    face.v(i).which_coords(which_vertices, &x,&y,&z);
    d[i][0] = x;
    d[i][1] = y;
    d[i][2] = z;
  }
  
  return face_barycentric_coords(d[0], d[1], d[2], cx, cy, cz, pl1, pl2, pl3);
}
