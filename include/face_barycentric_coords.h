#pragma once
/**
 * @file  face_barycentric_coords.h
 * @brief MRI_SURFACE utilities.
 *
 * Utilities, constants and structure definitions for manipulation
 * and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl, extracted by Bevin Brett
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2017/02/16 19:42:54 $
 *    $Revision: 1.391 $
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
  auto v0   = face.v(0);
  auto v1   = face.v(1);
  auto v2   = face.v(2);
  
  double V0[3], V1[3], V2[3];

  // Done this way because might be faster than v0.which_coords(...) ...
  // but should be pulled out of this specific function and used more widely
  
  switch (which_vertices) {
    default:
      ErrorExit(ERROR_BADPARM, "face_barycentric_coords: which %d not supported", which_vertices);
      break;
    case FLATTENED_VERTICES:
      V0[0] = v0.fx();
      V0[1] = v0.fy();
      V0[2] = v0.fz();
      V1[0] = v1.fx();
      V1[1] = v1.fy();
      V1[2] = v1.fz();
      V2[0] = v2.fx();
      V2[1] = v2.fy();
      V2[2] = v2.fz();
      break;
    case CURRENT_VERTICES:
      V0[0] = v0.x();
      V0[1] = v0.y();
      V0[2] = v0.z();
      V1[0] = v1.x();
      V1[1] = v1.y();
      V1[2] = v1.z();
      V2[0] = v2.x();
      V2[1] = v2.y();
      V2[2] = v2.z();
      break;
    case PIAL_VERTICES:
      V0[0] = v0.pialx();
      V0[1] = v0.pialy();
      V0[2] = v0.pialz();
      V1[0] = v1.pialx();
      V1[1] = v1.pialy();
      V1[2] = v1.pialz();
      V2[0] = v2.pialx();
      V2[1] = v2.pialy();
      V2[2] = v2.pialz();
      break;
    case CANONICAL_VERTICES:
      V0[0] = v0.cx();
      V0[1] = v0.cy();
      V0[2] = v0.cz();
      V1[0] = v1.cx();
      V1[1] = v1.cy();
      V1[2] = v1.cz();
      V2[0] = v2.cx();
      V2[1] = v2.cy();
      V2[2] = v2.cz();
      break;
  }

  return face_barycentric_coords(V0, V1, V2, cx, cy, cz, pl1, pl2, pl3);
}
