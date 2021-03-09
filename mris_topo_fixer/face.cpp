/*
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


#include "face.h"

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
