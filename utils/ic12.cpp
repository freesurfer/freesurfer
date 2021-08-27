/**
 * @brief ic0.tri - icosahedron
 *
 * in memory version of ic0.tri
 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include "diag.h"
#include "error.h"
#include "icosahedron.h"
#include "mrisurf.h"

#define ICO_NVERTICES 12
// clang-format off
IC_VERTEX ic12_vertices[ICO_NVERTICES]
    = { { .0000, .0000, 1.0000 },  { .8944, .0000, .4472 },   { .2764, .8507, .4472 },    { -.7236, .5257, .4472 },
        { -.7236, -.5257, .4472 }, { .2764, -.8507, .4472 },  { .7236, -.5257, -.4472 },  { .7236, .5257, -.4472 },
        { -.2764, .8507, -.4472 }, { -.8944, .0000, -.4472 }, { -.2764, -.8507, -.4472 }, { .0000, .0000, -1.0000 } };

#define ICO_NFACES 20

IC_FACE ic12_faces[ICO_NFACES]
    = { { { 1, 5, 4 } },   { { 1, 6, 5 } },    { { 1, 2, 6 } },   { { 1, 3, 2 } },   { { 1, 4, 3 } },
        { { 4, 9, 3 } },   { { 4, 10, 9 } },   { { 4, 5, 10 } },  { { 5, 11, 10 } }, { { 5, 6, 11 } },
        { { 6, 7, 11 } },  { { 6, 2, 7 } },    { { 2, 8, 7 } },   { { 2, 3, 8 } },   { { 3, 9, 8 } },
        { { 9, 10, 12 } }, { { 10, 11, 12 } }, { { 11, 7, 12 } }, { { 7, 8, 12 } },  { { 8, 9, 12 } } };
// clang-format on
MRI_SURFACE *ic12_make_surface(int max_vertices, int max_faces)
{
  static int first_time = 1;

  if (first_time) {
    first_time = 0;
    int vno, fno;
    for (fno = 0; fno < ICO_NFACES; fno++) {
      vno = ic12_faces[fno].vno[1];
      ic12_faces[fno].vno[1] = ic12_faces[fno].vno[2];
      ic12_faces[fno].vno[2] = vno;
    }
    /* position vertices */
    for (vno = 0; vno < ICO_NVERTICES; vno++) {
      ic12_vertices[vno].x *= 100;
      ic12_vertices[vno].y *= 100;
      ic12_vertices[vno].z *= 100;
    }
  }

  ICOSAHEDRON icos =
    {
        ICO_NVERTICES, 
        ICO_NFACES, 
        ic12_vertices,
        ic12_faces
    };
    
  return ICOtoMRIS(&icos, max_vertices, max_faces);
}
