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

#include <stdio.h>
#include "diag.h"
#include "error.h"
#include "icosahedron.h"
#include "mrisurf.h"
#include "utils.h"  // fgetl

#define ICO_NVERTICES 163842
#define ICO_NFACES 327680

static int static_read_icosahedron(void);
IC_VERTEX *ic163842_vertices = NULL;
IC_FACE *ic163842_faces = NULL;

MRI_SURFACE *ic163842_make_surface(int max_vertices, int max_faces)
{
  static_read_icosahedron();

  int fno;
  for (fno = 0; fno < ICO_NFACES; fno++) {
    int vno = ic163842_faces[fno].vno[1];
    ic163842_faces[fno].vno[1] = ic163842_faces[fno].vno[2];
    ic163842_faces[fno].vno[2] = vno;
  }

  /* position vertices */
  int vno;
  for (vno = 0; vno < ICO_NVERTICES; vno++) {
    ic163842_vertices[vno].x *= 100;
    ic163842_vertices[vno].y *= 100;
    ic163842_vertices[vno].z *= 100;
  }

  ICOSAHEDRON icos =
    {
        ICO_NVERTICES, 
        ICO_NFACES, 
        ic163842_vertices,
        ic163842_faces
    };
    
  MRIS* mris = ICOtoMRIS(&icos, max_vertices, max_faces);

  freeAndNULL(ic163842_vertices);
  freeAndNULL(ic163842_faces);

  return (mris);
}

static int static_read_icosahedron(void)
{
  FILE *fp;
  char line[200], *cp;
  int vno, fno, vno1, vno2, vno3, n;
  float x, y, z;
  IC_VERTEX *ic_vertices;
  IC_FACE *ic_faces;

  fp = fopen("ic163842.tri", "r");
  if (!fp) ErrorExit(ERROR_NOFILE, "static_read_icosahedron: could not open %s", "ic163842.tri");

  ic_vertices = ic163842_vertices = (IC_VERTEX *)calloc(ICO_NVERTICES, sizeof(IC_VERTEX));
  if (!ic163842_vertices) ErrorExit(ERROR_NOMEMORY, "static_read_ico: could not allocate vertex list");
  ic_faces = ic163842_faces = (IC_FACE *)calloc(ICO_NFACES, sizeof(IC_FACE));
  if (!ic163842_faces) ErrorExit(ERROR_NOMEMORY, "static_read_ico: could not allocate vertex list");

  fgetl(line, 150, fp); /* discard # of vertices */

  /* first read vertices */
  n = 0;
  while ((cp = fgetl(line, 150, fp)) != NULL) {
    if (sscanf(cp, "%d %f %f %f\n", &vno, &x, &y, &z) < 4) break;
    ic_vertices[vno - 1].x = x;
    ic_vertices[vno - 1].y = y;
    ic_vertices[vno - 1].z = z;
    if (++n >= ICO_NVERTICES) break;
  }
  n = 0;
  fgetl(line, 150, fp); /* discard # of faces */
  while ((cp = fgetl(line, 150, fp)) != NULL) {
    if (sscanf(cp, "%d %d %d %d\n", &fno, &vno1, &vno2, &vno3) < 4) break;
    ic_faces[fno - 1].vno[0] = vno1;
    ic_faces[fno - 1].vno[1] = vno2;
    ic_faces[fno - 1].vno[2] = vno3;
    if (++n >= ICO_NFACES) break;
  }
  fclose(fp);
  return (NO_ERROR);
}
