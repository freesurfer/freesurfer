/**
 * @brief surface2surface test routine
 *
 */
/*
 * Original Author: Y. Tosa
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

#include <iomanip>
#include <iostream>

#include "matrix.h"
#include "mri.h"
#include "mrisurf.h"
#include "transform.h"

const char *Progname = "surf2surf";

#define V4_LOAD(v, x, y, z, r)                                                 \
  (VECTOR_ELT(v, 1) = x, VECTOR_ELT(v, 2) = y, VECTOR_ELT(v, 3) = z,           \
   VECTOR_ELT(v, 4) = r);

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "surf2surf <surf> <dstvol> <xfm> <surfdst>" << std::endl;
    return -1;
  }
  LTA *lta = 0;
  std::cout << "Reading xfm     from " << argv[3] << std::endl;
  if (!stricmp(argv[3], "ID")) // no xfm is given
  {
    std::cout << "Just load identity ras-to-ras" << std::endl;
    lta       = LTAalloc(1, NULL);
    lta->type = LINEAR_RAS_TO_RAS;
  } else
    lta = LTAreadEx(argv[3]);

  std::cout << "Reading surface from " << argv[1] << std::endl;
  MRIS *mris = MRISread(argv[1]);
  if (mris->vg.valid == 0) {
    std::cout << "could not get the original volume info" << std::endl;
    return -1;
  }
  MRIS *mrisDst = MRISclone(mris);
  // build src volume
  MRI *src = MRIallocHeader(mris->vg.width, mris->vg.height, mris->vg.depth,
                            MRI_VOLUME_TYPE_UNKNOWN, 1);
  useVolGeomToMRI(&mris->vg, src);

  std::cout << "Reading dstvol  from " << argv[2] << std::endl;
  MRI *dst = MRIread(argv[2]);
  // modify geometry info
  getVolGeom(dst, &mrisDst->vg);

  VECTOR *sX = VectorAlloc(4, MATRIX_REAL);
  VECTOR *dX = VectorAlloc(4, MATRIX_REAL);

  MATRIX *surf2surf = surfaceRASFromSurfaceRAS_(dst, src, lta);
  MRIfree(&src);
  MRIfree(&dst);

  // now get all the vertex points and change them to the corresponding dst
  // surface vertices
  for (int i = 0; i < mris->nvertices; i++) {
    V4_LOAD(sX, mris->vertices[i].x, mris->vertices[i].y, mris->vertices[i].z,
            1.);
    MatrixMultiply(surf2surf, sX, dX);
    mrisDst->vertices[i].x = VECTOR_ELT(dX, 1);
    mrisDst->vertices[i].y = VECTOR_ELT(dX, 2);
    mrisDst->vertices[i].z = VECTOR_ELT(dX, 3);
  }
  // oh my I have to do the following?
  MRIScomputeNormals(mrisDst);
  mris->radius = MRISaverageRadius(mrisDst);
  MRIScomputeMetricProperties(mrisDst);
  MRISstoreCurrentPositions(mrisDst);
  // write it out
  std::cout << "Writing surface to " << argv[4] << std::endl;
  MRISwrite(mrisDst, argv[4]);

  MatrixFree(&surf2surf);
  MRISfree(&mris);
  MRISfree(&mrisDst);
  VectorFree(&sX);
  VectorFree(&dX);
}
