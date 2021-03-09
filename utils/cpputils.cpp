/**
 * @brief various utilities which make use of C++ classes
 *
 */
/*
 * Original Author: Krish Subramaniam
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

#include "cma.h"
#include "diag.h"
#include "macros.h"
#include "mri.h"
#include "mrisurf.h"

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

#include "MRISOBBTree.h"
#include "MRISdistancefield.h"

typedef Math::Point< int > Pointd;
typedef Math::Point< float > Point_f;
using namespace std;

/* This function finds the closest distance of voxels with the distance d of volume mri_dist to the surface mris.
 */
MRI *MRISsignedFixedDistanceTransform(MRI_SURFACE *mris, MRI *mri_dist, double distance)
{
  // this volume is used to track visited voxels
  MRI* mri_visited = MRIcloneDifferentType(mri_dist, MRI_INT);

  // Save before converting, restored later
  MRISsavedXYZ* savedXYZ = MRISsaveXYZ(mris);

  // Convert surface vertices to vox space
  Math::ConvertSurfaceRASToVoxel(mris, mri_dist);

  // Find the distance field
  MRISDistanceField *distfield = new MRISDistanceField(mris, mri_dist);
  distfield->SetMaxDistance(distance);
  distfield->Generate();  // mri_dist now has the distancefield

  // Construct the OBB Tree
  MRISOBBTree *OBBTree = new MRISOBBTree(mris);
  OBBTree->ConstructTree();

  std::queue< Pointd * > ptsqueue;
  // iterate through all the volume points
  // and apply sign
  for (int i = 0; i < mri_dist->width; i++) {
    // std::cerr << i <<" ";
    for (int j = 0; j < mri_dist->height; j++) {
      for (int k = 0; k < mri_dist->depth; k++) {
        if (MRIIvox(mri_visited, i, j, k)) continue;
        int res = OBBTree->PointInclusionTest(i, j, k);
        Pointd *pt = new Pointd;
        pt->v[0] = i;
        pt->v[1] = j;
        pt->v[2] = k;
        ptsqueue.push(pt);

        // First serve all the points in the queue before going to the next voxel
        while (!ptsqueue.empty()) {
          // serve the front and pop it
          Pointd *p = ptsqueue.front();
          const int x = p->v[0];
          const int y = p->v[1];
          const int z = p->v[2];
          delete p;
          ptsqueue.pop();

          if (MRIIvox(mri_visited, x, y, z)) continue;
          MRIIvox(mri_visited, x, y, z) = res;
          const float dist = MRIFvox(mri_dist, x, y, z);
          MRIFvox(mri_dist, x, y, z) = dist * res;

          // mark its 6 neighbors if distance > 1 ( triangle inequality )
          if (dist > 1) {
            // left neighbor in x
            if (x > 0 && !MRIIvox(mri_visited, x - 1, y, z)) {
              Pointd *ptemp = new Pointd;
              ptemp->v[0] = x - 1;
              ptemp->v[1] = y;
              ptemp->v[2] = z;
              ptsqueue.push(ptemp);
            }
            // bottom neighbor in y
            if (y > 0 && !MRIIvox(mri_visited, x, y - 1, z)) {
              Pointd *ptemp = new Pointd;
              ptemp->v[0] = x;
              ptemp->v[1] = y - 1;
              ptemp->v[2] = z;
              ptsqueue.push(ptemp);
            }
            // front neighbor in z
            if (z > 0 && !MRIIvox(mri_visited, x, y, z - 1)) {
              Pointd *ptemp = new Pointd;
              ptemp->v[0] = x;
              ptemp->v[1] = y;
              ptemp->v[2] = z - 1;
              ptsqueue.push(ptemp);
            }
            // right neighbor in x
            if (x < mri_visited->width - 1 && !MRIIvox(mri_visited, x + 1, y, z)) {
              Pointd *ptemp = new Pointd;
              ptemp->v[0] = x + 1;
              ptemp->v[1] = y;
              ptemp->v[2] = z;
              ptsqueue.push(ptemp);
            }
            // top neighbor in y
            if (y < mri_visited->height - 1 && !MRIIvox(mri_visited, x, y + 1, z)) {
              Pointd *ptemp = new Pointd;
              ptemp->v[0] = x;
              ptemp->v[1] = y + 1;
              ptemp->v[2] = z;
              ptsqueue.push(ptemp);
            }
            // back neighbor in z
            if (z < mri_visited->depth - 1 && !MRIIvox(mri_visited, x, y, z + 1)) {
              Pointd *ptemp = new Pointd;
              ptemp->v[0] = x;
              ptemp->v[1] = y;
              ptemp->v[2] = z + 1;
              ptsqueue.push(ptemp);
            }
          }
        }
      }
    }
  }
  
  MRISpopXYZ(mris,&savedXYZ);

  delete OBBTree;
  delete distfield;
  MRIfree(&mri_visited);
  return (mri_dist);
}


/* This uses the more accurate signed distance transform from a surface ( uses OBB Trees in turn ). It takes more time (
 * around 4 to 5 minutes ) but it's very accurate. returned MRI structure (mri_out) is an MRI_INT MRI volume with voxels
 * inside having the value 1 and voxels outside with value 0.
 */
MRI *MRISfillInterior2(MRI_SURFACE *mris, MRI *mri_interior)
{
  MRI *_mridist, *mri_out;

  _mridist = MRIcloneDifferentType(mri_interior, MRI_FLOAT);
  mri_out = MRIclone(mri_interior, NULL);

  MRISsignedFixedDistanceTransform(mris, _mridist, 3.0);

  for (int i = 0; i < _mridist->width; i++) {
    for (int j = 0; j < _mridist->height; j++) {
      for (int k = 0; k < _mridist->depth; k++) {
        if (MRIFvox(_mridist, i, j, k) > 0.0)
          MRIsetVoxVal(mri_out, i, j, k, 0, 1);
        else
          MRIsetVoxVal(mri_out, i, j, k, 0, 0);
      }
    }
  }

  MRIfree(&_mridist);
  return (mri_out);
}
