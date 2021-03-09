/**
 * @brief Computes accurate volume fractions remaining within a surface
 *
 * This program computes an accurate estimate of the fraction of the volume
 * remaining wihtin a surface.
 */
/*
 * Original Author: Ender Konukoglu
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
#include <stdlib.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "cmdargs.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "macros.h"
#include "mris_compVolFrac.h"
#include "mrisurf.h"
#include "mrisurf_metricProperties.h"
#include "utils.h"
#include "version.h"

MRI *MRIcomputeVolumeFractionFromSurface(MRI_SURFACE *mris, double acc, MRI *mri_src, MRI *mri_fractions)
{
  const int width = mri_src->width;
  const int height = mri_src->height;
  const int depth = mri_src->depth;
  int x, y, z, vno;
  double xs, ys, zs, dist;
  MRIS_HASH_TABLE *mht;

  /* preparing the output */
  printf("preparing the output\n");
  if (mri_fractions == NULL) {
    mri_fractions = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_fractions);
  }
  MRI *mri_shell, *mri_interior;
  /* creating a shell from the surface */
  printf("computing the shell\n");
  mri_shell = MRIclone(mri_src, NULL);
  mri_shell = MRISshell(mri_src, mris, mri_shell, 1);
  /* creating an interior image from the surface */
  printf("computing an interior image\n");
  mri_interior = MRIclone(mri_src, NULL);
  MRIclear(mri_interior);
  mri_interior = MRISfillInterior(mris, 0.0, mri_interior);
  /* creating the hash table related to the surface vertices */
  printf("computing the hash table\n");
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 10);
  /* looping over the nonzero elements of the shell */
  printf("computing the fractions\n");
  volFraction frac;
  octTreeVoxel V;
  double vox[3], vsize[3];
  vsize[0] = mri_src->xsize;
  vsize[1] = mri_src->ysize;
  vsize[2] = mri_src->zsize;
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (MRIgetVoxVal(mri_shell, x, y, z, 0) > 125.0) {
          /* change of coordinates from image to surface domain */
          MRIvoxelToSurfaceRAS(mri_shell, x, y, z, &xs, &ys, &zs);
          /* find the closest vertex to the point */
          MHTfindClosestVertexGeneric(mht, xs, ys, zs, 10, 2, &vno, &dist);
          /* creating the oct tree voxel structure */
          vox[0] = xs - vsize[0] / 2.0;
          vox[1] = ys - vsize[1] / 2.0;
          vox[2] = zs - vsize[2] / 2.0;
          V = octTreeVoxelCreate(vox, vsize);
          /* compute the volume fraction of this voxel */
          frac = MRIcomputeVoxelFractions(V, vno, acc, 1, mris);
          MRIsetVoxVal(mri_fractions, x, y, z, 0, frac.frac);
        }
        else if (MRIgetVoxVal(mri_interior, x, y, z, 0) > 0.0)
          MRIsetVoxVal(mri_fractions, x, y, z, 0, 1.0);
      }
    }
  }
  return mri_fractions;
}
volFraction MRIcomputeVoxelFractions(octTreeVoxel V, int vno, double acc, int current_depth, MRI_SURFACE *mris)
{
  /* inputs:
     V: voxel element, which includes center, corners and the vsize
     v: closest vertex to this voxel
     acc: the accuracy requirement - max error tolerable per voxel
     current_depth: current depth in the oct-tree structure
     mris: mri surface required for face normal information.
   */
   
  // HACK until the caller can provide the vno.  Can't right now.
  //
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
  
  int fnum, fiter, j, k;
  double meanNorm[3], meanVert[3];
  /* get the faces the closest vertex to the point (xs,ys,zs) is a part of */
  fnum = vt->num;
  /* compute a mean normal based on these faces */
  /* underlying assumption is that this mean normal is a good vector
     to decide if a point is inside or outside a surface
     CAUTION: there might very odd surfaces that can violate this assumption */
  meanNorm[0] = 0.0;
  meanNorm[1] = 0.0;
  meanNorm[2] = 0.0;
  for (fiter = 0; fiter < fnum; fiter++) {
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, vt->f[fiter]);
    meanNorm[0] = meanNorm[0] + fNorm->nx / (double)fnum;
    meanNorm[1] = meanNorm[1] + fNorm->ny / (double)fnum;
    meanNorm[2] = meanNorm[2] + fNorm->nz / (double)fnum;
  }
  meanVert[0] = v->x;
  meanVert[1] = v->y;
  meanVert[2] = v->z;
  /* detemining if the voxel is completely in or out */
  /* completely in means all the corners and the center of the voxel is inside the surface */
  int allin = 1;
  int allout = 1;
  double dotProduct;
  for (j = 0; j < 8; j++) {
    dotProduct = meanNorm[0] * (V.corn[j][0] - meanVert[0]) + meanNorm[1] * (V.corn[j][1] - meanVert[1]) +
                 meanNorm[2] * (V.corn[j][2] - meanVert[2]);
    if (dotProduct < 0.0) {
      allout = 0;
    }
    else {
      allin = 0;
    }
  }

  volFraction frac;
  frac.frac = 0;
  frac.err = 0;
  /* relative volume of the voxel to the whole */
  double relativeVolume =
      1.0 / (double)(pow(2.0, current_depth - 1) * pow(2.0, current_depth - 1) * pow(2.0, current_depth - 1));
  if (allin == 1) /* all corners are in return the relative volume of the voxel */
  {
    frac.frac = relativeVolume;
    frac.err = 0;
  }
  else if (allout == 1) /* all corners are out return 0 */
  {
    frac.frac = 0;
    frac.err = 0;
  }
  else /* some in some out */
  {
    if (relativeVolume < acc) /* the error we are making is small enough */
    {
      frac.frac = relativeVolume / 2;
      frac.err = relativeVolume / 2;
    }
    else /* the error we are making is too big we will redivide and recurse */
    {
      for (k = 0; k < 8; k++) {
        volFraction frac_new = MRIcomputeVoxelFractions(octTreeVoxelDivide(k + 1, V), vno, acc, current_depth + 1, mris);
        frac.frac += frac_new.frac;
        frac.err += frac_new.err;
      }
    }
  }
  return frac;
}
octTreeVoxel octTreeVoxelCreate(double *vox, double *vsize)
{
  octTreeVoxel v;
  int k;
  /*
  for (k = 0; k < 3; k++)
    { v.vox[k] = vox[k]; v.vsize[k] = vsize[k]; }
  */
  for (k = 0; k < 3; k++) {
    v.vox[k] = vox[k];
    v.vsize[k] = vsize[k];
  }
  /*center*/
  v.cent[0] = v.vox[0] + v.vsize[0] / 2.0;
  v.cent[1] = v.vox[1] + v.vsize[1] / 2.0;
  v.cent[2] = v.vox[2] + v.vsize[2] / 2.0;
  /* 000 */
  v.corn[0][0] = v.vox[0];
  v.corn[0][1] = v.vox[1];
  v.corn[0][2] = v.vox[2];
  /* 100 */
  v.corn[1][0] = v.vox[0] + v.vsize[0];
  v.corn[1][1] = v.vox[1];
  v.corn[1][2] = v.vox[2];
  /* 010 */
  v.corn[2][0] = v.vox[0];
  v.corn[2][1] = v.vox[1] + v.vsize[1];
  v.corn[2][2] = v.vox[2];
  /* 001 */
  v.corn[3][0] = v.vox[0];
  v.corn[3][1] = v.vox[1];
  v.corn[3][2] = v.vox[2] + v.vsize[2];
  /* 110 */
  v.corn[4][0] = v.vox[0] + v.vsize[0];
  v.corn[4][1] = v.vox[1] + v.vsize[1];
  v.corn[4][2] = v.vox[2];
  /* 011 */
  v.corn[5][0] = v.vox[0];
  v.corn[5][1] = v.vox[1] + v.vsize[1];
  v.corn[5][2] = v.vox[2] + v.vsize[2];
  /* 101 */
  v.corn[6][0] = v.vox[0] + v.vsize[0];
  v.corn[6][1] = v.vox[1];
  v.corn[6][2] = v.vox[2] + v.vsize[2];
  /* 111 */
  v.corn[7][0] = v.vox[0] + v.vsize[0];
  v.corn[7][1] = v.vox[1] + v.vsize[1];
  v.corn[7][2] = v.vox[2] + v.vsize[2];

  return v;
}
octTreeVoxel octTreeVoxelDivide(int type, octTreeVoxel v)
{
  double vsize_new[3];
  double vox_new[3] = {0, 0, 0};
  vsize_new[0] = v.vsize[0] / 2.0;
  vsize_new[1] = v.vsize[1] / 2.0;
  vsize_new[2] = v.vsize[2] / 2.0;
  switch (type) {
    case 1: /* 000 */
      vox_new[0] = v.vox[0];
      vox_new[1] = v.vox[1];
      vox_new[2] = v.vox[2];
      break;
    case 2: /* 100 */
      vox_new[0] = v.vox[0] + vsize_new[0];
      vox_new[1] = v.vox[1];
      vox_new[2] = v.vox[2];
      break;
    case 3: /* 010 */
      vox_new[0] = v.vox[0];
      vox_new[1] = v.vox[1] + vsize_new[1];
      vox_new[2] = v.vox[2];
      break;
    case 4: /* 001 */
      vox_new[0] = v.vox[0];
      vox_new[1] = v.vox[1];
      vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    case 5: /* 110 */
      vox_new[0] = v.vox[0] + vsize_new[0];
      vox_new[1] = v.vox[1] + vsize_new[1];
      vox_new[2] = v.vox[2];
      break;
    case 6: /* 101 */
      vox_new[0] = v.vox[0] + vsize_new[0];
      vox_new[1] = v.vox[1];
      vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    case 7: /* 011 */
      vox_new[0] = v.vox[0];
      vox_new[1] = v.vox[1] + vsize_new[1];
      vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    case 8: /* 111 */
      vox_new[0] = v.vox[0] + vsize_new[0];
      vox_new[1] = v.vox[1] + vsize_new[1];
      vox_new[2] = v.vox[2] + vsize_new[2];
      break;
    default:
      break;
  }
  return octTreeVoxelCreate(vox_new, vsize_new);
}
