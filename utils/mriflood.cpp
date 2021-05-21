/**
 * @brief flood and erode voxels
 *
 */
/*
 * Original Author: Andre van der Kouwe
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cma.h"
#include "diag.h"
#include "fsenv.h"
#include "macros.h"
#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_metricProperties.h"
#include "region.h"
#include "timer.h"

#define IMGSIZE 256

#define subvoxmask(i, j, k) (1 << ((i) + ((j) << 1) + ((k) << 2)))

#define NEIGHBOURHALFSIDE 2

int checkx = 152;
int checky = 127;
int checkz = 133;

unsigned char SubFill(unsigned char vox, int i, int j, int k);
MRI *MRIbitwiseand(MRI *mri1, MRI *mri2, MRI *mri_dst);
MRI *MRIbitwiseor(MRI *mri1, MRI *mri2, MRI *mri_dst);
MRI *MRIbitwisenot(MRI *mri_src, MRI *mri_dst);
MRI *MRImajority(MRI *mri_src, MRI *mri_dst);
MRI *MRImergecortexwhitecma(MRI *mri_cortex, MRI *mri_white, MRI *mri_cma, MRI *mri_left, MRI *mri_dst);
int HemisphereVote(MRI *mri_cma, int i, int j, int k, int halfside);
float distance(float x, float y, float z);
void MRIerodecerebralcortex(MRI *mri_masked, MRI *mri_cma, MRI *mri_white, MRI *mri_left);
int IllegalCorticalNeighbour(MRI *mri_masked, MRI *mri_white, int i, int j, int k);
void MRIcorrecthippocampus(MRI *mri_masked, MRI *mri_dst);

// return number of bits on ( possible values are 0 through 8 )
static int countBits(MRI *mri, int i, int j, int k)
{
  int nvox;
  nvox = (int)MRIgetVoxVal(mri, i, j, k, 0);
  return (((nvox >> 7) & 1) + ((nvox >> 6) & 1) + ((nvox >> 5) & 1) + ((nvox >> 4) & 1) + ((nvox >> 3) & 1) +
          ((nvox >> 2) & 1) + ((nvox >> 1) & 1) + (nvox & 1));
}
// return true if count of the right voxel > 4, else false
static int likely(MRI *mri, int i, int j, int k)
{
  if (countBits(mri, i, j, k) > 4)
    return 1;
  else
    return 0;
}

static void likelinessHistogram(MRI *mri, const char *msg)
{
  int i, j, k;
  long Hist[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (k = 0; k < mri->depth; ++k)
    for (j = 0; j < mri->height; ++j)
      for (i = 0; i < mri->width; ++i) Hist[countBits(mri, i, j, k)]++;

  printf("\nHistogram for %s\n", msg);
  for (i = 0; i < 9; ++i) printf(" %d : %ld\n", i, Hist[i]);
  printf("\n");
}

#ifndef __OPTIMIZE__
void DebugVoxel(char *msg, MRI *mri, int x, int y, int z)
{
  printf("================================================================\n");
  printf("%s (%d,%d,%d) = %d\n", msg, x, y, z, MRIvox(mri, x, y, z));
}
#endif

/* MRIribbon determines the space between the inner
   and outer MRI surfaces provided, */
/* and creates a volume in mri_dst corresponding to the input format mri_src */
MRI *MRISribbon(MRI_SURFACE *inner_mris, MRI_SURFACE *outer_mris, MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_inter;

  /* Allocate new MRI structures as needed */
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("MRISribbon: Creating new (_inter)MRI of %d, %d,%d...\n", mri_src->width, mri_src->height, mri_src->depth);
  mri_inter = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, mri_src->type);
  MRIcopyHeader(mri_src, mri_inter);
  if (!mri_dst) {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("MRISribbon: Creating new (_dst)MRI...\n");
    mri_dst = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Creating volume inside outer shell...\n");
  /* Create volume inside outer shell */
  /* Create shell corresponding to surface
     in MRI volume (includes outer shell in surface) */
  /////////////////////////////////////////////////////////
  // you should not combine MRISpartialshell() with MRIfloodoutside(), since
  // partialshell produces the shell value of 1 to 255.
  // However, floodoutside uses
  // 1 as the filled value, not shell value.  For example,
  //
  //      *   255  (MRISshell)        * 254  (MRISpartialshell)
  //     255   X                      1  X
  //
  // In the former, X is never set to 1,
  // but the latter X becomes 1 by floodoutside
  // Thus floodouside never becomes really floodoutside.
  //////////////////////////////////////////////////////////////
  MRISshell(mri_src, outer_mris, mri_inter, 1);
  MRISfloodoutside(mri_inter, mri_inter);
  MRISaccentuate(mri_inter, mri_inter, 1, 254);
  MRIcomplement(mri_inter, mri_inter);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Creating volume outside inner shell...\n");
  /* Create volume outside inner shell */
  MRISshell(mri_src, inner_mris, mri_dst, 1);
  MRISfloodoutside(mri_dst, mri_dst);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Finding intersection of volumes...\n");
  /* Find intersection of volumes to create ribbon */
  MRIintersect(mri_inter, mri_dst, mri_dst);
  MRISaccentuate(mri_dst, mri_dst, 1, 255);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Done with intersection...\n");
  MRIfree(&mri_inter);

  return mri_dst;
}

/* MRISshell needs mri_src as an example of
   the format for the output, to match size and type.  The */
/* surface is recreated in the MRI space (mri_dst)
   from the tesselated surface (mris) as voxels of 255 */
MRI *MRISshell(MRI *mri_src, MRI_SURFACE *mris, MRI *mri_dst, int clearflag)
{
  int width, height, depth, i, j, imnr, fno, numu, numv, u, v;
  // int imnr0;
  // float ps, st, xx0, xx1, yy0, yy1, zz0, zz1;
  float x0, y0, z0, x1, y1, z1, x2, y2, z2, d0, d1, d2, dmax;
  float px0, py0, pz0, px1, py1, pz1, px, py, pz;
  double fi, fj, fimnr;
  VERTEX *v_0, *v_1, *v_2;
  FACE *f;

  // imnr0 =  mri_src->imnr0;
  // st = mri_src->thick; /* slice thickness */
  // ps = mri_src->ps;
  // xx0 = mri_src->xstart;
  // xx1 = mri_src->xend;
  // yy0 = mri_src->ystart;
  // yy1 = mri_src->yend;
  // zz0 = mri_src->zstart;
  // zz1 = mri_src->zend;

  /* Create new blank MRI or clear existing destination MRI */
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    /*    printf("MRISshell: Creating new (_dst)MRI...\n");*/
    mri_dst = MRIalloc(width, height, depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (clearflag) MRIclear(mri_dst);

  /* Fill each face in MRI volume */
  for (fno = 0; fno < mris->nfaces; fno++) {
    /* Calculate (x,y,z) for each vertex for face */
    f = &mris->faces[fno];
    v_0 = &mris->vertices[f->v[0]];
    v_1 = &mris->vertices[f->v[1]];
    v_2 = &mris->vertices[f->v[2]];
    if (f->v[0] == Gdiag_no || f->v[1] == Gdiag_no || f->v[2] == Gdiag_no) DiagBreak();
    x0 = v_0->x;
    y0 = v_0->y;
    z0 = v_0->z;
    x1 = v_1->x;
    y1 = v_1->y;
    z1 = v_1->z;
    x2 = v_2->x;
    y2 = v_2->y;
    z2 = v_2->z;

    /* Calculate triangle side lengths */
    d0 = sqrt(SQR(x1 - x0) + SQR(y1 - y0) + SQR(z1 - z0)) / mri_dst->xsize;
    d1 = sqrt(SQR(x2 - x1) + SQR(y2 - y1) + SQR(z2 - z1)) / mri_dst->ysize;
    d2 = sqrt(SQR(x0 - x2) + SQR(y0 - y2) + SQR(z0 - z2)) / mri_dst->zsize;
    /* Divide space between sides into numv parallel lines */
    dmax = (d0 >= d1 && d0 >= d2) ? d0 : (d1 >= d0 && d1 >= d2) ? d1 : d2;
    numu = ceil(2 * d0);
    numv = ceil(2 * dmax);
    /* Fill each line in MRI volume */
    for (v = 0; v <= numv; v++) {
      px0 = x0 + (x2 - x0) * (float)v / (float)numv;
      py0 = y0 + (y2 - y0) * (float)v / (float)numv;
      pz0 = z0 + (z2 - z0) * (float)v / (float)numv;
      px1 = x1 + (x2 - x1) * (float)v / (float)numv;
      py1 = y1 + (y2 - y1) * (float)v / (float)numv;
      pz1 = z1 + (z2 - z1) * (float)v / (float)numv;
      /* Fill each voxel on line in MRI volume */
      for (u = 0; u <= numu; u++) {
        px = px0 + (px1 - px0) * (float)u / (float)numu;
        py = py0 + (py1 - py0) * (float)u / (float)numu;
        pz = pz0 + (pz1 - pz0) * (float)u / (float)numu;
/* Note mapping (x,y,z)<->(i,j,k) */
//        imnr = (int)((py-yy0)/st+1.5-imnr0);
//        i = (int)((xx1-px)/ps+0.5);
//        j = (int)((zz1-pz)/ps+1.0);
#if 0
        if (mris->useRealRAS)
          MRIworldToVoxel(mri_src,px,py,pz,&fi,&fj,&fimnr);
        else
          MRIsurfaceRASToVoxel(mri_src,px,py,pz,&fi,&fj,&fimnr);
#else
        MRISsurfaceRASToVoxelCached(mris, mri_src, px, py, pz, &fi, &fj, &fimnr);
#endif
        i = nint(fi);
        j = nint(fj);
        imnr = nint(fimnr);
        if (i >= 0 && i < mri_dst->width && j >= 0 && j < mri_dst->height && imnr >= 0 && imnr < mri_dst->depth)
          MRIsetVoxVal(mri_dst, i, j, imnr, 0, 255);
      }
    }
  }

  return mri_dst;
}

/* Floods MRI volume from outermost corners inward */
/* Fill with 1, boundary is anything but 0 and 1 */
// mri_src is just a dummy
MRI *MRISfloodoutside(MRI *mri_src, MRI *mri_dst)
{
  int newfilled, width, height, depth, i, j, k, v1, v2, v3;

  mri_dst = MRIcopy(mri_src, mri_dst);

  /* Set MRI size */
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

/* Set seed voxel in corner of box */
/*  MRIvox(mri_dst,1,1,1)=1;

newfilled=1;
while (newfilled>0)
{
newfilled=0;

for (i=1;i<width-1;i++)
for (j=1;j<height-1;j++)
for (k=1;k<depth-1;k++)
if ((int)MRIgetVoxVal(mri_dst,i,j,k,0)==0)
if ((int)MRIgetVoxVal(mri_dst,i,j,k-1,0)==1||
(int)MRIgetVoxVal(mri_dst,i-1,j,k,0)==1||
(int)MRIgetVoxVal(mri_dst,i,j-1,k)==1) {
(int)MRIgetVoxVal(mri_dst,i,j,k,0)=1;
newfilled++;
}
for (i=width-2;i>=1;i--)
for (j=height-2;j>=1;j--)
for (k=depth-2;k>=1;k--)
if ((int)MRIgetVoxVal(mri_dst,i,j,k,0)==0)
if ((int)MRIgetVoxVal(mri_dst,i,j,k+1,0)==1||
(int)MRIgetVoxVal(mri_dst,i+1,j,k,0)==1||
(int)MRIgetVoxVal(mri_dst,i,j+1,k,0)==1) {
(int)MRIgetVoxVal(mri_dst,i,j,k,0)=1;
newfilled++;
}
}*/

#define FILL_VAL 1
  MRIsetVoxVal(mri_dst, 0, 0, 0, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, width - 1, 0, 0, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, 0, height - 1, 0, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, 0, 0, depth - 1, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, width - 1, height - 1, 0, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, width - 1, 0, depth - 1, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, 0, height - 1, depth - 1, 0, FILL_VAL);
  MRIsetVoxVal(mri_dst, width - 1, height - 1, depth - 1, 0, FILL_VAL);

  newfilled = 1;
  while (newfilled > 0) {
    newfilled = 0;

    /* corner: 0 0 0 */
    for (i = 0; i < width; i++)
      for (j = 0; j < height; j++)
        for (k = 0; k < depth; k++) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i - 1 + ((i == 0) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j - 1 + ((j == 0) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k - 1 + ((k == 0) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: w 0 0 */
    for (i = width - 1; i >= 0; i--)
      for (j = 0; j < height; j++)
        for (k = 0; k < depth; k++) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i + 1 - ((i == width - 1) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j - 1 + ((j == 0) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k - 1 + ((k == 0) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: 0 h 0 */
    for (i = 0; i < width; i++)
      for (j = height - 1; j >= 0; j--)
        for (k = 0; k < depth; k++) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i - 1 + ((i == 0) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j + 1 - ((j == height - 1) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k - 1 + ((k == 0) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: 0 0 d */
    for (i = 0; i < width; i++)
      for (j = 0; j < height; j++)
        for (k = depth - 1; k >= 0; k--) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i - 1 + ((i == 0) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j - 1 + ((j == 0) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k + 1 - ((k == depth - 1) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: w h 0 */
    for (i = width - 1; i >= 0; i--)
      for (j = height - 1; j >= 0; j--)
        for (k = 0; k < depth; k++) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i + 1 - ((i == width - 1) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j + 1 - ((j == height - 1) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k - 1 + ((k == 0) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: w 0 d */
    for (i = width - 1; i >= 0; i--)
      for (j = 0; j < height; j++)
        for (k = depth - 1; k >= 0; k--) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i + 1 - ((i == width - 1) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j - 1 + ((j == 0) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k + 1 - ((k == depth - 1) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: 0 h d */
    for (i = 0; i < width; i++)
      for (j = height - 1; j >= 0; j--)
        for (k = depth - 1; k >= 0; k--) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i - 1 + ((i == 0) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j + 1 - ((j == height - 1) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k + 1 - ((k == depth - 1) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
    /* corner: w h d */
    for (i = width - 1; i >= 0; i--)
      for (j = height - 1; j >= 0; j--)
        for (k = depth - 1; k >= 0; k--) {
          if ((int)MRIgetVoxVal(mri_dst, i, j, k, 0) == 0) {
            v1 = (int)MRIgetVoxVal(mri_dst, i + 1 - ((i == width - 1) ? 1 : 0), j, k, 0);
            v2 = (int)MRIgetVoxVal(mri_dst, i, j + 1 - ((j == height - 1) ? 1 : 0), k, 0);
            v3 = (int)MRIgetVoxVal(mri_dst, i, j, k + 1 - ((k == depth - 1) ? 1 : 0), 0);
            if (v1 == FILL_VAL || v2 == FILL_VAL || v3 == FILL_VAL) {
              if (i == Gx && j == Gy && k == Gz) DiagBreak();
              MRIsetVoxVal(mri_dst, i, j, k, 0, FILL_VAL);
              newfilled++;
            }
          }
        }
  }

  return mri_dst;
}

MRI *MRISaccentuate(MRI *mri_src, MRI *mri_dst, int lo_thresh, int hi_thresh)
{
  int width, height, depth, i, j, k;
  float val;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    /*    printf("MRISaccentuate: Creating new (_dst)MRI...\n");*/
    mri_dst = MRIalloc(width, height, depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (k = 0; k < depth; k++)
    for (j = 0; j < height; j++)
      for (i = 0; i < width; i++) {
        int vox = (int)MRIgetVoxVal(mri_src, i, j, k, 0);
        val = ((vox >= lo_thresh) && (vox <= hi_thresh)) ? 255 : 0;
        MRIsetVoxVal(mri_dst, i, j, k, 0, val);
      }
  return mri_dst;
}

/* Set mri_dst voxel to 255 for every mri_src voxel for which the majority
   of subvoxels is set. */
MRI *MRImajority(MRI *mri_src, MRI *mri_dst)
{
  int width, height, depth, i, j, k, vox;
  long counts[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  float val;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (k = 0; k < depth; k++)
    for (j = 0; j < height; j++)
      for (i = 0; i < width; i++) {
        counts[countBits(mri_src, i, j, k)]++;
        vox = MRIgetVoxVal(mri_src, i, j, k, 0);
        val = ((((vox >> 7) & 1) + ((vox >> 6) & 1) + ((vox >> 5) & 1) + ((vox >> 4) & 1) + ((vox >> 3) & 1) +
                ((vox >> 2) & 1) + ((vox >> 1) & 1) + (vox & 1)) > 4)
                  ? 255
                  : 0;
        MRIsetVoxVal(mri_dst, i, j, k, 0, val);
      }
  printf("Inside the majority\n");
  for (i = 0; i < 9; ++i) printf(" %d : %ld\n", i, counts[i]);
  printf("\n");

  return mri_dst;
}

MRI *MRIbitwisenot(MRI *mri_src, MRI *mri_dst)
{
  int width, height, depth, i, j, k, vox;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (k = 0; k < depth; k++)
    for (j = 0; j < height; j++)
      for (i = 0; i < width; i++) {
        vox = (int)(MRIgetVoxVal(mri_src, i, j, k, 0));
        MRIsetVoxVal(mri_dst, i, j, k, 0, ~vox);
      }
  return mri_dst;
}

/* The following is the same as above, but adapted for 2 x 2 x 2 resolution */
/* per voxel. */

/* Partial shell fills the voxel space with a superresolution surface - */
/* each voxel is divided into 8 subvoxels represented by the 8 bits of the */
/* unsigned char in the MRI structure. */
MRI *MRISpartialshell(MRI *mri_src, MRI_SURFACE *mris, MRI *mri_dst, int clearflag)
{
  int width, height, depth, i, j, imnr, isub, jsub, isubmnr, fno, numu, numv, u, v;
  // int imnr0;
  // float ps, st, xx0, xx1, yy0, yy1, zz0, zz1;
  float x0, y0, z0, x1, y1, z1, x2, y2, z2, d0, d1, d2, dmax;
  float px0, py0, pz0, px1, py1, pz1, px, py, pz;
  double fi, fj, fimnr;
  VERTEX *v_0, *v_1, *v_2;
  FACE *f;
  int val;

  // imnr0 = mri_src->imnr0;
  // st = mri_src->thick; /* slice thickness */
  // ps = mri_src->ps;
  // xx0 = mri_src->xstart;
  // xx1 = mri_src->xend;
  // yy0 = mri_src->ystart;
  // yy1 = mri_src->yend;
  // zz0 = mri_src->zstart;
  // zz1 = mri_src->zend;

  /* Create new blank MRI or clear existing destination MRI */
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("MRISshell: Creating new (_dst)MRI...\n");
    mri_dst = MRIalloc(width, height, depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }
  if (clearflag) MRIclear(mri_dst);

  /* Fill each face in MRI volume */
  for (fno = 0; fno < mris->nfaces; fno++) {
    /* Calculate (x,y,z) for each vertex for face */
    f = &mris->faces[fno];
    v_0 = &mris->vertices[f->v[0]];
    v_1 = &mris->vertices[f->v[1]];
    v_2 = &mris->vertices[f->v[2]];
    x0 = v_0->x;
    y0 = v_0->y;
    z0 = v_0->z;
    x1 = v_1->x;
    y1 = v_1->y;
    z1 = v_1->z;
    x2 = v_2->x;
    y2 = v_2->y;
    z2 = v_2->z;

    /* Calculate triangle side lengths */
    d0 = sqrt(SQR(x1 - x0) + SQR(y1 - y0) + SQR(z1 - z0));
    d1 = sqrt(SQR(x2 - x1) + SQR(y2 - y1) + SQR(z2 - z1));
    d2 = sqrt(SQR(x0 - x2) + SQR(y0 - y2) + SQR(z0 - z2));
    /* Divide space between sides into numv parallel lines */
    dmax = (d0 >= d1 && d0 >= d2) ? d0 : (d1 >= d0 && d1 >= d2) ? d1 : d2;
    numu = ceil(2 * d0);
    numv = ceil(2 * dmax);
    /* Fill each line in MRI volume */
    for (v = 0; v <= numv * 2; v++) {
      px0 = x0 + (x2 - x0) * (float)v / (float)numv / 2;
      py0 = y0 + (y2 - y0) * (float)v / (float)numv / 2;
      pz0 = z0 + (z2 - z0) * (float)v / (float)numv / 2;
      px1 = x1 + (x2 - x1) * (float)v / (float)numv / 2;
      py1 = y1 + (y2 - y1) * (float)v / (float)numv / 2;
      pz1 = z1 + (z2 - z1) * (float)v / (float)numv / 2;
      /* Fill each voxel on line in MRI volume */
      for (u = 0; u <= numu * 2; u++) {
        px = px0 + (px1 - px0) * (float)u / (float)numu / 2;
        py = py0 + (py1 - py0) * (float)u / (float)numu / 2;
        pz = pz0 + (pz1 - pz0) * (float)u / (float)numu / 2;
        /* Note mapping (x,y,z)<->(i,j,k) */
        /* Increasing the offset of 1.5 shifts shell
           in the anterior direction */
        //        imnr = (int)((py-yy0)/st+1.5-imnr0);
        /* Increasing the offset of 0.5 shifts shell to the right */
        //        i = (int)((xx1-px)/ps+0.5);
        /* Increasing the offset of 1.0 shifts shell in
           the inferior direction */
        //        j = (int)((zz1-pz)/ps+1.0);
        // MRIworldToVoxel(mri_src,px,py,pz,&fi,&fj,&fimnr);
        MRISsurfaceRASToVoxelCached(mris, mri_src, px, py, pz, &fi, &fj, &fimnr);
        i = nint(fi);
        j = nint(fj);
        imnr = nint(fimnr);
        if (i >= 0 && i < IMGSIZE && j >= 0 && j < IMGSIZE && imnr >= 0 && imnr < depth) {
          /* Each voxel has 8 subvoxels, represented by
             the 8 bits of the unsigned char. */
          //          isubmnr = ((int)(((py-yy0)/st+1.5-imnr0)*2))-
          // ((int)((py-yy0)/st+1.5-imnr0))*2;
          //          isub = ((int)(((xx1-px)/ps+0.5)*2))-((int)((xx1-px)/ps+0.5))*2;
          //          jsub = ((int)(((zz1-pz)/ps+1.0)*2))-((int)((zz1-pz)/ps+1.0))*2;
          isubmnr = (int)((fimnr - nint(fimnr)) * 2 + 1);
          isub = (int)((fi - nint(fi)) * 2 + 1);
          jsub = (int)((fj - nint(fj)) * 2 + 1);
          /* (isubmnr, isub, jsub) should be in the range (0..1, 0..1, 0..1) */
          /* Assume that the initial value for all voxels is zero */
          val = (int)MRIgetVoxVal(mri_dst, i, j, imnr, 0) | subvoxmask(isub, jsub, isubmnr);
          MRIsetVoxVal(mri_dst, i, j, imnr, 0, val);
        }
      }
    }
  }

  return mri_dst;
}

MRI *MRIbitwiseor(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int width, height, depth, x, y, z;
  // BUFTYPE *p1, *p2, *pdst;
  BUFTYPE v1, v2;

  width = mri1->width;
  height = mri1->height;
  depth = mri1->depth;

  if (!mri_dst) mri_dst = MRIclone(mri1, NULL);

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      // pdst = &MRIvox(mri_dst, 0, y, z);
      // p1 = &MRIvox(mri1, 0, y, z);
      // p2 = &MRIvox(mri2, 0, y, z);
      for (x = 0; x < width; x++) {
        v1 = (int)MRIgetVoxVal(mri1, x, y, z, 0);
        v2 = (int)MRIgetVoxVal(mri2, x, y, z, 0);
        MRIsetVoxVal(mri_dst, x, y, z, 0, v1 | v2);
      }
    }
  }
  return (mri_dst);
}

MRI *MRIbitwiseand(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int width, height, depth, x, y, z;
  // BUFTYPE *p1, *p2, *pdst;
  BUFTYPE v1, v2;

  width = mri1->width;
  height = mri1->height;
  depth = mri1->depth;

  if (!mri_dst) mri_dst = MRIclone(mri1, NULL);

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      // pdst = &MRIvox(mri_dst, 0, y, z);
      // p1 = &MRIvox(mri1, 0, y, z);
      // p2 = &MRIvox(mri2, 0, y, z);
      for (x = 0; x < width; x++) {
        v1 = (int)MRIgetVoxVal(mri1, x, y, z, 0);
        v2 = (int)MRIgetVoxVal(mri2, x, y, z, 0);
        MRIsetVoxVal(mri_dst, x, y, z, 0, v1 & v2);
      }
    }
  }
  return (mri_dst);
}

/* Floods MRI volume from outermost corners inward. */
/* Upon return, mri_dst contains the filled volume NOT including the shell. */
/* mri_src and mri_dst cannot be the same volume! */
MRI *MRISpartialfloodoutside(MRI *mri_src, MRI *mri_dst)
{
  int newfilled, width, height, depth, i, j, k, is, js, ks, isub, jsub, ksub;
  int val;

  /* Set MRI size */
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  /* Set seed voxel in corner of voxel in corner of box */
  if (mri_dst == NULL) mri_dst = MRIcopy(mri_src, NULL);
  val = (int)MRIgetVoxVal(mri_dst, 1, 1, 1, 0) | subvoxmask(0, 0, 0);
  MRIsetVoxVal(mri_dst, 1, 1, 1, 0, val);

  newfilled = 1;
  while (newfilled > 0) {
    newfilled = 0;

    printf("    (left to right)\n");
    for (is = 2; is < 2 * width - 1; is++) {
      for (js = 2; js < 2 * height - 1; js++)
        for (ks = 2; ks < 2 * depth - 1; ks++) {
          i = is / 2;
          j = js / 2;
          k = ks / 2;
          isub = is % 2;
          jsub = js % 2;
          ksub = ks % 2;
          /*printf("i,j,k,isub,jsub,ksub =
            %d,%d,%d,%d,%d,%d\n",i,j,k,isub,jsub,ksub);
            printf("curr. vox =
            %d,%d\n",(MRIvox(mri_src,i,j,k)&subvoxmask(isub,jsub,ksub)),
            ((int)MRIgetVoxVal(mri_dst,i,j,k,0)&subvoxmask(isub,jsub,ksub)));
            printf("adj. voxels = %d,%d,%d\n",
            ((int)MRIgetVoxVal(mri_dst,i,j,(ks-1)/2,0)&subvoxmask(isub,jsub,(ks-1)%2)),
            ((int)MRIgetVoxVal(mri_dst,(is-1)/2,j,k,0)&subvoxmask((is-1)%2,jsub,ksub)),
            ((int)MRIgetVoxVal(mri_dst,(is-1)/2,j,k,0)&subvoxmask((is-1)%2,jsub,ksub)));*/
          if ((((int)MRIgetVoxVal(mri_src, i, j, k, 0) & subvoxmask(isub, jsub, ksub)) == 0) &&
              (((int)MRIgetVoxVal(mri_dst, i, j, k, 0) & subvoxmask(isub, jsub, ksub)) == 0)) {
            if ((((int)MRIgetVoxVal(mri_dst, i, j, (ks - 1) / 2, 0) & subvoxmask(isub, jsub, (ks - 1) % 2)) > 0) ||
                (((int)MRIgetVoxVal(mri_dst, (is - 1) / 2, j, k, 0) & subvoxmask((is - 1) % 2, jsub, ksub)) > 0) ||
                (((int)MRIgetVoxVal(mri_dst, i, (js - 1) / 2, k, 0) & subvoxmask(isub, (js - 1) % 2, ksub)) > 0)) {
              val = (int)MRIgetVoxVal(mri_dst, i, j, k, 0) | subvoxmask(isub, jsub, ksub);
              MRIsetVoxVal(mri_dst, i, j, k, 0, val);
              newfilled++;
            }
          }
        }
    }
    printf("    (right to left)\n");
    for (is = 2 * width - 2; is >= 1; is--) {
      for (js = 2 * height - 2; js >= 1; js--)
        for (ks = 2 * depth - 2; ks >= 1; ks--) {
          i = is / 2;
          j = js / 2;
          k = ks / 2;
          isub = is % 2;
          jsub = js % 2;
          ksub = ks % 2;
          if ((((int)MRIgetVoxVal(mri_src, i, j, k, 0) & subvoxmask(isub, jsub, ksub)) == 0) &&
              (((int)MRIgetVoxVal(mri_dst, i, j, k, 0) & subvoxmask(isub, jsub, ksub)) == 0)) {
            if ((((int)MRIgetVoxVal(mri_dst, i, j, (ks + 1) / 2, 0) & subvoxmask(isub, jsub, (ks + 1) % 2)) > 0) ||
                (((int)MRIgetVoxVal(mri_dst, (is + 1) / 2, j, k, 0) & subvoxmask((is + 1) % 2, jsub, ksub)) > 0) ||
                (((int)MRIgetVoxVal(mri_dst, i, (js + 1) / 2, k, 0) & subvoxmask(isub, (js + 1) % 2, ksub)) > 0)) {
              val = (int)MRIgetVoxVal(mri_dst, i, j, k, 0) | subvoxmask(isub, jsub, ksub);
              MRIsetVoxVal(mri_dst, i, j, k, 0, val);
              newfilled++;
            }
          }
        }
    }
    printf("    (filled %d voxels)\n", newfilled);
  }
  // so far we have touched the faces of the volume partially.
  // we fill them with 255 (there are 6 faces)
  val = MRIgetVoxVal(mri_dst, 1, 1, 1, 0);
  for (is = 0; is < width; ++is)
    for (js = 0; js < height; ++js) {
      MRIsetVoxVal(mri_dst, is, js, 0, 0, val);
      MRIsetVoxVal(mri_dst, is, js, depth - 1, 0, val);
    }
  for (is = 0; is < width; ++is)
    for (ks = 0; ks < depth; ++ks) {
      MRIsetVoxVal(mri_dst, is, 0, ks, 0, val);
      MRIsetVoxVal(mri_dst, is, height - 1, ks, 0, val);
    }
  for (ks = 0; ks < depth; ++ks)
    for (js = 0; js < height; ++js) {
      MRIsetVoxVal(mri_dst, 0, js, ks, 0, val);
      MRIsetVoxVal(mri_dst, width - 1, js, ks, 0, val);
    }

  return mri_dst;
}

MRI *MRISpartialribbon(MRI_SURFACE *inner_mris_lh,
                       MRI_SURFACE *outer_mris_lh,
                       MRI_SURFACE *inner_mris_rh,
                       MRI_SURFACE *outer_mris_rh,
                       MRI *mri_src,
                       MRI *mri_dst,
                       MRI *mri_mask)
{
  MRI *mri_inter1, *mri_inter2, *mri_inter3;

  /* Allocate new MRI structures as needed */
  // MRIalloc uses "calloc" which sets the memory region to be filled with zero
  mri_inter1 = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, mri_src->type);
  MRIcopyHeader(mri_src, mri_inter1);
  mri_inter2 = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, mri_src->type);
  MRIcopyHeader(mri_src, mri_inter2);
  mri_inter3 = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, mri_src->type);
  MRIcopyHeader(mri_src, mri_inter3);
  if (!mri_dst) {
    mri_dst = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, mri_src->type);
    MRIcopyHeader(mri_src, mri_dst);
  }

  printf("Creating partial volume inside outer shell...\n");
  /* Create volume inside outer shell */
  /* Create shell corresponding to
     surface in MRI volume (includes outer shell in surface) */
  MRISpartialshell(mri_src, outer_mris_lh, mri_inter1, 1); /* partial shell
                                                              in mri_inter1 */
  MRISpartialshell(mri_src, outer_mris_rh, mri_inter1, 0); /* partial shell
                                                              in mri_inter1 */
  // so far inter2 not used and thus filled with 0
  MRISpartialfloodoutside(mri_inter1, mri_inter2); /* flooded outside shell
                                                      in mri_inter2 */
  MRIbitwisenot(mri_inter2, mri_inter2);           /* flooded inside shell and shell
                                                      in mri_inter2 */

  printf("Creating partial volume outside inner shell...\n");
  /* Create volume outside inner shell */
  MRISpartialshell(mri_src, inner_mris_lh, mri_inter1, 1); /* 1 => clear first */
  MRISpartialshell(mri_src, inner_mris_rh, mri_inter1, 0);
  // so far dst not used and thus filled with 0
  MRISpartialfloodoutside(mri_inter1, mri_dst);
  // save this results in inter3
  MRIcopy(mri_dst, mri_inter3);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("  - finding union of inner shell and outward-filled volume...\n");
  MRIbitwiseor(mri_inter1, mri_dst, mri_dst);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf(
        "Finding intersection of outershell-inside "
        "volume and innershell-outside volume\n");
  /* Find bitwise and of volumes to create ribbon */
  MRIbitwiseand(mri_inter2, mri_dst, mri_dst);
  /* Clear up the partial edges. */
  // printf("Converting volume to original resolution...\n");
  // MRImajority(mri_dst,mri_dst);
  // printf("majority cortex results =========================\n");
  // likelinessHistogram(mri_dst, "majority cortex");

  /* The problem with the excessive white matter must
     be within this if statement.!!!!! */
  /* If masked, change to CMA labels, add white
     matter label and apply to mask. */
  if (mri_mask) {
    printf("Using CMA labels to fine tune results\n");
    //
    // change made by tosa
    //
    // We decided to use the partial shell for better classification
    // the original routine uses just the full shell to do the cma comparison
    // By doing this, we save the separate calculation time for the full shell
    // and reuse the full shell calculation
    //
    // original routines
    // printf("Creating full volume outside inner shells...\n");
    // MRISshell(mri_src,inner_mris_lh,mri_inter1,1); // must clear
    // MRISshell(mri_src,inner_mris_rh,mri_inter1,0);
    // MRISfloodoutside(mri_inter1,mri_inter1); // src is just a dummy
    // printf("  - inverting volume...\n");
    // MRISaccentuate(mri_inter1,mri_inter2,1,254); // we really need this
    // MRIbitwisenot(mri_inter2,mri_inter2);
    // printf("fullvoxel white results\n");
    // likelinessHistogram(mri_inter2, "full voxel white matter");

    /* Create white matter volume in mri_inter1, including shell */
    // printf("Creating full volume outside inner shells (subvoxels...) **************\n");
    //////////////////////////////////////////////////////////////////////////
    // we saved the previous calculation in inter3
    // MRISpartialshell(mri_src,inner_mris_lh,mri_inter1,1);
    // 1=clear, lh white matter surface
    // MRISpartialshell(mri_src,inner_mris_rh,mri_inter1,0);
    // rh white matter surface
    // MRIclear(mri_inter2); // clear (needed)
    // MRISpartialfloodoutside(mri_inter1,mri_inter2);
    // src dst must be different
    // printf("  - inverting volume...\n");
    printf("Creating the inner shell volume\n");
    MRIclear(mri_inter2);
    MRIbitwisenot(mri_inter3, mri_inter2);

    /* Create volume inside left outer surface as reference,
       mri_inter3 contains left reference. */
    printf("Creating full reference 'left' volume...\n");
    MRISshell(mri_src, outer_mris_lh, mri_inter3, 1);  // clear flag
    MRISfloodoutside(mri_inter3, mri_inter3);          // src is dummy
    // printf("  - inverting volume...\n");
    MRISaccentuate(mri_inter3, mri_inter3, 1, 254);  // we really need this
    MRIbitwisenot(mri_inter3, mri_inter3);

    /* mri_dst contains cerebral cortex,
       mri_inter1 contains left side voxels,
       mri_inter2 contains white matter and some of the gray inner shell */

    MRIcopy(mri_dst, mri_inter1);
    MRIclear(mri_dst);  // clear

    // printf("Here are the values set for (%d,%d,%d)\n", checkx, checky, checkz);
    // DebugVoxel("before merge: cortex: ", mri_inter1, checkx, checky, checkz);
    // DebugVoxel("            : white : ", mri_inter2, checkx, checky, checkz);
    // DebugVoxel("            : lh    : ", mri_inter3, checkx, checky, checkz);
    // DebugVoxel("            : cma   : ", mri_mask, checkx, checky, checkz);
    printf("Initial classification of voxel likeliness\n");
    likelinessHistogram(mri_inter1, "cortex");
    likelinessHistogram(mri_inter2, "white matter");
    // likelinessHistogram(mri_inter3, "lh");
    //              cortex   white      cma         lh     labeled out
    printf("Merging labelled volumes...\n");
    MRImergecortexwhitecma(mri_inter1, mri_inter2, mri_mask, mri_inter3, mri_dst);
  }

  MRIclear(mri_inter1);
  MRIcopy(mri_dst, mri_inter1);
  // DebugVoxel("after merge", mri_inter1, checkx, checky, checkz);

  // note that mri_mask can be null
  if (mri_mask) {
    printf("Illegality check on cortex voxels...\n");
    MRIerodecerebralcortex(mri_inter1, mri_mask, mri_inter2, mri_inter3);
    // DebugVoxel("after erode", mri_inter1, checkx, checky, checkz);
  }
  printf("Illegality check on cortex voxels near hippocampus...\n");
  MRIcorrecthippocampus(mri_inter1, mri_dst);
  // DebugVoxel("after hippo", mri_dst, checkx, checky, checkz);

  MRIfree(&mri_inter1);
  MRIfree(&mri_inter2);
  MRIfree(&mri_inter3);

  return mri_dst;
}

MRI *MRImergecortexwhitecma(MRI *mri_cortex, MRI *mri_white, MRI *mri_cma, MRI *mri_left, MRI *mri_dst)
{
  /* mri_cortex = cerebral cortex is labelled as 255, all else is 0.
     mri_white  = white matter and some of the inner
     gray matter shell is labelled as 255, all else is 0.
     mri_cma    = CMA labelled volume.

     This function takes the mri_cma volume and replaces all instances of:
       Left_Cerebral_Cortex;
       Left_Cerebral_White_Matter;
       Right_Cerebral_Cortex;
       Right_Cerebral_White_Matter;
       Unknown,
     with:
       Left/Right_Cerebral_Cortex if labelled in mri_cortex;
       Left/Right_Cerebral_White_Matter if labelled
       in mri_white (and not also labelled in mri_cortex),
       where left/right is given by the CMA label
       side or the nearest neighbour CMA label vote;
       Unknown otherwise.
  */

  int width, height, depth, i, j, k, vox;
  int countunknownwhite, countunknowncortex, countunknownunknown;
  int countcortexunknown, countwhiteunknown;
  int countBitsWhite;
  int countBitsCortex;
  int likelyCortex;
  int likelyWhite;
  countunknownwhite = 0;
  countunknowncortex = 0;
  countunknownunknown = 0;
  countcortexunknown = 0;
  countwhiteunknown = 0;
  width = mri_cma->width;
  height = mri_cma->height;
  depth = mri_cma->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, mri_cma->type);
    MRIcopyHeader(mri_cma, mri_dst);
  }

  for (k = 0; k < depth; k++)
    for (j = 0; j < height; j++)
      for (i = 0; i < width; i++) {
        vox = (int)MRIgetVoxVal(mri_cma, i, j, k, 0);
        // first set the values, using cma
        MRIsetVoxVal(mri_dst, i, j, k, 0, vox);
        // cache the values
        countBitsCortex = countBits(mri_cortex, i, j, k);
        countBitsWhite = countBits(mri_white, i, j, k);
        likelyCortex = likely(mri_cortex, i, j, k);
        likelyWhite = likely(mri_white, i, j, k);
        ///////////////////////////////////////////////////////////
        if ((vox == Left_Cerebral_Cortex) || (vox == Left_Cerebral_White_Matter)) {
          if (likelyWhite && (countBitsWhite >= countBitsCortex))
            MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_White_Matter);
          else if (likelyCortex)
            MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_Cortex);
          else {
            MRIsetVoxVal(mri_dst, i, j, k, 0, Unknown);
            if (vox == Left_Cerebral_Cortex)
              countcortexunknown++;
            else
              countwhiteunknown++;
          }
        }
        else if ((vox == Right_Cerebral_Cortex) || (vox == Right_Cerebral_White_Matter)) {
          if (likelyWhite && (countBitsWhite >= countBitsCortex))
            MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_White_Matter);
          else if (likelyCortex)
            MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_Cortex);
          else {
            MRIsetVoxVal(mri_dst, i, j, k, 0, Unknown);
            if (vox == Right_Cerebral_Cortex)
              countcortexunknown++;
            else
              countwhiteunknown++;
          }
        }
        else if (vox == Unknown) {
          // cases are
          //              countBitsWhite >= countBitsCortex >=0
          //              countBitsCortex > countBitsWhite >=0
          // But we make sure that the voxel
          // must be likely one (i.e. > 4 bits set)
          //
          // is it white ?  more likely to be white than cortex?
          // if white and cortex has the same votes, then pick white (=)
          //
          // this cover the case of 1. likelyWhite >= likelyCortex
          //                        2. likelyWhite only
          if (likelyWhite && (countBitsWhite >= countBitsCortex)) {
            switch (HemisphereVote(mri_cma, i, j, k, NEIGHBOURHALFSIDE)) {
              case -1:
                // if (MRIvox(mri_left,i,j,k)==255)
                if (likely(mri_left, i, j, k))
                  MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_White_Matter);
                else
                  MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_White_Matter);
                break;
              case 0:
                MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_White_Matter);
                break;
              case 1:
                MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_White_Matter);
                break;
            }
            countunknownwhite++;
          }
          // we checked white, is it cortex then?
          // this covers the case of 1. likelyCortex > likelyWhite
          //                         2. likelyCortex only
          else if (likelyCortex) {
            // hemisphere vote returns only 1 (left) or 0 (right)
            switch (HemisphereVote(mri_cma, i, j, k, NEIGHBOURHALFSIDE)) {
              case -1:  // left = 0 and right = 0
                // if (MRIvox(mri_left,i,j,k)==255)
                if (likely(mri_left, i, j, k))
                  MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_Cortex);
                else
                  MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_Cortex);
                break;
              case 0:
                MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_Cortex);
                break;
              case 1:
                MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_Cortex);
                break;
            }
            countunknowncortex++;
          }
          else
            countunknownunknown++;
          // we don't change unknown state
        }  // vox unknown
      }    // loop

  printf("After merge\n");
  printf("\tcma cortex  changed to unknown: %8d\n", countcortexunknown);
  printf("\tcma white   changed to unknown: %8d\n", countwhiteunknown);
  printf("\tcma unknown changed to white  : %8d\n", countunknownwhite);
  printf("\tcma unknown changed to cortex : %8d\n", countunknowncortex);
  printf("\tcma unknwon remained unknown  : %8d\n", countunknownunknown);

  return mri_dst;
}

/* Return 1 for left,
   0 for right (searched a cube of sidelength halfside*2+1). */
int HemisphereVote(MRI *mri_cma, int i, int j, int k, int halfside)
{
  int x, y, z, vox;
  float leftvote, rightvote;
  int width, height, depth;

  width = mri_cma->width;
  height = mri_cma->height;
  depth = mri_cma->depth;

  leftvote = 0.;
  rightvote = 0.;

  for (x = i - halfside; x <= i + halfside; x++)
    for (y = j - halfside; y <= j + halfside; y++)
      for (z = k - halfside; z <= k + halfside; z++) {
        if ((x > 0) && (y > 0) && (z > 0) && (x < width) && (y < height) && (z < depth)) {
          vox = (int)MRIgetVoxVal(mri_cma, x, y, z, 0);
          if ((vox == Left_Cerebral_Cortex) || (vox == Left_Cerebral_White_Matter) || (vox == Left_Cerebral_Exterior) ||
              (vox == Left_Lateral_Ventricle) || (vox == Left_Inf_Lat_Vent) || (vox == Left_Cerebellum_Exterior) ||
              (vox == Left_Cerebellum_White_Matter) || (vox == Left_Cerebellum_Cortex) || (vox == Left_Thalamus) ||
              (vox == Left_Caudate) || (vox == Left_Putamen) ||
              (vox == Left_Pallidum) || (vox == Left_Hippocampus) || (vox == Left_Amygdala) || (vox == Left_Insula) ||
              (vox == Left_Operculum) || (vox == Left_Lesion) || (vox == Left_Accumbens_area) ||
              (vox == Left_Substancia_Nigra) || (vox == Left_VentralDC) || (vox == Left_undetermined) ||
              (vox == Left_vessel) || (vox == Left_choroid_plexus) || (vox == Left_F3orb) || (vox == Left_lOg) ||
              (vox == Left_aOg) || (vox == Left_mOg) || (vox == Left_pOg) || (vox == Left_Stellate) ||
              (vox == Left_Porg) || (vox == Left_Aorg))
            leftvote += 1 / distance((float)(x - i), (float)(y - j), (float)(z - k));
          if ((vox == Right_Cerebral_Cortex) || (vox == Right_Cerebral_White_Matter) ||
              (vox == Right_Cerebral_Exterior) || (vox == Right_Lateral_Ventricle) || (vox == Right_Inf_Lat_Vent) ||
              (vox == Right_Cerebellum_Exterior) || (vox == Right_Cerebellum_White_Matter) ||
              (vox == Right_Cerebellum_Cortex) || (vox == Right_Thalamus) ||
              (vox == Right_Caudate) || (vox == Right_Putamen) || (vox == Right_Pallidum) ||
              (vox == Right_Hippocampus) || (vox == Right_Amygdala) || (vox == Right_Insula) ||
              (vox == Right_Operculum) || (vox == Right_Lesion) || (vox == Right_Accumbens_area) ||
              (vox == Right_Substancia_Nigra) || (vox == Right_VentralDC) || (vox == Right_undetermined) ||
              (vox == Right_vessel) || (vox == Right_choroid_plexus) || (vox == Right_F3orb) || (vox == Right_lOg) ||
              (vox == Right_aOg) || (vox == Right_mOg) || (vox == Right_pOg) || (vox == Right_Stellate) ||
              (vox == Right_Porg) || (vox == Right_Aorg))
            rightvote += 1 / distance((float)(x - i), (float)(y - j), (float)(z - k));
        }
      }
  if ((leftvote == rightvote) && (leftvote == 0)) return -1;
  if (leftvote == rightvote) printf("Ambiguous voxel (%d, %d, %d) labelled right (%1.2f votes).\n", i, j, k, leftvote);

  return leftvote > rightvote;
}

float distance(float x, float y, float z) { return sqrt(x * x + y * y + z * z); }

void MRIerodecerebralcortex(MRI *mri_masked, MRI *mri_cma, MRI *mri_white, MRI *mri_left)
{
  int width, height, depth, i, j, k, vox, erodedvoxelcount, olderodedvoxelcount, unknowncount;
  int erodewhitecount;
  int erodecortexcma;

  width = mri_cma->width;
  height = mri_cma->height;
  depth = mri_cma->depth;

  olderodedvoxelcount = 0;
  erodecortexcma = 0;
  erodedvoxelcount = -1;
  unknowncount = 0;
  erodewhitecount = 0;
  // repreat danger
  while ((erodedvoxelcount != 0) && (erodedvoxelcount != olderodedvoxelcount)) {
    olderodedvoxelcount = erodedvoxelcount;
    erodedvoxelcount = 0;
    for (k = 0; k < depth; k++)
      for (j = 0; j < height; j++)
        for (i = 0; i < width; i++) {
          vox = MRIgetVoxVal(mri_masked, i, j, k, 0);
          /* If this voxel is not cerebral cortex,
             copy it directly to the destination volume */
          if ((vox == Left_Cerebral_Cortex) || (vox == Right_Cerebral_Cortex)) {
            if (IllegalCorticalNeighbour(mri_cma, mri_white, i, j, k)) {
              if ((int)MRIgetVoxVal(mri_cma, i, j, k, 0) != Unknown) {
                MRIsetVoxVal(mri_masked, i, j, k, 0, MRIgetVoxVal(mri_cma, i, j, k, 0));
                erodecortexcma++;  // illegal and thus take the cma value again
              }
              else { /* if the voxel needs to be eroded,
                        but the CMA didn't label it,
                        check if it's in the white matter volume. */
                // if (MRIvox(mri_white,i,j,k)==255)
                if (likely(mri_white, i, j, k)) {
                  switch (HemisphereVote(mri_cma, i, j, k, NEIGHBOURHALFSIDE)) {
                    case -1:
                      // if ((int)MRIgetVoxVal(mri_left,i,j,k,0)==255)
                      if (likely(mri_left, i, j, k))
                        MRIsetVoxVal(mri_masked, i, j, k, 0, Left_Cerebral_White_Matter);
                      else
                        MRIsetVoxVal(mri_masked, i, j, k, 0, Right_Cerebral_White_Matter);
                      break;
                    case 0:
                      MRIsetVoxVal(mri_masked, i, j, k, 0, Right_Cerebral_White_Matter);
                      break;
                    case 1:
                      MRIsetVoxVal(mri_masked, i, j, k, 0, Left_Cerebral_White_Matter);
                      break;
                  }
                  // printf("labelled as white (%d,%d,%d)\n", i,j,k);
                  erodewhitecount++;
                }
                else {
                  // printf("Voxel labelled as cortex,
                  // not in white matter volume (%d,%d,%d)\n",i,j,k);
                  MRIsetVoxVal(mri_masked, i, j, k, 0, Unknown);
                  unknowncount++;
                }
              }
              erodedvoxelcount++; /* only if value changed!
                                     then don't need old value */
            }
          }
        }
    printf("After illegality check on cortex labelled voxels\n");
    printf("\ttotal illegal voxels   : %8d\n", erodedvoxelcount);
  }
  printf("\tcortex became cma value: %8d\n", erodecortexcma);
  printf("\tcma unknown voxels changed as\n");
  printf("\tcortex became white    : %8d\n", erodewhitecount);
  printf("\tcortex became unknown  : %8d\n", unknowncount);
}

int IllegalCorticalNeighbour(MRI *mri_masked, MRI *mri_white, int i, int j, int k)
{
  int width, height, depth, x, y, z, vox, illegalflag;
  int nvox, ii, jj, kk;
  width = mri_masked->width;
  height = mri_masked->height;
  depth = mri_masked->depth;

  illegalflag = 0;
  for (x = i - 2; x <= i + 2; x++)
    for (y = j - 2; y <= j + 2; y++)
      for (z = k - 2; z <= k + 2; z++) {
        if ((x > 0) && (y > 0) && (z > 0) && (x < width) && (y < height) && (z < depth)) {
          vox = (int)MRIgetVoxVal(mri_masked, x, y, z, 0);
          /* caudate: Left_Caudate; Right_Caudate
             lateral ventricles:
             Left_Lateral_Ventricle; Right_Lateral_Ventricle
             inferior lateral ventricle: Left_Inf_Lat_Vent; Right_Inf_Lat_Vent
             thalamus: Left_Thalamus;
             Right_Thalamus */
          if ((vox == Left_Lateral_Ventricle) || (vox == Right_Lateral_Ventricle) || (vox == Left_Inf_Lat_Vent) ||
              (vox == Right_Inf_Lat_Vent)) {
            int dfx = abs(i - x);
            int dfy = abs(j - y);
            int dfz = abs(k - z);
            // there is a chance that there is a white matter in-between
            // Get the in-between position closest to the
            // ventricle and the voxel to check.
            // There is a choice to make, since the choice
            // can be A or B.  I pick A
            //                   V  A  *
            //                   *  B  x
            //
            if (dfx > 1 || dfy > 1 || dfz) {
              if (dfx > 1)
                ii = ((i - x) > 0) ? (x + 1) : (x - 1);
              else
                ii = x;
              if (dfy > 1)
                jj = ((j - y) > 0) ? (y + 1) : (y - 1);
              else
                jj = y;
              if (dfz > 1)
                kk = ((k - z) > 0) ? (z + 1) : (z - 1);
              else
                kk = z;
              nvox = (int)MRIgetVoxVal(mri_masked, ii, jj, kk, 0);
              if ((nvox == Left_Cerebral_White_Matter) || (nvox == Right_Cerebral_White_Matter)) {
                // printf("A WM (%d,%d,%d) between
                // (%d,%d,%d) ventricle and (%d,%d,%d) cortical\n",
                //       ii,jj,kk, x,y,z, i, j, k);
                continue;
              }
              else if (nvox == Unknown)  // CMA labelled as unknown,
                                         // then check if it is in white volume
              {
                if (likely(mri_white, ii, jj, kk))
                  // nvox = MRIvox(mri_white, ii, jj, kk);
                  // if (nvox == 255) // it is marked
                  // as white in white only volume
                  continue;
                else
                  illegalflag++;
              }
            }
            else
              illegalflag++;
          }
          else if ((vox == Left_Caudate) || (vox == Right_Caudate) || (vox == Left_Thalamus) ||
                   (vox == Right_Thalamus)) {
            illegalflag++;
          }
          // debug illegal check voxel values for voxels near i,j,k
          // if ((i==checkx) && (j==checky) && (k==checkz))
          //  printf("illegal check at (%d,%d,%d) vox = %d\n", x, y, z, vox);
        }
      }

  // the total voxel areas are 5x5x5 = 125
  return ((illegalflag > 0) ? 1 : 0);
}

void MRIcorrecthippocampus(MRI *mri_masked, MRI *mri_dst)
{
  /* mri_dst must differ from mri_masked. */
  int width, height, depth, i, j, k, vox, hippocount;

  width = mri_masked->width;
  height = mri_masked->height;
  depth = mri_masked->depth;
  hippocount = 0;

  for (k = 0; k < depth; k++)
    for (j = 0; j < height; j++)
      for (i = 0; i < width; i++) {
        vox = (int)MRIgetVoxVal(mri_masked, i, j, k, 0);
        MRIsetVoxVal(mri_dst, i, j, k, 0, vox);
        /* If this voxel is not cerebral cortex,
           copy it directly to the destination volume */
        if ((vox == Left_Cerebral_Cortex) || (vox == Right_Cerebral_Cortex)) {
          if ((((int)MRIgetVoxVal(mri_masked, i, j + 1, k, 0) == Left_Hippocampus) ||
               ((int)MRIgetVoxVal(mri_masked, i, j + 1, k, 0) == Right_Hippocampus) ||
               ((int)MRIgetVoxVal(mri_masked, i, j + 1, k, 0) == Left_Cerebral_White_Matter) ||
               ((int)MRIgetVoxVal(mri_masked, i, j + 1, k, 0) == Right_Cerebral_White_Matter)
               ////////////////////////////////////////////////////
               || ((int)MRIgetVoxVal(mri_masked, i, j + 2, k, 0) == Left_Hippocampus) ||
               ((int)MRIgetVoxVal(mri_masked, i, j + 2, k, 0) == Right_Hippocampus) ||
               ((int)MRIgetVoxVal(mri_masked, i, j + 2, k, 0) == Left_Cerebral_White_Matter) ||
               ((int)MRIgetVoxVal(mri_masked, i, j + 2, k, 0) == Right_Cerebral_White_Matter)) &&
              (((int)MRIgetVoxVal(mri_masked, i, j - 1, k, 0) == Left_Hippocampus) ||
               (MRIgetVoxVal(mri_masked, i, j - 1, k, 0) == Right_Hippocampus)
               ////////////////////////////////////////////////////
               || ((int)MRIgetVoxVal(mri_masked, i, j - 2, k, 0) == Left_Hippocampus) ||
               ((int)MRIgetVoxVal(mri_masked, i, j - 2, k, 0) == Right_Hippocampus))) {
            if (vox == Left_Cerebral_Cortex)
              MRIsetVoxVal(mri_dst, i, j, k, 0, Left_Cerebral_White_Matter);
            else
              MRIsetVoxVal(mri_dst, i, j, k, 0, Right_Cerebral_White_Matter);
            hippocount++;
          }
        }
      }
  printf("After hippocampus check\n");
  printf("\tcortex became white : %8d\n", hippocount);
}

MRI *MRISfillInteriorOld(MRI_SURFACE *mris, double resolution, MRI *mri_interior)
{
  int width, height, depth, x, y, z, val, saved_use_Real_RAS, interior_alloced;
  MATRIX *m_vox2ras;
  MRI *mri_shell, *mri_outside;

  saved_use_Real_RAS = mris->useRealRAS;
  mris->useRealRAS = 1;  // MRISshell needs this

  MRIScomputeMetricProperties(mris);

  if (mri_interior) {
    interior_alloced = 0;
    MRIclear(mri_interior);
    m_vox2ras = MRIgetVoxelToRasXform(mri_interior);
    mri_shell = MRIcloneDifferentType(mri_interior, MRI_FLOAT);
  }
  else {
    interior_alloced = 1;
    width = ceil((mris->xhi - mris->xlo) / resolution);
    height = ceil((mris->yhi - mris->ylo) / resolution);
    depth = ceil((mris->zhi - mris->zlo) / resolution);

    mri_shell = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIsetResolution(mri_shell, resolution, resolution, resolution);

    m_vox2ras = MatrixIdentity(4, NULL);
    *MATRIX_RELT(m_vox2ras, 1, 1) = resolution;
    *MATRIX_RELT(m_vox2ras, 2, 2) = resolution;
    *MATRIX_RELT(m_vox2ras, 3, 3) = resolution;

    *MATRIX_RELT(m_vox2ras, 1, 4) = mris->xlo + mris->vg.c_r;
    *MATRIX_RELT(m_vox2ras, 2, 4) = mris->ylo + mris->vg.c_a;
    *MATRIX_RELT(m_vox2ras, 3, 4) = mris->zlo + mris->vg.c_s;

    MRIsetVoxelToRasXform(mri_shell, m_vox2ras);
    mri_interior = MRIclone(mri_shell, NULL);
  }
  MRISshell(mri_shell, mris, mri_shell, 1);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_shell, "shell.mgz");

  mri_outside = MRISfloodoutside(mri_shell, NULL);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_outside, "out.mgz");
  for (x = 0; x < mri_interior->width; x++)
    for (y = 0; y < mri_interior->height; y++)
      for (z = 0; z < mri_interior->depth; z++) {
        val = (int)MRIgetVoxVal(mri_outside, x, y, z, 0);
        MRIsetVoxVal(mri_interior, x, y, z, 0, val == 1 ? 0 : 1);
      }

  mris->useRealRAS = saved_use_Real_RAS;

  if (interior_alloced) {
    MRI *mri_tmp;
    MATRIX *m, *m_invertz;
    double val;

    // invert one dimension (z) so that ras2vox xform has negative determinant (standard)
    m_invertz = MatrixIdentity(4, NULL);
    *MATRIX_RELT(m_invertz, 3, 3) = -1;
    *MATRIX_RELT(m_invertz, 3, 4) = (mri_interior->depth - 1);
    mri_tmp = MRIclone(mri_interior, NULL);
    m = MatrixMultiply(m_vox2ras, m_invertz, NULL);
    MRIsetVoxelToRasXform(mri_tmp, m);
    MatrixFree(&m);
    MatrixFree(&m_invertz);

    for (x = 0; x < mri_interior->width; x++)
      for (y = 0; y < mri_interior->height; y++)
        for (z = 0; z < mri_interior->depth; z++) {
          val = MRIgetVoxVal(mri_interior, x, y, z, 0);
          MRIsetVoxVal(mri_tmp, x, y, (mri_tmp->depth - 1) - z, 0, val);
        }
    MRIfree(&mri_interior);
    mri_interior = mri_tmp;
  }
  MRIfree(&mri_outside);
  MRIfree(&mri_shell);
  MatrixFree(&m_vox2ras);
  return (mri_interior);
}


/*!
  \brief Creates a bounding volume that encompasses a surface.

  \param mris        Input surface
  \param resolution  Resolution of bounding volume
*/
MRI *MRISmakeBoundingVolume(MRIS *mris, double resolution)
{
  // make surface source volume for reference
  MRI *src = MRIallocFromVolGeom(&mris->vg, MRI_UCHAR, 1, 1);

  // get crop coordinates of surface in the "source" voxel space
  double clo, rlo, slo, chi, rhi, shi;
  MRISsurfaceRASToVoxel(mris, src, mris->xlo, mris->ylo, mris->zlo, &clo, &rlo, &slo);
  MRISsurfaceRASToVoxel(mris, src, mris->xhi, mris->yhi, mris->zhi, &chi, &rhi, &shi);

  // the low surf ras coordinate might not map directly to the voxel low coordinate
  int real_clo = floor(std::min(clo, chi));
  int real_rlo = floor(std::min(rlo, rhi));
  int real_slo = floor(std::min(slo, shi));
  int real_chi = ceil(std::max(clo, chi));
  int real_rhi = ceil(std::max(rlo, rhi));
  int real_shi = ceil(std::max(slo, shi));
  int cropped_width  = real_chi - real_clo + 1;
  int cropped_height = real_rhi - real_rlo + 1;
  int cropped_depth  = real_shi - real_slo + 1;

  // make bounding volume with correct shape
  int width  = ceil(src->xsize * cropped_width  / resolution);
  int height = ceil(src->ysize * cropped_height / resolution);
  int depth  = ceil(src->zsize * cropped_depth  / resolution);
  MRI *dst = MRIalloc(width, height, depth, MRI_FLOAT);

  // set voxel size
  dst->x_r = src->x_r;
  dst->x_a = src->x_a;
  dst->x_s = src->x_s;
  dst->y_r = src->y_r;
  dst->y_a = src->y_a;
  dst->y_s = src->y_s;
  dst->z_r = src->z_r;
  dst->z_a = src->z_a;
  dst->z_s = src->z_s;

  // set voxel size
  dst->xsize = resolution;
  dst->ysize = resolution;
  dst->zsize = resolution;

  // get ras at center of original bounding box and set as center ras
  MATRIX *vox2ras = MRIxfmCRS2XYZ(src, 0);
  MATRIX *c_crs = MatrixAlloc(4, 1, MATRIX_REAL);
  c_crs->rptr[1][1] = real_clo + cropped_width  / 2.0;
  c_crs->rptr[2][1] = real_rlo + cropped_height / 2.0;
  c_crs->rptr[3][1] = real_slo + cropped_depth  / 2.0;
  c_crs->rptr[4][1] = 1;
  MATRIX *c_ras = MatrixMultiply(vox2ras, c_crs, nullptr);
  dst->c_r = c_ras->rptr[1][1];
  dst->c_a = c_ras->rptr[2][1];
  dst->c_s = c_ras->rptr[3][1];

  MatrixFree(&vox2ras);
  MatrixFree(&c_crs);
  MatrixFree(&c_ras);
  MRIfree(&src);

  return dst;
}


/*!
\fn MRI *MRISfillInterior(MRI_SURFACE *mris, double resolution, MRI *mri_dst)
\brief Fills in the interior of a surface by creating a "watertight"
shell and filling everything outside of the shell. This is much faster
but slightly less accurate than a ray-tracing algorithm.  See also
MRISfillInteriorOld() and MRISfillInteriorRibbonTest().
\param mris - input surface
\param resolution - only used if mri_dst is NULL
\param mri_dst - output
*/
MRI *MRISfillInterior(MRI_SURFACE *mris, double resolution, MRI *mri_dst)
{
  int col, row, slc, fno, numu, numv, u, v, nhits, width, height, depth;
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, d0, d1, d2, dmax;
  double px0, py0, pz0, px1, py1, pz1, px, py, pz;
  double fcol, frow, fslc, dcol, drow, dslc, val, val2;
  double vx, vy, vz, vlen, ux, uy, uz, cosa;
  VERTEX *v_0, *v_1, *v_2;
  FACE *f;
  MATRIX *crs, *xyz = NULL, *vox2sras = NULL, *m_vox2ras;
  MRI *mri_cosa, *mri_vlen, *mri_shell, *shellbb, *outsidebb;
  MRI_REGION *region;
  Timer start;

  MRIScomputeMetricProperties(mris);

  if (!mri_dst) {
    // not sure this will work
    // ATH: it doesn't when the surface source geometry differs
    // from resolution or when geometry is not LIA. In the future,
    // this whole section should be replaced by MRISmakeBoundingVolume(),
    // which just needs to be tested more first
    width = ceil((mris->xhi - mris->xlo) / resolution);
    height = ceil((mris->yhi - mris->ylo) / resolution);
    depth = ceil((mris->zhi - mris->zlo) / resolution);
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIsetResolution(mri_dst, resolution, resolution, resolution);
    m_vox2ras = MatrixIdentity(4, NULL);
    *MATRIX_RELT(m_vox2ras, 1, 1) = resolution;
    *MATRIX_RELT(m_vox2ras, 2, 2) = resolution;
    *MATRIX_RELT(m_vox2ras, 3, 3) = resolution;
    *MATRIX_RELT(m_vox2ras, 1, 4) = mris->xlo + mris->vg.c_r;
    *MATRIX_RELT(m_vox2ras, 2, 4) = mris->ylo + mris->vg.c_a;
    *MATRIX_RELT(m_vox2ras, 3, 4) = mris->zlo + mris->vg.c_s;
    MRIsetVoxelToRasXform(mri_dst, m_vox2ras);
  }
  MRIclear(mri_dst);

  dcol = mri_dst->xsize;
  drow = mri_dst->ysize;
  dslc = mri_dst->zsize;

  crs = MatrixAlloc(4, 1, MATRIX_REAL);
  crs->rptr[4][1] = 1;

  mri_cosa = MRIconst(mri_dst->width, mri_dst->height, mri_dst->depth, 1, 0, NULL);
  if (mri_cosa == NULL) {
    printf("ERROR: alloc fill interior cosa\n");
    return (NULL);
  }
  MRIcopyHeader(mri_dst, mri_cosa);

  mri_vlen = MRIconst(mri_dst->width, mri_dst->height, mri_dst->depth, 1, 10.0, NULL);
  if (mri_vlen == NULL) {
    printf("ERROR: alloc fill interior vlen\n");
    return (NULL);
  }
  MRIcopyHeader(mri_dst, mri_vlen);

  // printf("  creating shell\n");
  mri_shell = MRIalloc(mri_dst->width, mri_dst->height, mri_dst->depth, MRI_INT);
  if (mri_shell == NULL) {
    printf("ERROR: alloc fill interior shell\n");
    return (NULL);
  }
  MRIcopyHeader(mri_dst, mri_shell);

  MRIS_SurfRAS2VoxelMap* map = MRIS_makeRAS2VoxelMap(mri_dst, mris);
  MRIS_loadRAS2VoxelMap(map, mri_dst, mris);
    
  /* Create a "watertight" shell of the surface by filling each face in MRI volume */
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];

    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);

    v_0 = &mris->vertices[f->v[0]];
    v_1 = &mris->vertices[f->v[1]];
    v_2 = &mris->vertices[f->v[2]];
    if (f->v[0] == Gdiag_no || f->v[1] == Gdiag_no || f->v[2] == Gdiag_no) DiagBreak();
    /* (x,y,z) for each vertex for face */
    x0 = v_0->x;
    y0 = v_0->y;
    z0 = v_0->z;
    x1 = v_1->x;
    y1 = v_1->y;
    z1 = v_1->z;
    x2 = v_2->x;
    y2 = v_2->y;
    z2 = v_2->z;

    /* Calculate triangle side lengths in voxels */
    d0 = sqrt(SQR((x1 - x0) / dcol) + SQR((y1 - y0) / drow) + SQR((z1 - z0) / dslc));
    d1 = sqrt(SQR((x2 - x1) / dcol) + SQR((y2 - y1) / drow) + SQR((z2 - z1) / dslc));
    d2 = sqrt(SQR((x0 - x2) / dcol) + SQR((y0 - y2) / drow) + SQR((z0 - z2) / dslc));

    /* Divide space between sides into numv parallel lines */
    dmax = (d0 >= d1 && d0 >= d2) ? d0 : (d1 >= d0 && d1 >= d2) ? d1 : d2;
    numv = ceil(2 * dmax);
    numu = ceil(2 * d0);
    /* Fill each line in MRI volume */
    for (v = 0; v <= numv; v++) {
      px0 = x0 + (x2 - x0) * (double)v / (double)numv;
      py0 = y0 + (y2 - y0) * (double)v / (double)numv;
      pz0 = z0 + (z2 - z0) * (double)v / (double)numv;
      px1 = x1 + (x2 - x1) * (double)v / (double)numv;
      py1 = y1 + (y2 - y1) * (double)v / (double)numv;
      pz1 = z1 + (z2 - z1) * (double)v / (double)numv;
      /* Fill each voxel on line in MRI volume */
      for (u = 0; u <= numu; u++) {
        px = px0 + (px1 - px0) * (double)u / (double)numu;
        py = py0 + (py1 - py0) * (double)u / (double)numu;
        pz = pz0 + (pz1 - pz0) * (double)u / (double)numu;
        
        MRIS_useRAS2VoxelMap(map, mri_dst, px, py, pz, &fcol, &frow, &fslc);
        if (vox2sras == NULL) vox2sras = MatrixInverse(map->sras2vox, NULL);
 
        col = nint(fcol);
        row = nint(frow);
        slc = nint(fslc);
        if (col >= 0 && col < mri_dst->width && row >= 0 && row < mri_dst->height && slc >= 0 && slc < mri_dst->depth) {
          MRIsetVoxVal(mri_shell, col, row, slc, 0, 255);
          /* Compute the cosine of the angle to determine whether the
          center of the voxel falls outside of the surface by
          computing the cosine of the angle of the normal and the
          vector that connects the center of the voxel to the closest
          point on the surface. This method is not perfect. It can
          fail when there is a sharp turn in the voxel.  These voxels
          must be part of the shell or else it will not be
          watertight. The "old" method was not able to remove any of
          these voxels and so the interior was biased to be too
          large. */
          crs->rptr[1][1] = col;
          crs->rptr[2][1] = row;
          crs->rptr[3][1] = slc;
          xyz = MatrixMultiply(vox2sras, crs, xyz);  // xyz in surface space
          vx = xyz->rptr[1][1];
          vy = xyz->rptr[2][1];
          vz = xyz->rptr[3][1];
          vlen = sqrt(SQR(px - vx) + SQR(py - vy) + SQR(pz - vz));
          if (vlen > MRIgetVoxVal(mri_vlen, col, row, slc, 0)) continue;  // closest?
          ux = (px - vx) / vlen;
          uy = (py - vy) / vlen;
          uz = (pz - vz) / vlen;
          cosa = ux * fNorm->nx + uy * fNorm->ny + uz * fNorm->nz;
          MRIsetVoxVal(mri_cosa, col, row, slc, 0, cosa);
          MRIsetVoxVal(mri_vlen, col, row, slc, 0, vlen);
        }
      }
    }
  }
  MRIS_freeRAS2VoxelMap(&map);
  MatrixFree(&crs);
  MatrixFree(&xyz);
  if (vox2sras != NULL) MatrixFree(&vox2sras);
  MRIfree(&mri_vlen);

  // Reduce the volume size to speed things up
  region = REGIONgetBoundingBox(mri_shell, 2);
  if (Gdiag_no > 0) {
    printf("  Bounding box around shell: ");
    REGIONprint(stdout, region);
  }

  // Copy shell into bounding box
  shellbb = MRIalloc(region->dx, region->dy, region->dz, MRI_UCHAR);
  if (shellbb == NULL) {
    printf("ERROR: alloc fill shellbb\n");
    return (NULL);
  }
  for (col = region->x; col < region->x + region->dx; col++) {
    for (row = region->y; row < region->y + region->dy; row++) {
      for (slc = region->z; slc < region->z + region->dz; slc++) {
        val = MRIgetVoxVal(mri_shell, col, row, slc, 0);
        MRIsetVoxVal(shellbb, col - region->x, row - region->y, slc - region->z, 0, val);
      }
    }
  }
  MRIfree(&mri_shell);

  // Flood the outside (outside includes shell). This is the part that takes the longest.
  // printf("  flooding outside  ");fflush(stdout);
  outsidebb = MRISfloodoutside(shellbb, NULL);

  /* Note: it is possible that there are voxels outside of the shell
     that do not get flooded because they form a hole and the flood
     cannot reach it. This can particularly happen when the resolution
     is 1mm. It is possible to fix this by inverting the shell and
     finding clusters. The two big clusters are for inside and outside
     the shell and the smaller ones will be holes that are actually
     outside the shell. Of course, once this is done, the flooding
     is not needed. The clustering takes much longer than the flooding.
     The clustering could also be done in a mask made by dilating
     the shell a few times. I have not implemented this fix because
     I don't think it happens that often, if ever, for resolutions
     of 0.5 or finer.
   */

  // Add the shell to the interior, remove voxels that are mostly on
  // the outside of the surface
  nhits = 0;
  for (col = region->x; col < region->x + region->dx; col++) {
    for (row = region->y; row < region->y + region->dy; row++) {
      for (slc = region->z; slc < region->z + region->dz; slc++) {
        val = MRIgetVoxVal(shellbb, col - region->x, row - region->y, slc - region->z, 0);
        val2 = MRIgetVoxVal(outsidebb, col - region->x, row - region->y, slc - region->z, 0);
        if (val < 0.5 && val2 > 0.5) continue;
        // Don't include if voxel is on the outside. Not perfect
        cosa = MRIgetVoxVal(mri_cosa, col, row, slc, 0);
        if (cosa >= 0) {
          MRIsetVoxVal(mri_dst, col, row, slc, 0, 1);
          nhits++;
        }
      }
    }
  }
  MRIfree(&mri_cosa);
  MRIfree(&shellbb);
  MRIfree(&outsidebb);
  free(region);
  // printf("  Found %d voxels in interior, volume = %g\n",nhits,
  // nhits*mri_dst->xsize*mri_dst->ysize*mri_dst->zsize);
  if (Gdiag_no > 0) printf("  MRISfillInterior t = %g\n", start.seconds());

  return (mri_dst);
}

/*!
\fn int MRISfillInteriorRibbonTest(char *subject, int UseNew, FILE *fp)
\brief Runs a test on MRISfillInterior() by comparing its results to
the ribbon.mgz file.  The ribbon.mgz file is the gold standard
generated using a ray tracing algorithm. Typical results are that the
new MRISfillInterior() will overlap ribbon.mgz to better than 99.5%
but is on the order of 20 times faster.
*/
int MRISfillInteriorRibbonTest(char *subject, int UseNew, FILE *fp)
{
  FSENV *fsenv;
  char tmpstr[4000];
  const char *hemistr = NULL, *surfname = NULL;
  MRI *ribbon;
  int hemi, surftype, c, r, s, nfp, nfn, ntp, v, vrib, wmval = 0, ctxval = 0;
  MRIS *surf;
  MRI *mri;

  if (UseNew)
    setenv("USE_NEW_FILL_INTERIOR", "1", 1);
  else
    setenv("USE_NEW_FILL_INTERIOR", "0", 1);

  fsenv = FSENVgetenv();

  // ribbon.mgz is 3 or 42 if within the ribbon or 2 or 41 if within the wm surface
  sprintf(tmpstr, "%s/%s/mri/ribbon.mgz", fsenv->SUBJECTS_DIR, subject);
  ribbon = MRIread(tmpstr);
  if (ribbon == NULL) {
    printf("ERROR: cannot find %s\n", tmpstr);
    return (1);
  }

  mri = MRIcopy(ribbon, NULL);

  fprintf(fp, "%s   ", subject);
  for (surftype = 0; surftype < 2; surftype++) {
    if (surftype == 0) surfname = "white";
    if (surftype == 1) surfname = "pial";

    for (hemi = 0; hemi < 2; hemi++) {
      if (hemi == 0) {
        hemistr = "lh";
        wmval = 2;
        ctxval = 3;
      }
      if (hemi == 1) {
        hemistr = "rh";
        wmval = 41;
        ctxval = 42;
      }

      // printf("%s %s %s\n",subject,hemistr,surfname); fflush(stdout);
      sprintf(tmpstr, "%s/%s/surf/%s.%s", fsenv->SUBJECTS_DIR, subject, hemistr, surfname);
      surf = MRISread(tmpstr);
      if (surf == NULL) {
        printf("ERROR: cannot load %s\n", tmpstr);
        return (1);
      }
      MRIclear(mri);
      if (UseNew) MRISfillInterior(surf, 1, mri);      // resolution = 1
      if (!UseNew) MRISfillInteriorOld(surf, 1, mri);  // resolution = 1

      nfp = 0;  // false positive - not in ribbon but in interior
      nfn = 0;  // false negative - in ribbon but not in interior
      ntp = 0;  // total number within the surface in ribbon.mgz
      for (c = 0; c < mri->width; c++) {
        for (r = 0; r < mri->height; r++) {
          for (s = 0; s < mri->depth; s++) {
            vrib = MRIgetVoxVal(ribbon, c, r, s, 0);
            v = MRIgetVoxVal(mri, c, r, s, 0);
            if (surftype == 0) {  // white
              // only look for ribbon voxels that are wmval (2 or 41)
              if (vrib == wmval) ntp++;
              if (vrib == wmval && v == 1) continue;
              if (vrib != wmval && v == 1) nfp++;
              if (vrib == wmval && v == 0) nfn++;
            }
            if (surftype == 1) {  // pial
              // look for ribbon voxels that are either wmval (2 or 41)
              // or cortex value (3 or 42)
              if (vrib == wmval || vrib == ctxval) ntp++;
              if ((vrib == wmval || vrib == ctxval) && v == 1) continue;
              if ((vrib != wmval && vrib != ctxval) && v == 1) nfp++;
              if ((vrib == wmval || vrib == ctxval) && v == 0) nfn++;
            }
          }
        }
      }
      printf("#P# %s %s %s %3d %3d %6d\n", subject, hemistr, surfname, nfp, nfn, ntp);
      fflush(stdout);
      fprintf(fp, "%4d %4d %6d   ", nfp, nfn, ntp);

    }  // surf type
  }    // hemi
  fprintf(fp, "\n");
  fflush(fp);

  MRIfree(&mri);
  return (0);
}
