/**
 * @brief utilities for computing connected components. See also volcluster.cpp
 *
 */
/*
 * Original Author: Florent Segonne
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

/* connectcomp.c */

#include "connectcomp.h"
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mri.h"

static int xoff6[6] = {1, 0, 0, -1, 0, 0};
static int yoff6[6] = {0, 1, 0, 0, -1, 0};
static int zoff6[6] = {0, 0, 1, 0, 0, -1};

static int xoff26[26] = {0, 1, 0, 0, -1, 0, 1, 1, -1, -1, 1, 1, -1, -1, 0, 0, 0, 0, -1, -1, 1, 1, -1, -1, 1, 1};
static int yoff26[26] = {1, 0, 0, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1};
static int zoff26[26] = {0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1};

void RemoveHoles(MRI *orivol)
{
  /* This function assumes the object is disconnected to the volume boundary.
     It first finds the bkground CC that connected with the volume boundary,
     then set all the voxels of the volume to object value(1) except for this CC.
     See also MRIremoveSliceHoles(),  MRIremoveVolumeHoles(), MRIremoveVolumeIslands()
     and other functions in volcluster.cpp.
   */

  MRI *tmpvol;
  MRI *Label;
  int i, j, k, curSize;
  POINTI seed;
  int minX, minY, minZ, maxX, maxY, maxZ;
  int XN, YN, ZN;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  Label = MRIalloc(orivol->width, orivol->height, orivol->depth, MRI_INT);
  MRIcopyHeader(orivol, Label);

  tmpvol = MRIalloc(orivol->width, orivol->height, orivol->depth, MRI_UCHAR);
  MRIcopyHeader(orivol, tmpvol);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        MRIsetVoxVal(Label, j, i, k,0,0); /* Initialization */

        /* Invert the volume inorder to do Connected-Component labelling on
           background */
        if (MRIgetVoxVal(orivol, j, i, k, 0) <= 0)
          MRIvox(tmpvol, j, i, k) = 1;
        else
          MRIvox(tmpvol, j, i, k) = 0;
      }

  /* Find a seed for the boundary CC. Here we use the boundary of X-axis */
  for (j = 0; j < XN; j++) {
    if (MRIvox(tmpvol, j, 0, 0) != 0 && MRIgetVoxVal(Label, j, 0, 0,0) == 0) {
      seed.x = j;
      seed.y = 0;
      seed.z = 0;
      GrassFire6(tmpvol, Label, 1, &seed, &curSize, &minX, &maxX, &minY, &maxY, &minZ, &maxZ);
    }
  }

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        if (MRIgetVoxVal(Label, j, i, k,0) == 0) MRIsetVoxVal(orivol, j, i, k, 0, 1);
      }

  MRIfree(&Label);
  MRIfree(&tmpvol);

  return;
}

void GrassFire(MRI *orivol,
               MRI *Label,
               int label,
               POINTI *Pt,
               int *curSize,
               int *minX,
               int *maxX,
               int *minY,
               int *maxY,
               int *minZ,
               int *maxZ)
{
  /* This function does binary region growing from seed Pt.
     It assumes that object has value != 0, while bkground = 0.
     minX,maxX,...,maxZ denote the boundary of the current region.
     curSize denotes the size (#ofpoints) of the current region.
     label is the LABEL assigned to all points of current region.
     If a point has LABEL = 0, then it's a unlabelled point or bkground point.
   */

  POINTI cPt, nPt;
  MYqueue NeiQ;
  int ci, cj, ck, ni, nj, nk;
  int ioff, joff, koff;
  int XN, YN, ZN;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  (*minX) = Pt->x;
  (*maxX) = Pt->x;
  (*minY) = Pt->y;
  (*maxY) = Pt->y;
  (*minZ) = Pt->z;
  (*maxZ) = Pt->z;

  NeiQ = myQueue(sizeof(POINTI));
  MRIsetVoxVal(Label, Pt->x, Pt->y, Pt->z,0, label);

  myQueuePush(NeiQ, Pt);

  (*curSize) = 0;

  while (!myQueueIsEmpty(NeiQ)) {
    myQueuePop(NeiQ, &cPt);
    (*curSize)++;
    ci = cPt.y;
    cj = cPt.x;
    ck = cPt.z;

    if ((*minX) > cj) (*minX) = cj;
    if ((*maxX) < cj) (*maxX) = cj;
    if ((*minY) > ci) (*minY) = ci;
    if ((*maxY) < ci) (*maxY) = ci;
    if ((*minZ) > ck) (*minZ) = ck;
    if ((*maxZ) < ck) (*maxZ) = ck;

    for (ioff = -1; ioff <= 1; ioff++)
      for (joff = -1; joff <= 1; joff++)
        for (koff = -1; koff <= 1; koff++) { /* 26 connected neighbourhood */
          ni = ci + ioff;
          nj = cj + joff;
          nk = ck + koff;
          if (ni >= 0 && ni < YN && nj >= 0 && nj < XN && nk >= 0 && nk < ZN) {
            if (MRIgetVoxVal(Label, nj, ni, nk,0) == 0 && MRIgetVoxVal(orivol, nj, ni, nk, 0) > 0) {
              /* Unlabelled object point found */
              nPt.x = nj;
              nPt.y = ni;
              nPt.z = nk;
              MRIsetVoxVal(Label, nj, ni, nk,0, label);
              myQueuePush(NeiQ, &nPt);
            }
          }
        }
  }

  myQueueDelete(NeiQ);
  return;
}

void GrassFire6(MRI *orivol,
                MRI *Label,
                int label,
                POINTI *Pt,
                int *curSize,
                int *minX,
                int *maxX,
                int *minY,
                int *maxY,
                int *minZ,
                int *maxZ)
{
  /* This function does binary region growing from seed Pt.
     It assumes that object has value != 0, while bkground = 0.
     minX,maxX,...,maxZ denote the boundary of the current region.
     curSize denotes the size (#ofpoints) of the current region.
     label is the LABEL assigned to all points of current region.
     If a point has LABEL = 0, then it's a unlabelled point or bkground point.
   */

  POINTI cPt, nPt;
  MYqueue NeiQ;
  int ci, cj, ck, ni, nj, nk;
  int XN, YN, ZN;
  int index;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  (*minX) = Pt->x;
  (*maxX) = Pt->x;
  (*minY) = Pt->y;
  (*maxY) = Pt->y;
  (*minZ) = Pt->z;
  (*maxZ) = Pt->z;

  NeiQ = myQueue(sizeof(POINTI));
  MRIsetVoxVal(Label, Pt->x, Pt->y, Pt->z,0, label);

  myQueuePush(NeiQ, Pt);

  (*curSize) = 0;

  while (!myQueueIsEmpty(NeiQ)) {
    myQueuePop(NeiQ, &cPt);
    (*curSize)++;
    ci = cPt.y;
    cj = cPt.x;
    ck = cPt.z;

    if ((*minX) > cj) (*minX) = cj;
    if ((*maxX) < cj) (*maxX) = cj;
    if ((*minY) > ci) (*minY) = ci;
    if ((*maxY) < ci) (*maxY) = ci;
    if ((*minZ) > ck) (*minZ) = ck;
    if ((*maxZ) < ck) (*maxZ) = ck;

    for (index = 0; index < 6; index++) {
      ni = ci + yoff6[index];
      nj = cj + xoff6[index];
      nk = ck + zoff6[index];

      if (ni >= 0 && ni < YN && nj >= 0 && nj < XN && nk >= 0 && nk < ZN) {
        if (MRIgetVoxVal(Label, nj, ni, nk,0) == 0 && MRIgetVoxVal(orivol, nj, ni, nk, 0) > 0) {
          /* Unlabelled object point found */
          nPt.x = nj;
          nPt.y = ni;
          nPt.z = nk;
          MRIsetVoxVal(Label, nj, ni, nk,0, label);
          myQueuePush(NeiQ, &nPt);
        }
      }
    }
  }

  myQueueDelete(NeiQ);
  return;
}

void GrassFire18(MRI *orivol,
                 MRI *Label,
                 int label,
                 POINTI *Pt,
                 int *curSize,
                 int *minX,
                 int *maxX,
                 int *minY,
                 int *maxY,
                 int *minZ,
                 int *maxZ)
{
  /* This function does binary region growing from seed Pt.
     It assumes that object has value != 0, while bkground = 0.
     minX,maxX,...,maxZ denote the boundary of the current region.
     curSize denotes the size (#ofpoints) of the current region.
     label is the LABEL assigned to all points of current region.
     If a point has LABEL = 0, then it's a unlabelled point or bkground point.
   */

  POINTI cPt, nPt;
  MYqueue NeiQ;
  int ci, cj, ck, ni, nj, nk;
  int XN, YN, ZN;
  int index;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  (*minX) = Pt->x;
  (*maxX) = Pt->x;
  (*minY) = Pt->y;
  (*maxY) = Pt->y;
  (*minZ) = Pt->z;
  (*maxZ) = Pt->z;

  NeiQ = myQueue(sizeof(POINTI));
  MRIsetVoxVal(Label, Pt->x, Pt->y, Pt->z, 0, label);

  myQueuePush(NeiQ, Pt);

  (*curSize) = 0;

  while (!myQueueIsEmpty(NeiQ)) {
    myQueuePop(NeiQ, &cPt);
    (*curSize)++;
    ci = cPt.y;
    cj = cPt.x;
    ck = cPt.z;

    if ((*minX) > cj) (*minX) = cj;
    if ((*maxX) < cj) (*maxX) = cj;
    if ((*minY) > ci) (*minY) = ci;
    if ((*maxY) < ci) (*maxY) = ci;
    if ((*minZ) > ck) (*minZ) = ck;
    if ((*maxZ) < ck) (*maxZ) = ck;

    for (index = 0; index < 18; index++) {
      ni = ci + yoff26[index];
      nj = cj + xoff26[index];
      nk = ck + zoff26[index];

      if (ni >= 0 && ni < YN && nj >= 0 && nj < XN && nk >= 0 && nk < ZN) {
        if (MRIgetVoxVal(Label, nj, ni, nk,0) == 0 && MRIgetVoxVal(orivol, nj, ni, nk, 0) > 0) {
          /* Unlabelled object point found */
          nPt.x = nj;
          nPt.y = ni;
          nPt.z = nk;
          MRIsetVoxVal(Label, nj, ni, nk, 0, label);
          myQueuePush(NeiQ, &nPt);
        }
      }
    }
  }

  myQueueDelete(NeiQ);
  return;
}

void GetLargestCC6(MRI *orivol)
{
  /* This function keeps the largest CC, and reset all other CC to bgvalue (0) */
  MRI *Label;
  int i, j, k;
  int maxSize, maxLabel, curSize, curLabel;
  int minX, minY, minZ, maxX, maxY, maxZ;
  int XN, YN, ZN;
  POINTI Pt;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  Label = MRIalloc(XN, YN, ZN, MRI_INT);
  MRIcopyHeader(orivol, Label);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) MRIsetVoxVal(Label, j, i, k,0, 0);

  curLabel = 1;
  maxSize = 0;
  maxLabel = 1;

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        if (MRIgetVoxVal(orivol, j, i, k, 0) > 0 && MRIgetVoxVal(Label, j, i, k,0) == 0) {
          Pt.x = j;
          Pt.y = i;
          Pt.z = k;
          GrassFire6(orivol, Label, curLabel, &Pt, &curSize, &minX, &maxX, &minY, &maxY, &minZ, &maxZ);
          if (maxSize < curSize) {
            maxSize = curSize;
            maxLabel = curLabel;
          }
          curLabel++;
        }
      }

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        if (MRIgetVoxVal(Label, j, i, k,0) != maxLabel) MRIsetVoxVal(orivol, j, i, k, 0, 0);
      }

  MRIfree(&Label);
  return;
}

void GetLargestCC18(MRI *orivol)
{
  /* This function keeps the largest CC, and reset all other CC to bgvalue (0) */
  MRI *Label;
  int i, j, k;
  int maxSize, maxLabel, curSize, curLabel;
  int minX, minY, minZ, maxX, maxY, maxZ;
  int XN, YN, ZN;
  POINTI Pt;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  Label = MRIalloc(XN, YN, ZN, MRI_INT);
  MRIcopyHeader(orivol, Label);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) MRIsetVoxVal(Label, j, i, k,0,0);

  curLabel = 1;
  maxSize = 0;
  maxLabel = 1;

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        if (MRIgetVoxVal(orivol, j, i, k, 0) > 0 && MRIgetVoxVal(Label, j, i, k,0) == 0) {
          Pt.x = j;
          Pt.y = i;
          Pt.z = k;
          GrassFire18(orivol, Label, curLabel, &Pt, &curSize, &minX, &maxX, &minY, &maxY, &minZ, &maxZ);
          if (maxSize < curSize) {
            maxSize = curSize;
            maxLabel = curLabel;
          }
          curLabel++;
        }
      }

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        if (MRIgetVoxVal(Label, j, i, k,0) != maxLabel) MRIsetVoxVal(orivol, j, i, k, 0, 0);
      }

  MRIfree(&Label);
  return;
}

MRI *Dilation6(MRI *ori, MRI *out, int R)
{
  int i, j, k, index, ci, cj, ck, count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  MRIcopyHeader(ori, tmpvol);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori, j, i, k, 0);
      }

  for (count = 1; count <= R; count++) {
    for (i = 0; i < YN; i++)
      for (j = 0; j < XN; j++)
        for (k = 0; k < ZN; k++) {
          if (MRIvox(tmpvol, j, i, k) == 1) {
            MRIsetVoxVal(out, j, i, k, 0, 1);
            continue;
          }

          MRIsetVoxVal(out, j, i, k, 0, 0);

          for (index = 0; index < 6; index++) {
            ci = i + yoff6[index];
            cj = j + xoff6[index];
            ck = k + zoff6[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN) continue;
            if (MRIvox(tmpvol, cj, ci, ck) == 1) {
              MRIsetVoxVal(out, j, i, k, 0, 1);
              break;
            }
          }
        }

    if (count < R) {
      for (i = 0; i < YN; i++)
        for (j = 0; j < XN; j++)
          for (k = 0; k < ZN; k++) {
            MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(out, j, i, k, 0);
          }
    }
  }

  MRIfree(&tmpvol);
  return (out);
}

MRI *Erosion6(MRI *ori, MRI *out, int R)
{
  int i, j, k, index, ci, cj, ck, count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  MRIcopyHeader(ori, tmpvol);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori, j, i, k, 0);
      }

  for (count = 1; count <= R; count++) {
    for (i = 0; i < YN; i++)
      for (j = 0; j < XN; j++)
        for (k = 0; k < ZN; k++) {
          if (MRIvox(tmpvol, j, i, k) == 0) {
            MRIsetVoxVal(out, j, i, k, 0, 0);
            continue;
          }

          MRIsetVoxVal(out, j, i, k, 0, 1);

          for (index = 0; index < 6; index++) {
            ci = i + yoff6[index];
            cj = j + xoff6[index];
            ck = k + zoff6[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN) continue;

            if (MRIvox(tmpvol, cj, ci, ck) == 0) {
              MRIsetVoxVal(out, j, i, k, 0, 0);
              break;
            }
          }
        }

    if (count < R) {
      for (i = 0; i < YN; i++)
        for (j = 0; j < XN; j++)
          for (k = 0; k < ZN; k++) {
            MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(out, j, i, k, 0);
          }
    }
  }

  MRIfree(&tmpvol);
  return (out);
}

MRI *Dilation26(MRI *ori, MRI *out, int R)
{
  int i, j, k, index, ci, cj, ck, count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  if (tmpvol == NULL) {
    printf("Unable to allocate memory. Exit.\n");
    exit(0);
  }
  MRIcopyHeader(ori, tmpvol);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori, j, i, k, 0);
      }

  for (count = 1; count <= R; count++) {
    for (i = 0; i < YN; i++)
      for (j = 0; j < XN; j++)
        for (k = 0; k < ZN; k++) {
          if (MRIvox(tmpvol, j, i, k) == 1) {
            MRIsetVoxVal(out, j, i, k, 0, 1);
            continue;
          }

          MRIsetVoxVal(out, j, i, k, 0, 0);

          for (index = 0; index < 26; index++) {
            ci = i + yoff26[index];
            cj = j + xoff26[index];
            ck = k + zoff26[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN) continue;
            if (MRIvox(tmpvol, cj, ci, ck) == 1) {
              MRIsetVoxVal(out, j, i, k, 0, 1);
              break;
            }
          }
        }

    if (count < R) {
      for (i = 0; i < YN; i++)
        for (j = 0; j < XN; j++)
          for (k = 0; k < ZN; k++) {
            MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(out, j, i, k, 0);
          }
    }
  }

  MRIfree(&tmpvol);
  return (out);
}

MRI *Erosion26(MRI *ori, MRI *out, int R)
{
  int i, j, k, index, ci, cj, ck, count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  MRIcopyHeader(ori, tmpvol);

  for (i = 0; i < YN; i++)
    for (j = 0; j < XN; j++)
      for (k = 0; k < ZN; k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori, j, i, k, 0);
      }

  for (count = 1; count <= R; count++) {
    for (i = 0; i < YN; i++)
      for (j = 0; j < XN; j++)
        for (k = 0; k < ZN; k++) {
          if (MRIvox(tmpvol, j, i, k) == 0) {
            MRIsetVoxVal(out, j, i, k, 0, 0);
            continue;
          }

          MRIsetVoxVal(out, j, i, k, 0, 1);

          for (index = 0; index < 26; index++) {
            ci = i + yoff26[index];
            cj = j + xoff26[index];
            ck = k + zoff26[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN) continue;

            if (MRIvox(tmpvol, cj, ci, ck) == 0) {
              MRIsetVoxVal(out, j, i, k, 0, 0);
              break;
            }
          }
        }

    if (count < R) {
      for (i = 0; i < YN; i++)
        for (j = 0; j < XN; j++)
          for (k = 0; k < ZN; k++) {
            MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(out, j, i, k, 0);
          }
    }
  }

  MRIfree(&tmpvol);
  return (out);
}

MRI *BinaryOpen6(MRI *ori, MRI *out, int R)
{
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;

  tmpvol = Erosion6(ori, tmpvol, R);
  out = Dilation6(tmpvol, out, R);

  MRIfree(&tmpvol);

  return (out);
}

MRI *BinaryOpen26(MRI *ori, MRI *out, int R)
{
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;

  printf("Erosion...\n");
  tmpvol = Erosion26(ori, tmpvol, R);
  printf("Dilation...\n");
  out = Dilation26(tmpvol, out, R);

  MRIfree(&tmpvol);

  return (out);
}

MRI *BinaryClose6(MRI *ori, MRI *out, int R)
{
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;

  tmpvol = Dilation6(ori, tmpvol, R);
  out = Erosion6(tmpvol, out, R);

  MRIfree(&tmpvol);

  return (out);
}

MRI *BinaryClose26(MRI *ori, MRI *out, int R)
{
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;

  tmpvol = Dilation26(ori, tmpvol, R);
  out = Erosion26(tmpvol, out, R);

  MRIfree(&tmpvol);

  return (out);
}

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

/*----------------------------------------------------------------------------
//
//      File: myutil.c
//      A utility library
//
//--------------------------------------------------------------------------*/

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*---------------------------------------------------------------------------
// Construct an empty MYlist with specified capacity and capacity increment
//-------------------------------------------------------------------------*/
MYlist myList2(int elementSize, int capacity, int capacityIncrement)
{
  MYlist list;
  void *data;

  if (elementSize < 1) {
    fprintf(stderr, "myList(): elementSize must be a postive integer!\n");
    exit(1);
  }
  if (capacity < 0) {
    fprintf(stderr, "myList(): capacity must not be negative!\n");
    exit(1);
  }

  list = (MYlist)myMalloc(sizeof(MYlistStruct));

  myListSetElementSize(list, elementSize);
  myListSetSize(list, 0);
  myListSetCapacity(list, capacity);
  myListSetCapacityIncrement(list, capacityIncrement);
  if (capacity == 0)
    myListSetData(list, NULL);
  else {
    data = (void *)myMalloc(elementSize * capacity);
    myListSetData(list, data);
  }

  return (list);
}

/*---------------------------------------------------------------------------
// Construct an empty MYlist with specified capacity and default
// capacity increment as 100
//-------------------------------------------------------------------------*/
MYlist myList1(int elementSize, int capacity)
{
  MYlist list;

  list = myList2(elementSize, capacity, 100);

  return (list);
}

/*---------------------------------------------------------------------------
// Construct an empty MYlist with default capacity as 0 and
// capacity increment as 100
//-------------------------------------------------------------------------*/
MYlist myList(int elementSize)
{
  MYlist list;

  list = myList2(elementSize, 0, 100);

  return (list);
}

/*---------------------------------------------------------------------------
// Construct an empty MYlist with specified size, all the elements are set to
// zero
//-------------------------------------------------------------------------*/
MYlist myListOfSize(int size, int elementSize)
{
  MYlist list;
  char *data;
  int i;
  int capacity, capacityIncrement;

  if (size < 0) {
    fprintf(stderr, "myListOfSize(): size must not be negative!\n");
    exit(1);
  }

  capacity = size;
  capacityIncrement = 100;
  list = myList2(elementSize, capacity, capacityIncrement);
  myListSetSize(list, size);
  data = (char *)myListData(list);
  for (i = 0; i < elementSize * size; i++) data[i] = 0;

  return (list);
}

/*---------------------------------------------------------------------------
// Delete this list
//-------------------------------------------------------------------------*/
void myListDelete(MYlist list)
{
  void *data;

  data = myListData(list);
  free(data);
  free(list);
}

/*---------------------------------------------------------------------------
// Add an element to this list
//-------------------------------------------------------------------------*/
void myListAddElement(MYlist list, void *element)
{
  int size, capacity, elementSize, capacityIncrement;
  void *data;

  size = myListSize(list);
  capacity = myListCapacity(list);
  elementSize = myListElementSize(list);
  data = myListData(list);
  if (size >= capacity) {
    capacityIncrement = myListCapacityIncrement(list);
    capacity += capacityIncrement;
    myListSetCapacity(list, capacity);
    if (data == NULL) {
      /* initial list */
      data = (void *)myMalloc(elementSize * capacity);
    }
    else {
      /* allocate a larger list */
      data = (void *)myRealloc(data, elementSize * capacity);
    }
    myListSetData(list, data);
  }

  memcpy((char *)data + size * elementSize, (char *)element, elementSize);
  myListSetSize(list, size + 1);
}

/*---------------------------------------------------------------------------
// Add an integer to this list (must be a list consists of only integers)
//-------------------------------------------------------------------------*/
void myListAddInt(MYlist list, int element)
{
  int size, capacity, elementSize, capacityIncrement;
  int *data;

  size = myListSize(list);
  capacity = myListCapacity(list);
  elementSize = myListElementSize(list);
  data = (int *)myListData(list);

  if (size >= capacity) {
    capacityIncrement = myListCapacityIncrement(list);
    capacity += capacityIncrement;
    myListSetCapacity(list, capacity);
    if (data == NULL) {
      /* initial list */
      data = (int *)myMalloc(elementSize * capacity);
    }
    else {
      /* allocate a larger list */
      data = (int *)myRealloc(data, elementSize * capacity);
    }
    myListSetData(list, data);
  }

  data[size] = element;
  myListSetSize(list, size + 1);
}

/*---------------------------------------------------------------------------
// Add an array to this list
//-------------------------------------------------------------------------*/
void myListAddArray(MYlist list, void *array, int num)
{
  int size, capacity, elementSize, capacityIncrement, actualIncrement;
  void *data;

  size = myListSize(list);
  capacity = myListCapacity(list);
  elementSize = myListElementSize(list);
  data = myListData(list);
  if (size + num > capacity) {
    capacityIncrement = myListCapacityIncrement(list);
    actualIncrement = (capacityIncrement > num) ? capacityIncrement : num;
    capacity += actualIncrement;
    myListSetCapacity(list, capacity);

    if (data == NULL) {
      /* initial list */
      data = (void *)myMalloc(elementSize * capacity);
    }
    else {
      /* allocate a larger list */
      data = (void *)myRealloc(data, elementSize * capacity);
    }
    myListSetData(list, data);
  }

  memcpy((char *)data + size * elementSize, (char *)array, num * elementSize);
  myListSetSize(list, size + num);
}

/*---------------------------------------------------------------------------
// Insert an element into the list at the specified index
//-------------------------------------------------------------------------*/
int myListInsertElementAt(MYlist list, int index, void *element)
{
  int size, elementSize;
  void *data;
  void *tempPtr;
  char *currentPtr, *nextPtr;
  int i;

  size = myListSize(list);
  elementSize = myListElementSize(list);

  if (index < 0 || index > size - 1) {
    return (-1); /* out of bound error */
  }

  tempPtr = (void *)myMalloc(elementSize);
  myListAddElement(list, tempPtr);

  data = myListData(list);

  for (i = size - 1; i >= index; i--) {
    currentPtr = (char *)data + i * elementSize;
    nextPtr = (char *)currentPtr + elementSize;
    memcpy(nextPtr, currentPtr, elementSize);
  }

  memcpy((char *)data + index * elementSize, (char *)element, elementSize);

  return (0);
}

/*---------------------------------------------------------------------------
// Retrieve an element from this list at a given index
//-------------------------------------------------------------------------*/
int myListElementAt(MYlist list, int index, void *element)
{
  int size, elementSize;
  void *data;

  size = myListSize(list);
  elementSize = myListElementSize(list);
  data = myListData(list);

  if (index < 0 || index > size - 1) {
    return (-1); /* out of bound error */
  }
  memcpy((char *)element, (char *)data + index * elementSize, elementSize);

  return (0);
}

/*---------------------------------------------------------------------------
// Sets a list element at a given index
//-------------------------------------------------------------------------*/
int myListSetElementAt(MYlist list, int index, void *element)
{
  int size, elementSize;
  void *data;

  size = myListSize(list);
  elementSize = myListElementSize(list);
  data = myListData(list);

  if (index < 0 || index > size - 1) {
    return (-1); /* out of bound error */
  }

  memcpy((char *)data + index * elementSize, (char *)element, elementSize);

  return (0);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
int myListRemoveElementAt(MYlist list, int index)
{
  int size, elementSize;
  void *data;
  char *currentPtr, *nextPtr;
  int i;

  size = myListSize(list);
  elementSize = myListElementSize(list);
  data = myListData(list);

  if (index < 0 || index > size - 1) {
    return (-1); /* out of bound error */
  }

  for (i = index; i < size - 1; i++) {
    currentPtr = (char *)data + i * elementSize;
    nextPtr = (char *)currentPtr + elementSize;
    memcpy(currentPtr, nextPtr, elementSize);
  }

  myListSetSize(list, size - 1);

  return (0);
}

/*---------------------------------------------------------------------------
// Removes all elements from this list and sets its size to zero
//-------------------------------------------------------------------------*/
void myListRemoveAllElements(MYlist list) { myListSetSize(list, 0); }

/*---------------------------------------------------------------------------
// Trim this list to current size
//-------------------------------------------------------------------------*/
void myListTrim(MYlist list)
{
  void *data;
  int size, elementSize;

  size = myListSize(list);
  elementSize = myListElementSize(list);
  data = myListData(list);

  data = (void *)myRealloc(data, elementSize * size);
  myListSetData(list, data);
  myListSetCapacity(list, size);
}

void myListInfo(MYlist list)
{
  int elementSize, size, capacity, capacityIncrement;

  elementSize = myListElementSize(list);
  size = myListSize(list);
  capacity = myListCapacity(list);
  capacityIncrement = myListCapacityIncrement(list);

  printf("         elementSize = %d\n", elementSize);
  printf("                size = %d\n", size);
  printf("            capacity = %d\n", capacity);
  printf("   capacityIncrement = %d\n", capacityIncrement);
  printf("\n");
}

/*---------------------------------------------------------------------------
// pop out the top element from the stack
//-------------------------------------------------------------------------*/
void myStackPop(MYstack stack, void *element)
{
  int size, elementSize;
  void *data;

  size = myListSize(stack);
  elementSize = myListElementSize(stack);
  data = myListData(stack);

  memcpy((char *)element, (char *)data + (size - 1) * elementSize, elementSize);
  myListSetSize(stack, size - 1);
}

/*---------------------------------------------------------------------------
// Construct an empty MYqueue with specified capacity and capacity increment
//-------------------------------------------------------------------------*/
MYqueue myQueue2(int elementSize, int capacity, int capacityIncrement)
{
  MYqueue queue;

  queue = (MYqueue)myMalloc(sizeof(MYqueueStruct));
  queue->start = 0;
  queue->end = -1;
  queue->list = myList2(elementSize, capacity, capacityIncrement);

  return (queue);
}

/*---------------------------------------------------------------------------
// Construct an empty MYqueue with specified capacity and default
// capacity increment as 100
//-------------------------------------------------------------------------*/
MYqueue myQueue1(int elementSize, int capacity)
{
  MYqueue queue;

  queue = myQueue2(elementSize, capacity, 100);

  return (queue);
}

/*---------------------------------------------------------------------------
// default constructor
//-------------------------------------------------------------------------*/
MYqueue myQueue(int elementSize)
{
  MYqueue queue;

  queue = myQueue2(elementSize, 0, 100);

  return (queue);
}

/* private functinos */

void myQueueEnsureSize(MYqueue queue) { myListSetSize(queue->list, myQueueSize(queue)); }

/*---------------------------------------------------------------------------
// delete the queue and its memory
//-------------------------------------------------------------------------*/
void myQueueDelete(MYqueue queue)
{
  myListDelete(queue->list);
  free(queue);
}

/*---------------------------------------------------------------------------
// removes all elements from the queue without releasing memory
//-------------------------------------------------------------------------*/
void myQueueRemoveAllElements(MYqueue queue)
{
  queue->start = 0;
  queue->end = -1;
  myQueueEnsureSize(queue);
}

/* a private function for clean queue, ie. move things to the front */
void myQueueMoveToFront(MYqueue queue)
{
  void *data;
  void *s1, *s2;
  int elementSize, start, end;

  elementSize = myListElementSize(queue->list);
  data = myListData(queue->list);
  start = queue->start;
  end = queue->end;
  s2 = (char *)data + start * elementSize;
  s1 = data;
  memmove(s1, s2, (end - start + 1) * elementSize);
  queue->end = end - start;
  queue->start = 0;
  myQueueEnsureSize(queue);
}

/*---------------------------------------------------------------------------
// push an element to the end of the queue
//-------------------------------------------------------------------------*/
void myQueuePush(MYqueue queue, void *element)
{
  int size, capacity;
  int q = MY_QUEUE_Q;
  /* this is the factor determines the frequency of cleaning */
  /* a factor of n denotes that (n-1)*memory(queue) will be
  wasted */

  size = myQueueSize(queue);
  capacity = myListCapacity(queue->list);

  if (queue->end >= (capacity - 1) && q * size < capacity) /* move the block to front, and release more memory */
    myQueueMoveToFront(queue);

  /* just keep adding the element */
  queue->end = queue->end + 1;
  myListAddElement(queue->list, element);
}

/*---------------------------------------------------------------------------
// push an array to the end of the queue
//-------------------------------------------------------------------------*/
void myQueuePushArray(MYqueue queue, void *array, int num)
{
  int size, capacity;
  int q = MY_QUEUE_Q;
  /* this is the factor determines the frequency of cleaning */
  /* a factor of n denotes that (n-1)*memory(queue) will be
  wasted */

  size = myQueueSize(queue);
  capacity = myListCapacity(queue->list);

  if (queue->end >= (capacity - 1) && q * size < capacity) /* move the block to front, and release more memory */
    myQueueMoveToFront(queue);

  /* just keep adding the array */
  queue->end = queue->end + num;
  myListAddArray(queue->list, array, num);
}

/*---------------------------------------------------------------------------
// pop out the first element from the queue

//-------------------------------------------------------------------------*/
int myQueuePop(MYqueue queue, void *element)
{
  int start;

  start = queue->start;
  if (myQueueIsEmpty(queue) == 0) {
    myListElementAt(queue->list, start, element); /* get the first one */
    queue->start = start + 1;
    return (1); /* correctly got the result */
  }
  else {
    return (0); /* element is undefined */
  }
}

/*---------------------------------------------------------------------------
// trim the queue to its current size
//-------------------------------------------------------------------------*/
void myQueueTrim(MYqueue queue)
{
  myQueueMoveToFront(queue);

  myListTrim(queue->list);
}

/*---------------------------------------------------------------------------
// extract an array from the queue
//-------------------------------------------------------------------------*/
void *myQueueToArray(MYqueue queue)
{
  void *arr, *data;
  int size, elementSize;

  size = myQueueSize(queue);
  elementSize = myListElementSize(queue->list);

  if (size == 0) {
    arr = NULL;
  }
  else {
    data = myListData(queue->list);
    arr = (void *)myMalloc(elementSize * size);
    memcpy(arr, (char *)data + queue->start * elementSize, size * elementSize);
  }

  return (arr);
}

void myQueueInfo(MYqueue queue)
{
  printf("               start = %d\n", queue->start);
  printf("               end   = %d\n", queue->end);
  myListInfo(queue->list);
}

void *myMalloc(int size)
{
  void *value = (void *)malloc(size);
  if (value == 0) myError("Virtual memory exhausted!");

  return value;
}

void *myRealloc(void *ptr, int size)
{
  void *value = (void *)realloc(ptr, size);
  if (value == 0) myError("Virtual memory exhausted!");

  return value;
}

/* Prints error message */
void myError(const char error_text[])
{
  fprintf(stderr, "Utility run-time error:\n");
  fprintf(stderr, "%s\n", error_text);
  fprintf(stderr, "Now exiting to system.\n");
  exit(1);
}


