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


/* subroutines.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include "mri.h"
#include "subroutines.h"
#include "myutil.h"

static int xoff6[6] = {
                        1,0,0,-1, 0, 0
                      };
static int yoff6[6] = {
                        0,1,0, 0,-1, 0
                      };
static int zoff6[6] = {
                        0,0,1, 0, 0,-1
                      };

static int xoff26[26] = {
                          0,1,0, 0,-1, 0, 1, 1,-1,-1,1, 1,-1,-1,0, 0, 0, 0,-1,-1, 1, 1,-1,-1, 1, 1
                        };
static int yoff26[26] = {
                          1,0,0,-1, 0, 0, 1,-1, 1,-1,0, 0, 0, 0,1, 1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1
                        };
static int zoff26[26] = {
                          0,0,1, 0, 0,-1, 0, 0, 0, 0,1,-1, 1,-1,1,-1, 1,-1,-1, 1,-1, 1,-1, 1,-1, 1
                        };

void RemoveHoles(MRI *orivol) {
  /* This function assumes the object is disconnected to the volume boundary.
     It first finds the bkground CC that connected with the volume boundary,
     then set all the voxels of the volume to object value(1) except for this CC.
   */

  MRI *tmpvol;
  MRI *Label;
  int i,j,k, curSize;
  POINTI seed;
  int minX,minY,minZ,maxX,maxY,maxZ;
  int XN, YN, ZN;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;


  Label = MRIalloc(orivol->width, orivol->height, orivol->depth, MRI_INT);
  MRIcopyHeader(orivol, Label);

  tmpvol = MRIalloc(orivol->width, orivol->height, orivol->depth, MRI_UCHAR);
  MRIcopyHeader(orivol, tmpvol);

  for (i=0; i<YN;i++)
    for (j=0; j< XN;j++)
      for (k=0;k<ZN;k++) {
        MRIIvox(Label,j,i,k) = 0; /* Initialization */

        /* Invert the volume inorder to do Connected-Component labelling on
           background */
        if (MRIgetVoxVal(orivol,j,i,k,0) <= 0) MRIvox(tmpvol,j,i,k) = 1;
        else MRIvox(tmpvol, j, i, k) = 0;
      }

  /* Find a seed for the boundary CC. Here we use the boundary of X-axis */
  for (j=0;j<XN;j++) {
    if (MRIvox(tmpvol, j, 0, 0) != 0 && MRIIvox(Label,j,0,0) == 0) {
      seed.x = j;
      seed.y = 0;
      seed.z = 0;
      GrassFire6(tmpvol,Label,1,&seed,&curSize,&minX,&maxX,&minY,&maxY,&minZ,&maxZ);
    }
  }

  for (i=0; i<YN;i++)
    for (j=0; j< XN;j++)
      for (k=0;k<ZN;k++) {
        if (MRIIvox(Label,j,i,k) == 0) MRIsetVoxVal(orivol,j,i,k,0,1);
      }

  MRIfree(&Label);
  MRIfree(&tmpvol);

  return;
}

void  GrassFire(MRI *orivol, MRI *Label, int label, POINTI *Pt,
                int *curSize, int *minX, int *maxX, int *minY, int *maxY,
                int *minZ, int *maxZ) {
  /* This function does binary region growing from seed Pt.
     It assumes that object has value != 0, while bkground = 0.
     minX,maxX,...,maxZ denote the boundary of the current region.
     curSize denotes the size (#ofpoints) of the current region.
     label is the LABEL assigned to all points of current region.
     If a point has LABEL = 0, then it's a unlabelled point or bkground point.
   */

  POINTI cPt, nPt;
  MYqueue NeiQ;
  int ci,cj,ck,ni,nj,nk;
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
  MRIIvox(Label, Pt->x, Pt->y, Pt->z) = label;

  myQueuePush(NeiQ, Pt);

  (*curSize) = 0;

  while (!myQueueIsEmpty(NeiQ)) {
    myQueuePop(NeiQ, &cPt);
    (*curSize)++;
    ci = cPt.y;
    cj = cPt.x;
    ck = cPt.z;


    if ((*minX)>cj) (*minX) = cj;
    if ((*maxX)<cj) (*maxX) = cj;
    if ((*minY)>ci) (*minY) = ci;
    if ((*maxY)<ci) (*maxY) = ci;
    if ((*minZ)>ck) (*minZ) = ck;
    if ((*maxZ)<ck) (*maxZ) = ck;

    for (ioff =-1;ioff<=1; ioff++)
      for (joff=-1;joff<=1;joff++)
        for (koff=-1;koff<=1;koff++) { /* 26 connected neighbourhood */
          ni = ci + ioff;
          nj = cj + joff;
          nk = ck + koff;
          if (ni>=0 && ni<YN && nj>=0 && nj<XN && nk>=0 && nk<ZN) {
            if (MRIIvox(Label, nj, ni, nk) == 0 && MRIgetVoxVal(orivol,nj,ni,nk,0) > 0) {
              /* Unlabelled object point found */
              nPt.x = nj;
              nPt.y = ni;
              nPt.z = nk;
              MRIIvox(Label, nj, ni, nk) = label;
              myQueuePush(NeiQ, &nPt);
            }
          }
        }
  }

  myQueueDelete(NeiQ);
  return;
}



void  GrassFire6(MRI *orivol, MRI *Label, int label, POINTI *Pt,
                 int *curSize, int *minX, int *maxX, int *minY, int *maxY,
                 int *minZ, int *maxZ) {
  /* This function does binary region growing from seed Pt.
     It assumes that object has value != 0, while bkground = 0.
     minX,maxX,...,maxZ denote the boundary of the current region.
     curSize denotes the size (#ofpoints) of the current region.
     label is the LABEL assigned to all points of current region.
     If a point has LABEL = 0, then it's a unlabelled point or bkground point.
   */

  POINTI cPt, nPt;
  MYqueue NeiQ;
  int ci,cj,ck,ni,nj,nk;
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
  MRIIvox(Label, Pt->x, Pt->y, Pt->z) = label;

  myQueuePush(NeiQ, Pt);

  (*curSize) = 0;

  while (!myQueueIsEmpty(NeiQ)) {
    myQueuePop(NeiQ, &cPt);
    (*curSize)++;
    ci = cPt.y;
    cj = cPt.x;
    ck = cPt.z;


    if ((*minX)>cj) (*minX) = cj;
    if ((*maxX)<cj) (*maxX) = cj;
    if ((*minY)>ci) (*minY) = ci;
    if ((*maxY)<ci) (*maxY) = ci;
    if ((*minZ)>ck) (*minZ) = ck;
    if ((*maxZ)<ck) (*maxZ) = ck;

    for (index = 0; index < 6; index++) {
      ni = ci + yoff6[index];
      nj = cj + xoff6[index];
      nk = ck + zoff6[index];

      if (ni>=0 && ni<YN && nj>=0 && nj<XN && nk>=0 && nk<ZN) {
        if (MRIIvox(Label, nj, ni, nk) == 0 && MRIgetVoxVal(orivol,nj,ni,nk,0) > 0) {
          /* Unlabelled object point found */
          nPt.x = nj;
          nPt.y = ni;
          nPt.z = nk;
          MRIIvox(Label, nj, ni, nk) = label;
          myQueuePush(NeiQ, &nPt);
        }

      }
    }
  }

  myQueueDelete(NeiQ);
  return;

}

void  GrassFire18(MRI *orivol, MRI *Label, int label, POINTI *Pt,
                  int *curSize, int *minX, int *maxX, int *minY, int *maxY,
                  int *minZ, int *maxZ) {
  /* This function does binary region growing from seed Pt.
     It assumes that object has value != 0, while bkground = 0.
     minX,maxX,...,maxZ denote the boundary of the current region.
     curSize denotes the size (#ofpoints) of the current region.
     label is the LABEL assigned to all points of current region.
     If a point has LABEL = 0, then it's a unlabelled point or bkground point.
   */

  POINTI cPt, nPt;
  MYqueue NeiQ;
  int ci,cj,ck,ni,nj,nk;
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
  MRIIvox(Label, Pt->x, Pt->y, Pt->z) = label;

  myQueuePush(NeiQ, Pt);

  (*curSize) = 0;

  while (!myQueueIsEmpty(NeiQ)) {
    myQueuePop(NeiQ, &cPt);
    (*curSize)++;
    ci = cPt.y;
    cj = cPt.x;
    ck = cPt.z;


    if ((*minX)>cj) (*minX) = cj;
    if ((*maxX)<cj) (*maxX) = cj;
    if ((*minY)>ci) (*minY) = ci;
    if ((*maxY)<ci) (*maxY) = ci;
    if ((*minZ)>ck) (*minZ) = ck;
    if ((*maxZ)<ck) (*maxZ) = ck;

    for (index = 0; index < 18; index++) {
      ni = ci + yoff26[index];
      nj = cj + xoff26[index];
      nk = ck + zoff26[index];

      if (ni>=0 && ni<YN && nj>=0 && nj<XN && nk>=0 && nk<ZN) {
        if (MRIIvox(Label, nj, ni, nk) == 0 && MRIgetVoxVal(orivol,nj,ni,nk,0) > 0) {
          /* Unlabelled object point found */
          nPt.x = nj;
          nPt.y = ni;
          nPt.z = nk;
          MRIIvox(Label, nj, ni, nk) = label;
          myQueuePush(NeiQ, &nPt);
        }

      }
    }
  }

  myQueueDelete(NeiQ);
  return;

}

void GetLargestCC6(MRI *orivol) {
  /* This function keeps the largest CC, and reset all other CC to bgvalue (0) */
  MRI *Label;
  int i,j,k;
  int maxSize, maxLabel, curSize, curLabel;
  int minX,minY,minZ,maxX,maxY,maxZ;
  int XN, YN, ZN;
  POINTI Pt;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  Label = MRIalloc(XN, YN, ZN, MRI_INT);
  MRIcopyHeader(orivol, Label);

  for (i=0; i<YN;i++)
    for (j=0; j< XN;j++)
      for (k=0;k<ZN;k++)
        MRIIvox(Label, j, i, k) = 0;

  curLabel = 1;
  maxSize = 0;
  maxLabel = 1;

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        if (MRIgetVoxVal(orivol, j, i, k, 0) > 0 && MRIIvox(Label,j,i,k) == 0) {
          Pt.x = j;
          Pt.y = i;
          Pt.z = k;
          GrassFire6(orivol,Label,curLabel,&Pt,&curSize,&minX,&maxX,&minY,&maxY,&minZ,&maxZ);
          if (maxSize < curSize) {
            maxSize = curSize;
            maxLabel = curLabel;
          }
          curLabel++;
        }
      }

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        if (MRIIvox(Label,j,i,k) != maxLabel)
          MRIsetVoxVal(orivol, j, i, k, 0, 0);
      }

  MRIfree(&Label);
  return;
}

void GetLargestCC18(MRI *orivol) {
  /* This function keeps the largest CC, and reset all other CC to bgvalue (0) */
  MRI *Label;
  int i,j,k;
  int maxSize, maxLabel, curSize, curLabel;
  int minX,minY,minZ,maxX,maxY,maxZ;
  int XN, YN, ZN;
  POINTI Pt;

  XN = orivol->width;
  YN = orivol->height;
  ZN = orivol->depth;

  Label = MRIalloc(XN, YN, ZN, MRI_INT);
  MRIcopyHeader(orivol, Label);

  for (i=0; i<YN;i++)
    for (j=0; j< XN;j++)
      for (k=0;k<ZN;k++)
        MRIIvox(Label, j, i, k) = 0;

  curLabel = 1;
  maxSize = 0;
  maxLabel = 1;

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        if (MRIgetVoxVal(orivol, j, i, k, 0) > 0 && MRIIvox(Label,j,i,k) == 0) {
          Pt.x = j;
          Pt.y = i;
          Pt.z = k;
          GrassFire18(orivol,Label,curLabel,&Pt,&curSize,&minX,&maxX,&minY,&maxY,&minZ,&maxZ);
          if (maxSize < curSize) {
            maxSize = curSize;
            maxLabel = curLabel;
          }
          curLabel++;
        }
      }

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        if (MRIIvox(Label,j,i,k) != maxLabel)
          MRIsetVoxVal(orivol, j, i, k, 0, 0);
      }

  MRIfree(&Label);
  return;
}

MRI * Dilation6(MRI *ori, MRI *out, int R) {
  int i,j,k,index,ci,cj,ck,count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  MRIcopyHeader(ori, tmpvol);

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori,j, i, k,0);
      }

  for (count=1; count<=R;count++) {
    for (i=0; i<YN;i++)
      for (j=0;j<XN;j++)
        for (k=0;k<ZN;k++) {
          if (MRIvox(tmpvol,j,i,k) == 1) {
            MRIsetVoxVal(out, j,i,k,0, 1);
            continue;
          }

          MRIsetVoxVal(out, j,i,k,0, 0);

          for (index=0;index < 6; index++) {
            ci = i + yoff6[index];
            cj = j + xoff6[index];
            ck = k + zoff6[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN)
              continue;
            if (MRIvox(tmpvol, cj, ci, ck) == 1) {
              MRIsetVoxVal(out, j,i,k,0, 1);
              break;
            }
          }
        }

    if (count < R) {
      for (i=0; i<YN;i++)
        for (j=0;j<XN;j++)
          for (k=0;k<ZN;k++) {
            MRIvox(tmpvol,j,i,k) = MRIgetVoxVal(out,j,i,k,0);
          }
    }

  }

  MRIfree(&tmpvol);
  return (out);
}

MRI * Erosion6(MRI *ori, MRI *out, int R) {
  int i,j,k,index,ci,cj,ck,count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  MRIcopyHeader(ori, tmpvol);

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori,j, i, k,0);
      }

  for (count=1; count<=R;count++) {
    for (i=0; i<YN;i++)
      for (j=0;j<XN;j++)
        for (k=0;k<ZN;k++) {
          if (MRIvox(tmpvol,j,i,k) == 0) {
            MRIsetVoxVal(out, j,i,k,0, 0);
            continue;
          }

          MRIsetVoxVal(out, j,i,k,0, 1);

          for (index=0;index < 6; index++) {
            ci = i + yoff6[index];
            cj = j + xoff6[index];
            ck = k + zoff6[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN)
              continue;

            if (MRIvox(tmpvol, cj, ci, ck) == 0) {
              MRIsetVoxVal(out, j,i,k,0, 0);
              break;
            }
          }
        }

    if (count < R) {
      for (i=0; i<YN;i++)
        for (j=0;j<XN;j++)
          for (k=0;k<ZN;k++) {
            MRIvox(tmpvol,j,i,k) = MRIgetVoxVal(out,j,i,k,0);
          }
    }

  }

  MRIfree(&tmpvol);
  return (out);
}


MRI * Dilation26(MRI *ori, MRI *out, int R) {
  int i,j,k,index,ci,cj,ck,count, XN, YN, ZN;
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

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori,j, i, k,0);
      }

  for (count=1; count<=R;count++) {
    for (i=0; i<YN;i++)
      for (j=0;j<XN;j++)
        for (k=0;k<ZN;k++) {
          if (MRIvox(tmpvol,j,i,k) == 1) {
            MRIsetVoxVal(out, j,i,k,0, 1);
            continue;
          }

          MRIsetVoxVal(out, j,i,k,0, 0);

          for (index=0;index < 26; index++) {
            ci = i + yoff26[index];
            cj = j + xoff26[index];
            ck = k + zoff26[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN)
              continue;
            if (MRIvox(tmpvol, cj, ci, ck) == 1) {
              MRIsetVoxVal(out, j,i,k,0, 1);
              break;
            }
          }
        }

    if (count < R) {
      for (i=0; i<YN;i++)
        for (j=0;j<XN;j++)
          for (k=0;k<ZN;k++) {
            MRIvox(tmpvol,j,i,k) = MRIgetVoxVal(out,j,i,k,0);
          }
    }

  }

  MRIfree(&tmpvol);
  return (out);
}

MRI * Erosion26(MRI *ori, MRI *out, int R) {
  int i,j,k,index,ci,cj,ck,count, XN, YN, ZN;
  MRI *tmpvol;

  XN = ori->width;
  YN = ori->height;
  ZN = ori->depth;

  if (out == NULL) {
    out = MRIclone(ori, NULL);
  }

  tmpvol = MRIalloc(XN, YN, ZN, MRI_UCHAR);
  MRIcopyHeader(ori, tmpvol);

  for (i=0; i<YN;i++)
    for (j=0;j<XN;j++)
      for (k=0;k<ZN;k++) {
        MRIvox(tmpvol, j, i, k) = MRIgetVoxVal(ori,j, i, k,0);
      }

  for (count=1; count<=R;count++) {
    for (i=0; i<YN;i++)
      for (j=0;j<XN;j++)
        for (k=0;k<ZN;k++) {
          if (MRIvox(tmpvol,j,i,k) == 0) {
            MRIsetVoxVal(out, j,i,k,0, 0);
            continue;
          }

          MRIsetVoxVal(out, j,i,k,0, 1);

          for (index=0;index < 26; index++) {
            ci = i + yoff26[index];
            cj = j + xoff26[index];
            ck = k + zoff26[index];
            if (ci < 0 || ci >= YN || cj < 0 || cj >= XN || ck < 0 || ck >= ZN)
              continue;

            if (MRIvox(tmpvol, cj, ci, ck) == 0) {
              MRIsetVoxVal(out, j,i,k,0, 0);
              break;
            }
          }
        }

    if (count < R) {
      for (i=0; i<YN;i++)
        for (j=0;j<XN;j++)
          for (k=0;k<ZN;k++) {
            MRIvox(tmpvol,j,i,k) = MRIgetVoxVal(out,j,i,k,0);
          }
    }

  }

  MRIfree(&tmpvol);
  return (out);
}

MRI * BinaryOpen6(MRI *ori, MRI *out, int R) {
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;

  tmpvol = Erosion6(ori, tmpvol, R);
  out = Dilation6(tmpvol,out, R);

  MRIfree(&tmpvol);

  return (out);
}

MRI * BinaryOpen26(MRI *ori, MRI *out, int R) {
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;

  printf("Erosion...\n");
  tmpvol = Erosion26(ori, tmpvol, R);
  printf("Dilation...\n");
  out = Dilation26(tmpvol,out, R);

  MRIfree(&tmpvol);

  return (out);
}


MRI  * BinaryClose6(MRI *ori, MRI *out, int R) {
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;


  tmpvol = Dilation6(ori, tmpvol, R);
  out = Erosion6(tmpvol,out, R);

  MRIfree(&tmpvol);

  return (out);
}


MRI * BinaryClose26(MRI *ori, MRI *out, int R) {
  /* This function assumes input volume is binary. It performs openning
     on the input volume using a structure element defined as R-times
     convolution of the basic SE (6, or 18).
   */

  MRI *tmpvol = NULL;


  tmpvol = Dilation26(ori, tmpvol, R);
  out = Erosion26(tmpvol,out, R);

  MRIfree(&tmpvol);

  return (out);
}
