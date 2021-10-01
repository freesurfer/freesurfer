/*
 * Original Author: Doug Greve
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
#include <string.h>

#include <numerics.h>
#include "version.h"
#include "diag.h"
#include "matrix.h"
#include "mri.h"
#include "randomfields.h"
#include "resample.h"
#include "transform.h"
#include "utils.h"
#define VOLCLUSTER_SRC
#include "surfcluster.h"
#include "volcluster.h"


static int ConvertCRS2XYZ(int col, int row, int slc, MATRIX *CRS2XYZ, float *x, float *y, float *z);

/*----------------------------------------------------------------*/
VOLCLUSTER *clustAllocCluster(int nmembers)
{
  VOLCLUSTER *vc;

  vc = (VOLCLUSTER *)calloc(1, sizeof(VOLCLUSTER));

  if (nmembers == 0) return (vc);

  vc->nmembers = nmembers;
  vc->col = (int *)calloc(nmembers, sizeof(int));
  vc->row = (int *)calloc(nmembers, sizeof(int));
  vc->slc = (int *)calloc(nmembers, sizeof(int));
  vc->x = (float *)calloc(nmembers, sizeof(float));
  vc->y = (float *)calloc(nmembers, sizeof(float));
  vc->z = (float *)calloc(nmembers, sizeof(float));

  return (vc);
}

/*----------------------------------------------------------------*/
VOLCLUSTER **clustAllocClusterList(int nlist)
{
  VOLCLUSTER **vclist;
  vclist = (VOLCLUSTER **)calloc(nlist, sizeof(VOLCLUSTER *));
  if (vclist == NULL) {
    fprintf(stderr, "ERROR: clustAllocClusterList: could not alloc %d\n", nlist);
    return (NULL);
  }
  return (vclist);
}

/*----------------------------------------------------------------*/
int clustFreeCluster(VOLCLUSTER **ppvc)
{
  VOLCLUSTER *vc;
  vc = *ppvc;

  if (vc->col != NULL) {
    free(vc->col);
    vc->col = NULL;
  }
  if (vc->row != NULL) {
    free(vc->row);
    vc->row = NULL;
  }
  if (vc->slc != NULL) {
    free(vc->slc);
    vc->slc = NULL;
  }

  if (vc->x != NULL) {
    free(vc->x);
    vc->x = NULL;
  }
  if (vc->y != NULL) {
    free(vc->y);
    vc->y = NULL;
  }
  if (vc->z != NULL) {
    free(vc->z);
    vc->z = NULL;
  }

  free(*ppvc);
  *ppvc = NULL;

  return (0);
}

/*------------------------------------------------------------------------*/
int clustFreeClusterList(VOLCLUSTER ***pppvclist, int nlist)
{
  int n;
  VOLCLUSTER **vclist;
  vclist = *pppvclist;

  if (nlist == 0) return (0);

  for (n = 0; n < nlist; n++)
    if (vclist[n] != NULL) clustFreeCluster(&vclist[n]);

  free(**pppvclist);
  **pppvclist = NULL;

  return (0);
}

/*------------------------------------------------------------------------*/
int clustDumpCluster(FILE *fp, VOLCLUSTER *vc, MRI *vol, int frame)
{
  int n;
  float val;

  for (n = 0; n < vc->nmembers; n++) {
    fprintf(fp, "%4d  %3d %3d %3d ", n, vc->col[n], vc->row[n], vc->slc[n]);
    fprintf(fp, "%7.2f %7.2f %7.2f ", vc->x[n], vc->y[n], vc->z[n]);
    if (vol != NULL) {
      val = MRIgetVoxVal(vol, vc->col[n], vc->row[n], vc->slc[n], frame);
      fprintf(fp, "%12.5f\n", val);
    }
    else
      fprintf(fp, "\n");
  }
  return (0);
}

/*------------------------------------------------------------------------*/
int clustDumpClusterList(FILE *fp, VOLCLUSTER **vclist, int nlist, MRI *vol, int frame)
{
  int n;

  for (n = 0; n < nlist; n++) {
    fprintf(fp,
            "%3d %5d %5d %g-------------------------------\n",
            n,
            vclist[n]->nmembers,
            vclist[n]->maxmember,
            vclist[n]->maxval);
    clustDumpCluster(fp, vclist[n], vol, frame);
  }

  return (0);
}

/*------------------------------------------------------------------------*/
int clustValueInRange(float val, float thmin, float thmax, int thsign)
{
  if (thsign == 0) val = fabs(val);
  if (thsign == -1) val = -val;

  if (thmax > 0) {
    if (val >= thmin && val <= thmax) return 1;
  }
  else {
    if (val >= thmin) return 1;
  }

  return (0);
}

/*------------------------------------------------------------------------*/
MRI *clustInitHitMap(MRI *vol,
                     int frame,
                     float thmin,
                     float thmax,
                     int thsign,
                     int *nhits,
                     int **hitcol,
                     int **hitrow,
                     int **hitslc,
                     MRI *binmask,
                     int maskframe)
{
  MRI *HitMap;
  int row, col, slc;
  float val;
  int nh, *hrow, *hcol, *hslc;
  int maskval;

  if (Gdiag_no > 0) {
    printf("clustInitHitMap: \n");
    printf("   thmin  = %f\n", thmin);
    printf("   thmax  = %f\n", thmax);
    printf("   thsign = %d\n", thsign);
    printf("   maskframe = %d\n", maskframe);
  }

  /* allocate the hit map */
  HitMap = MRIalloc(vol->width, vol->height, vol->depth, MRI_INT);
  if (HitMap == NULL) {
    printf("ERROR: clustInitHitMap: could not alloc HitMap\n");
    return (NULL);
  }

  /* count the number of hits */
  nh = 0;
  for (col = 0; col < vol->width; col++) {
    for (row = 0; row < vol->height; row++) {
      for (slc = 0; slc < vol->depth; slc++) {
        if (binmask != NULL) {
          maskval = MRIgetVoxVal(binmask, col, row, slc, maskframe);
          if (maskval == 0) continue;
        }
        val = MRIgetVoxVal(vol, col, row, slc, frame);
        if (clustValueInRange(val, thmin, thmax, thsign)) nh++;
      }
    }
  }
  // printf("INFO: clustInitHitMap: found %d hits\n", nh );

  /* check that there are hits */
  if (nh == 0) {
    if (Gdiag_no > 0) printf("no hits in hit map\n");
    *nhits = 0;
    return (HitMap);
  }

  /* allocate the rows cols and slices */
  hcol = (int *)calloc(nh, sizeof(int));
  if (hcol == NULL) {
    printf("ERROR: clustInitHitMap: could not alloc hit cols\n");
    MRIfree(&HitMap);
    return (NULL);
  }
  hrow = (int *)calloc(nh, sizeof(int));
  if (hrow == NULL) {
    printf("ERROR: clustInitHitMap: could not alloc hit rows\n");
    free(hcol);
    MRIfree(&HitMap);
    return (NULL);
  }
  hslc = (int *)calloc(nh, sizeof(int));
  if (hslc == NULL) {
    printf("ERROR: clustInitHitMap: could not alloc hit slices\n");
    free(hcol);
    free(hrow);
    MRIfree(&HitMap);
    return (NULL);
  }

  /* Now go back through the volume to assign values */
  nh = 0;
  for (col = 0; col < vol->width; col++) {
    for (row = 0; row < vol->height; row++) {
      for (slc = 0; slc < vol->depth; slc++) {
        val = MRIgetVoxVal(vol, col, row, slc, frame);
        if (binmask != NULL) {
          maskval = MRIgetVoxVal(binmask, col, row, slc, maskframe);
          // printf("%2d %2d %2d  %d %g\n",col,row,slc,maskval,val);
          if (maskval == 0) {
            MRIsetVoxVal(HitMap, col, row, slc, 0, 1);
            continue;
          }
        }
        if (clustValueInRange(val, thmin, thmax, thsign)) {
          hcol[nh] = col;
          hrow[nh] = row;
          hslc[nh] = slc;
          MRIsetVoxVal(HitMap, col, row, slc, 0, 0);
          nh++;
        }
        else
          MRIsetVoxVal(HitMap, col, row, slc, 0, 1);
      }
    }
  }
  if (Gdiag_no > 1) printf("INFO: clustInitHitMap: found %d hits\n", nh);

  *hitcol = hcol;
  *hitrow = hrow;
  *hitslc = hslc;
  *nhits = nh;

  return (HitMap);
}

/*------------------------------------------------------------------------*/
int clustAddMember(VOLCLUSTER *vc, int col, int row, int slc)
{
  int nmemb;
  int *tmp;
  float *ftmp;

  nmemb = vc->nmembers;

  tmp = (int *)realloc(vc->col, (nmemb + 1) * sizeof(int));
  if (tmp == NULL) {
    printf("ERROR: clustAddMember: could not alloc %d\n", nmemb + 1);
    return (1);
  }
  vc->col = tmp;

  tmp = (int *)realloc(vc->row, (nmemb + 1) * sizeof(int));
  if (tmp == NULL) {
    printf("ERROR: clustAddMember: could not alloc %d\n", nmemb + 1);
    return (1);
  }
  vc->row = tmp;

  tmp = (int *)realloc(vc->slc, (nmemb + 1) * sizeof(int));
  if (tmp == NULL) {
    printf("ERROR: clustAddMember: could not alloc %d\n", nmemb + 1);
    return (1);
  }
  vc->slc = tmp;

  ftmp = (float *)realloc(vc->x, (nmemb + 1) * sizeof(float));
  if (ftmp == NULL) {
    printf("ERROR: clustAddMember: could not alloc %d\n", nmemb + 1);
    return (1);
  }
  vc->x = ftmp;

  ftmp = (float *)realloc(vc->y, (nmemb + 1) * sizeof(float));
  if (ftmp == NULL) {
    printf("ERROR: clustAddMember: could not alloc %d\n", nmemb + 1);
    return (1);
  }
  vc->y = ftmp;

  ftmp = (float *)realloc(vc->z, (nmemb + 1) * sizeof(float));
  if (ftmp == NULL) {
    printf("ERROR: clustAddMember: could not alloc %d\n", nmemb + 1);
    return (1);
  }
  vc->z = ftmp;

  vc->col[nmemb] = col;
  vc->row[nmemb] = row;
  vc->slc[nmemb] = slc;

  vc->nmembers = nmemb + 1;

  return (0);
}

/*------------------------------------------------------------------------*/
int clustGrowOneVoxel(VOLCLUSTER *vc, int col0, int row0, int slc0, MRI *HitMap, int AllowDiag)
{
  int col, row, slc;
  int dcol, drow, dslc, dsum;
  int nadded;

  nadded = 0;
  for (dcol = -1; dcol <= +1; dcol++) {
    for (drow = -1; drow <= +1; drow++) {
      for (dslc = -1; dslc <= +1; dslc++) {
        // Check for neighbor beyond the edge-of-volume
        col = col0 + dcol;
        if (col < 0 || col >= HitMap->width) continue;
        row = row0 + drow;
        if (row < 0 || row >= HitMap->height) continue;
        slc = slc0 + dslc;
        if (slc < 0 || slc >= HitMap->depth) continue;

        // If not allowing for diagonal connections
        if (!AllowDiag) {
          dsum = fabs(dcol) + fabs(drow) + fabs(dslc);
          if (dsum != 1) continue;
        }

        if (MRIgetVoxVal(HitMap, col, row, slc, 0)) continue;
        // printf("Adding %3d %3d %3d\n",col,row,slc);

        clustAddMember(vc, col, row, slc);
        MRIsetVoxVal(HitMap, col, row, slc, 0, 1);
        nadded++;
      }
    }
  }
  // printf("GrowOneFrom: %3d %3d %3d  %d\n",col0,row0,slc0,nadded);

  return (nadded);
}

/*------------------------------------------------------------------------*/
VOLCLUSTER *clustGrow(int col0, int row0, int slc0, MRI *HitMap, int AllowDiag, int npassesmax)
{
  VOLCLUSTER *vc;
  int nthmember, nmembers_now;
  int col, row, slc;
  int nadded, nthpass, n;

  vc = (VOLCLUSTER *)calloc(1, sizeof(VOLCLUSTER));

  /* put the seed point in the cluster */
  clustAddMember(vc, col0, row0, slc0);
  MRIsetVoxVal(HitMap, col0, row0, slc0, 0, 1);
  vc->voxsize = HitMap->xsize * HitMap->ysize * HitMap->zsize;

  nthpass = 0;
  nadded = 1;
  while (nadded > 0) {
    // printf("%4d  %5d  %d\n",nthpass,vc->nmembers,nadded);

    nadded = 0;
    nmembers_now = vc->nmembers;
    for (nthmember = 0; nthmember < nmembers_now; nthmember++) {
      col = vc->col[nthmember];
      row = vc->row[nthmember];
      slc = vc->slc[nthmember];
      n = clustGrowOneVoxel(vc, col, row, slc, HitMap, AllowDiag);
      // printf("Grown: %3d %3d %3d  %d\n",col,row,slc,n);
      nadded += n;
    }

    nthpass++;
    if(npassesmax > 0 && nthpass > npassesmax) break;
  }

  return (vc);
}

/*-------------------------------------------------------------------*/
int clustMaxMember(VOLCLUSTER *vc, MRI *vol, int frame, int thsign)
{
  int n;
  float val = 0.0, val0;

  vc->maxval = 0.0;
  for (n = 0; n < vc->nmembers; n++) {
    val0 = MRIgetVoxVal(vol, vc->col[n], vc->row[n], vc->slc[n], frame);
    if (thsign == 1) val = val0;
    if (thsign == 0) val = fabs(val0);
    if (thsign == -1) val = -val0;
    if (fabs(vc->maxval) < val) {
      vc->maxval = val0; /* keep orginal signed value */
      vc->maxmember = n;
    }
  }

  return (0);
}

/*------------------------------------------------------------------------*/
VOLCLUSTER **clustPruneBySize(VOLCLUSTER **vclist, int nlist, float voxsize, float sizethresh, int *nkeep)
{
  VOLCLUSTER **vcprune;
  int n;
  float clustersize;

  /* count the number to keep */
  *nkeep = 0;
  for (n = 0; n < nlist; n++) {
    clustersize = vclist[n]->nmembers * voxsize;
    if (clustersize >= sizethresh) (*nkeep)++;
  }

  vcprune = (VOLCLUSTER **)calloc(*nkeep, sizeof(VOLCLUSTER *));

  *nkeep = 0;
  for (n = 0; n < nlist; n++) {
    clustersize = vclist[n]->nmembers * voxsize;
    if (clustersize >= sizethresh) {
      vcprune[(*nkeep)] = clustCopyCluster(vclist[n]);
      (*nkeep)++;
    }
  }

  return (vcprune);
}

/*------------------------------------------------------------------------*/
VOLCLUSTER **clustPruneByDistance(VOLCLUSTER **vclist, int nlist, float distthresh, int *nkeep)
{
  VOLCLUSTER **vcprune;
  int n1, n2, nmax1, nmax2;
  float max1, max2;
  int keep;
  float x1, y1, z1;
  float x2, y2, z2;
  int pass;
  float d;

  vcprune = NULL;

  /* Two passes: (1) counts the number for alloc, (2) copies clusters */
  for (pass = 1; pass <= 2; pass++) {
    if (pass == 2) {
      if (*nkeep == 0) return (NULL);
      vcprune = (VOLCLUSTER **)calloc(*nkeep, sizeof(VOLCLUSTER *));
    }

    /*-- Go through each cluster -- */
    *nkeep = 0;
    for (n1 = 0; n1 < nlist; n1++) {
      nmax1 = vclist[n1]->maxmember;
      max1 = vclist[n1]->maxval;
      x1 = vclist[n1]->x[nmax1];
      y1 = vclist[n1]->y[nmax1];
      z1 = vclist[n1]->z[nmax1];

      /*-- Compare to every other cluster -- */
      keep = 1;
      for (n2 = 0; n2 < nlist; n2++) {
        if (n1 == n2) continue; /* dont compare to self */

        nmax2 = vclist[n2]->maxmember;
        max2 = vclist[n2]->maxval;
        x2 = vclist[n2]->x[nmax2];
        y2 = vclist[n2]->y[nmax2];
        z2 = vclist[n2]->z[nmax2];

        /* Compute the distance from the max of one to the max of the
           other */
        d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

        /* If the distance is less than threshold and the max of the
           first is less than the max of the second, throw out the
           first (dont worry about the second here */
        if (d < distthresh && max1 < max2) {
          // printf("Pruning %d: (%5.2f %5.2f %5.2f) (%5.2f %5.2f %5.2f) %g\n",
          // n1,x1,y1,z1,x2,y2,z2,d);
          keep = 0;
          break;
        }

      } /* end n2 loop */

      if (keep) {
        if (pass == 2) vcprune[(*nkeep)] = clustCopyCluster(vclist[n1]);
        (*nkeep)++;
      }

    } /* end n1 loop */

  } /* end pass loop */

  return (vcprune);
}

/*------------------------------------------------------------------------
  Note: Prunes by pvalue so must be LESS than cwpvalthresh  NOT -log10(p)
------------------------------------------------------------------------*/
VOLCLUSTER **clustPruneByCWPval(VOLCLUSTER **vclist, int nlist, double cwpvalthresh, int *nkeep)
{
  VOLCLUSTER **vcprune;
  int n;

  /* count the number to keep -- Note: pvalue, NOT -log10(p)*/
  *nkeep = 0;
  for (n = 0; n < nlist; n++)
    if (vclist[n]->pval_clusterwise <= cwpvalthresh) (*nkeep)++;
  vcprune = (VOLCLUSTER **)calloc(*nkeep, sizeof(VOLCLUSTER *));

  *nkeep = 0;
  for (n = 0; n < nlist; n++) {
    if (vclist[n]->pval_clusterwise <= cwpvalthresh) {
      vcprune[(*nkeep)] = clustCopyCluster(vclist[n]);
      (*nkeep)++;
    }
  }

  return (vcprune);
}

/*----------------------------------------------------------------*/
VOLCLUSTER *clustCopyCluster(VOLCLUSTER *vc)
{
  VOLCLUSTER *vc2;
  int ncopy;

  vc2 = clustAllocCluster(vc->nmembers);

  ncopy = vc->nmembers * sizeof(int);
  memmove(vc2->col, vc->col, ncopy);
  memmove(vc2->row, vc->row, ncopy);
  memmove(vc2->slc, vc->slc, ncopy);

  ncopy = vc->nmembers * sizeof(float);
  memmove(vc2->x, vc->x, ncopy);
  memmove(vc2->y, vc->y, ncopy);
  memmove(vc2->z, vc->z, ncopy);

  vc2->maxmember = vc->maxmember;
  vc2->maxval = vc->maxval;
  vc2->voxsize = vc->voxsize;

  vc2->pval_clusterwise = vc->pval_clusterwise;
  vc2->pval_clusterwise_low = vc->pval_clusterwise_low;
  vc2->pval_clusterwise_hi = vc->pval_clusterwise_hi;

  return (vc2);
}

/*----------------------------------------------------------------*/
VOLCLUSTER **clustCopyClusterList(VOLCLUSTER **vclist, int nlist, VOLCLUSTER **vclist2)
{
  int n;

  if (vclist2 == NULL) vclist2 = clustAllocClusterList(nlist);

  for (n = 0; n < nlist; n++) {
    vclist2[n] = clustCopyCluster(vclist[n]);
  }

  return (vclist2);
}

/*----------------------------------------------------------------*/
int clustCompareCluster(const void *a, const void *b)
{
  VOLCLUSTER *vc1, *vc2;

  vc1 = *((VOLCLUSTER **)a);
  vc2 = *((VOLCLUSTER **)b);

  // sort by most significant
  if (vc1->pval_clusterwise < vc2->pval_clusterwise) return (-1);
  if (vc1->pval_clusterwise > vc2->pval_clusterwise) return (+1);

  // sort by extent
  if (vc1->nmembers > vc2->nmembers) return (-1);
  if (vc1->nmembers < vc2->nmembers) return (+1);

  // now sort by maxval
  if (fabs(vc1->maxval) > fabs(vc2->maxval)) return (-1);
  if (fabs(vc1->maxval) < fabs(vc2->maxval)) return (+1);

  return (0);
}

/*----------------------------------------------------------------*/
VOLCLUSTER **clustSortClusterList(VOLCLUSTER **vclist, int nlist, VOLCLUSTER **vcsorted)
{
  if (vcsorted == NULL) vcsorted = clustAllocClusterList(nlist);

  if (vclist != vcsorted) clustCopyClusterList(vclist, nlist, vcsorted);

  qsort((void *)vcsorted, nlist, sizeof(VOLCLUSTER **), clustCompareCluster);

  return (vcsorted);
}

/*----------------------------------------------------------------
  clustComputeXYZ() - computes the xyz coordinate of each member
  of a cluster given the 4x4 matrix that transforms the col, row,
  and slice into x, y, and z.
  ----------------------------------------------------------------*/
int clustComputeXYZ(VOLCLUSTER *vc, MATRIX *CRS2XYZ)
{
  int n;

  for (n = 0; n < vc->nmembers; n++) {
    ConvertCRS2XYZ(vc->col[n], vc->row[n], vc->slc[n], CRS2XYZ, &(vc->x[n]), &(vc->y[n]), &(vc->z[n]));
  }

  return (0);
}

/*----------------------------------------------------------------
  clustComputeTal() - computes the talairach xyz coordinate of each
  member of a cluster given the 4x4 matrix that transforms the col,
  row, and slice into MNI coorinates. The MNI coordinates are
  transformed into talairach coordinates using a piece-wise linear
  transformation. See FixMNITal in transforms.c.
  ----------------------------------------------------------------*/
int clustComputeTal(VOLCLUSTER *vc, MATRIX *CRS2MNI)
{
  int n;

  for (n = 0; n < vc->nmembers; n++) {
    ConvertCRS2XYZ(vc->col[n], vc->row[n], vc->slc[n], CRS2MNI, &(vc->x[n]), &(vc->y[n]), &(vc->z[n]));
    if (VolClustFixMNI) FixMNITal(vc->x[n], vc->y[n], vc->z[n], &(vc->x[n]), &(vc->y[n]), &(vc->z[n]));
  }

  return (0);
}

/*----------------------------------------------------------------*/
MRI *clustClusterList2Vol(VOLCLUSTER **vclist, int nlist, MRI *tvol, int frame, int ValOption)
{
  MRI *vol;
  int nthvc, n, f;
  VOLCLUSTER *vc;
  float val;

  // Allocate vol and fill with all zeros
  if (ValOption)
    vol = MRIallocSequence(tvol->width, tvol->height, tvol->depth, tvol->type, tvol->nframes);
  else
    vol = MRIallocSequence(tvol->width, tvol->height, tvol->depth, tvol->type, 1);
  MRIcopyHeader(tvol, vol);

  // Go thru each cluster and assign values to the voxels in clusters
  for (nthvc = 0; nthvc < nlist; nthvc++) {
    vc = vclist[nthvc];
    for (n = 0; n < vc->nmembers; n++) {
      if (ValOption == 1) {
        // Copy the values from the input
        for (f = 0; f < tvol->nframes; f++) {
          val = MRIgetVoxVal(tvol, vc->col[n], vc->row[n], vc->slc[n], f);
          MRIsetVoxVal(vol, vc->col[n], vc->row[n], vc->slc[n], f, val);
        }
      }
      else
        // Set values to the cluster number
        MRIsetVoxVal(vol, vc->col[n], vc->row[n], vc->slc[n], 0, nthvc + 1);
    }
  }
  return (vol);
}

/*----------------------------------------------------------------*/
LABEL *clustCluster2Label(
    VOLCLUSTER *vc, MRI *vol, int frame, float colres, float rowres, float sliceres, MATRIX *FSA2Func)
{
  LABEL *label;
  MATRIX *CRS2Func, *Func2FSA;
  MATRIX *xyzFunc, *xyzFSA;
  int n;
  float val;

  CRS2Func = MRIxfmCRS2XYZtkreg(vol);
  clustComputeXYZ(vc, CRS2Func);

  /* Matrix to convert from volume to anatomical */
  Func2FSA = MatrixInverse(FSA2Func, NULL);

  /* Preallocate vectors */
  xyzFunc = MatrixAlloc(4, 1, MATRIX_REAL);
  xyzFunc->rptr[4][1] = 1;
  xyzFSA = MatrixAlloc(4, 1, MATRIX_REAL);

  /* Alloc the label */
  label = LabelAlloc(vc->nmembers, NULL, NULL);

  /* Assign label values */
  for (n = 0; n < vc->nmembers; n++) {
    val = MRIgetVoxVal(vol, vc->col[n], vc->row[n], vc->slc[n], frame);

    /* convert Functional XYZ FSA XYZ */
    xyzFunc->rptr[1][1] = vc->x[n];
    xyzFunc->rptr[2][1] = vc->y[n];
    xyzFunc->rptr[3][1] = vc->z[n];
    MatrixMultiply(Func2FSA, xyzFunc, xyzFSA);

    /* assign fields to label */
    label->lv[n].x = rint(xyzFSA->rptr[1][1]);
    label->lv[n].y = rint(xyzFSA->rptr[2][1]);
    label->lv[n].z = rint(xyzFSA->rptr[3][1]);
    label->lv[n].stat = val;
    label->lv[n].vno = -1;
  }
  label->n_points = vc->nmembers;

  MatrixFree(&xyzFunc);
  MatrixFree(&xyzFSA);
  MatrixFree(&Func2FSA);
  MatrixFree(&CRS2Func);

  return (label);
}

/*-------------------------------------------------------------*/
VOLCLUSTER **clustGetClusters(MRI *vol,
                              int frame,
                              float threshmin,
                              float threshmax,
                              int threshsign,
                              float minclustsizemm3,
                              MRI *binmask,
                              int *nClusters,
                              MATRIX *XFM)
{
  int nthhit, nclusters, nhits, *hitcol = NULL, *hitrow = NULL, *hitslc = NULL;
  int col, row, slc, allowdiag = 0, nprunedclusters;
  MRI *HitMap;
  VOLCLUSTER **ClusterList, **ClusterList2;
  float voxsizemm3, distthresh = 0;

  voxsizemm3 = vol->xsize * vol->ysize * vol->zsize;

  /* Initialize the hit map - this is a map of voxels that have been
     accounted for as either outside of a the threshold range or
     belonging to a cluster. The initialization is for thresholding */
  HitMap = clustInitHitMap(vol, frame, threshmin, threshmax, threshsign, &nhits, &hitcol, &hitrow, &hitslc, binmask, 0);
  if (HitMap == NULL) {
    // Nothing survived the first thresholding
    *nClusters = 0;
    return (NULL);
  }
  if (Gdiag_no > 0) printf("INFO: Found %d voxels in threhold range\n", nhits);

  /* Allocate an array of clusters equal to the number of hits -- this
     is the maximum number of clusters possible */
  ClusterList = clustAllocClusterList(nhits);
  if (ClusterList == NULL) {
    printf("ERROR: could not alloc %d clusters\n", nhits);
    return (NULL);
  }

  nclusters = 0;
  for (nthhit = 0; nthhit < nhits; nthhit++) {
    /* Determine whether this hit is still valid. It may
       not be if it was assigned to the cluster of a
       previous hit */
    col = hitcol[nthhit];
    row = hitrow[nthhit];
    slc = hitslc[nthhit];
    if (MRIgetVoxVal(HitMap, col, row, slc, 0)) continue;

    /* Grow cluster using this hit as a seed */
    ClusterList[nclusters] = clustGrow(col, row, slc, HitMap, allowdiag, -1);
    ClusterList[nclusters]->voxsize = voxsizemm3;

    /* Determine the member with the maximum value */
    clustMaxMember(ClusterList[nclusters], vol, frame, threshsign);

    if (XFM) clustComputeTal(ClusterList[nclusters], XFM);

    /* increment the number of clusters */
    nclusters++;
  }
  free(hitcol);
  free(hitrow);
  free(hitslc);

  if (Gdiag_no > 0) printf("INFO: Found %d clusters that meet threshold criteria\n", nclusters);

  /* Remove clusters that do not meet the minimum size requirement */
  ClusterList2 = clustPruneBySize(ClusterList, nclusters, voxsizemm3, minclustsizemm3, &nprunedclusters);
  clustFreeClusterList(&ClusterList, nclusters);
  nclusters = nprunedclusters;
  ClusterList = ClusterList2;

  if (Gdiag_no > 0) printf("INFO: Found %d clusters that meet size criteria\n", nclusters);

  /* Remove clusters that do not meet the minimum distance requirement */
  if (distthresh > 0.0) {
    if (Gdiag_no > 0) printf("INFO: pruning by distance %g\n", distthresh);
    ClusterList2 = clustPruneByDistance(ClusterList, nclusters, distthresh, &nprunedclusters);
    clustFreeClusterList(&ClusterList, nclusters);
    nclusters = nprunedclusters;
    ClusterList = ClusterList2;
  }

  /* Sort Clusters by MaxValue */
  ClusterList2 = clustSortClusterList(ClusterList2, nclusters, NULL);
  clustFreeClusterList(&ClusterList, nclusters);
  ClusterList = ClusterList2;

  MRIfree(&HitMap);

  if (Gdiag_no > 0) printf("INFO: Found %d final clusters\n", nclusters);
  *nClusters = nclusters;
  return (ClusterList);
}

/*-----------------------------------------------------------------
  clustMaxClusterCount() - returns the voxel count of the cluster
  with the largest count. Multiply this by the voxel size
  to get the max cluster volume.
  -------------------------------------------------------------*/
int clustMaxClusterCount(VOLCLUSTER **VolClustList, int nClusters)
{
  int n, MaxCount = 0;

  for (n = 0; n < nClusters; n++)
    if (VolClustList[n]->nmembers > MaxCount) MaxCount = VolClustList[n]->nmembers;
  return (MaxCount);
}

/*------------------------------------------------------------------------*/
int clustDumpSummary(FILE *fp, VOLCLUSTER **ClusterList, int nClusters)
{
  int n, col, row, slc;
  float x, y, z;

  fprintf(fp, "--------------------------------------------------\n");
  for (n = 0; n < nClusters; n++) {
    col = ClusterList[n]->col[ClusterList[n]->maxmember];
    row = ClusterList[n]->row[ClusterList[n]->maxmember];
    slc = ClusterList[n]->slc[ClusterList[n]->maxmember];
    x = ClusterList[n]->x[ClusterList[n]->maxmember];
    y = ClusterList[n]->y[ClusterList[n]->maxmember];
    z = ClusterList[n]->z[ClusterList[n]->maxmember];
    fprintf(fp,
            "%3d %4d  %7.1f  %3d %3d %3d  %g  %g %g %g   \n",
            n + 1,
            ClusterList[n]->nmembers,
            ClusterList[n]->nmembers * ClusterList[n]->voxsize,
            col,
            row,
            slc,
            x,
            y,
            z,
            ClusterList[n]->maxval);
  }
  fprintf(fp, "--------------------------------------------------\n");
  return (0);
}

/*-------------------------------------------------------------*/

/*----------------------------------------------------------------
  ConvertCRS2XYZ() - computes the xyz coordinate given the CRS and
  the transform matrix. This function just hides the matrix
  operations.
  ----------------------------------------------------------------*/
static int ConvertCRS2XYZ(int col, int row, int slc, MATRIX *CRS2XYZ, float *x, float *y, float *z)
{
  MATRIX *crs, *xyz;

  crs = MatrixAlloc(4, 1, MATRIX_REAL);
  crs->rptr[1][1] = (float)col;
  crs->rptr[2][1] = (float)row;
  crs->rptr[3][1] = (float)slc;
  crs->rptr[4][1] = 1;

  xyz = MatrixMultiply(CRS2XYZ, crs, NULL);

  *x = xyz->rptr[1][1];
  *y = xyz->rptr[2][1];
  *z = xyz->rptr[3][1];

  MatrixFree(&crs);
  MatrixFree(&xyz);

  return (0);
}

/*----------------------------------------------------*/
CHT *CHTalloc(int n_ithr, double ithr_lo, double ithr_hi, int n_sthr, double sthr_lo, double sthr_hi)
{
  CHT *cht;
  double dithr, dsthr;
  int i, v;

  cht = (CHT *)calloc(1, sizeof(CHT));

  cht->n_ithr = n_ithr;
  cht->ithr = (double *)calloc(n_ithr, sizeof(double));
  if (n_ithr != 1)
    dithr = (ithr_hi - ithr_lo) / (n_ithr - 1);
  else
    dithr = 0;
  for (i = 0; i < n_ithr; i++) cht->ithr[i] = ithr_lo + dithr * i;
  cht->ithr_lo = ithr_lo;
  cht->ithr_hi = ithr_hi;

  cht->n_sthr = n_sthr;
  cht->sthr = (double *)calloc(n_sthr, sizeof(double));
  if (n_sthr != 1)
    dsthr = (sthr_hi - sthr_lo) / (n_sthr - 1);
  else
    dsthr = 0;
  for (v = 0; v < n_sthr; v++) cht->sthr[v] = sthr_lo + dsthr * v;
  cht->sthr_lo = sthr_lo;
  cht->sthr_hi = sthr_hi;

  cht->hits = (int **)calloc(n_ithr, sizeof(int *));
  for (i = 0; i < n_ithr; i++) cht->hits[i] = (int *)calloc(n_sthr, sizeof(int));

  return (cht);
}

/*----------------------------------------------------*/
int CHTfree(CHT **ppcht)
{
  CHT *cht;
  int i;

  cht = *ppcht;
  for (i = 0; i < cht->n_ithr; i++) free(cht->hits[i]);
  free(cht->hits);
  free(cht->ithr);
  free(cht->sthr);
  free(*ppcht);
  *ppcht = NULL;

  return (0);
}

/*----------------------------------------------------*/
int CHTprint(FILE *fp, CHT *cht)
{
  int i, v;

  fprintf(fp, "# CHT 1\n");
  fprintf(fp, "# nsim        %d\n", cht->nsim);
  fprintf(fp, "# seed        %ld\n", cht->seed);
  fprintf(fp, "# nvox        %d\n", cht->nvox);
  fprintf(fp, "# totsize     %lf\n", cht->totsize);
  fprintf(fp, "# fwhm        %lf\n", cht->fwhm);
  fprintf(fp, "# nsmooth     %d\n", cht->nsmooth);
  fprintf(fp, "# n_ithr      %d\n", cht->n_ithr);
  fprintf(fp, "# ithr_lo     %lf\n", cht->ithr_lo);
  fprintf(fp, "# ithr_hi     %lf\n", cht->ithr_hi);
  fprintf(fp, "# ithr_sign   %s\n", cht->ithr_sign);
  fprintf(fp, "# n_sthr      %d\n", cht->n_sthr);
  fprintf(fp, "# sthr_lo     %lf\n", cht->sthr_lo);
  fprintf(fp, "# sthr_hi     %lf\n", cht->sthr_hi);
  // fprintf(fp,"# STARTDATA\n");

  fprintf(fp, "     ");
  for (v = 0; v < cht->n_sthr; v++) fprintf(fp, "%4.2lf ", cht->sthr[v]);
  fprintf(fp, "\n");

  for (i = 0; i < cht->n_ithr; i++) {
    fprintf(fp, "%4.2lf ", cht->ithr[i]);
    for (v = 0; v < cht->n_sthr; v++) fprintf(fp, "%5d ", cht->hits[i][v]);
    fprintf(fp, "\n");
  }

  return (0);
}

/*----------------------------------------------------*/
int CHTwrite(char *fname, CHT *cht)
{
  FILE *fp;

  fp = fopen(fname, "w");
  if (fp == NULL) {
    printf("ERROR: could not open %s for writing\n", fname);
    return (1);
  }

  CHTprint(fp, cht);

  fclose(fp);

  return (0);
}

/*----------------------------------------------------*/
CHT *CHTread(char *fname)
{
  CHT *cht;
  FILE *fp;
  char tag[1000];
  char tmpstr[1000];
  int nsim;                /* number of simulation runs to generate table */
  long int seed;           // Seed for random number generator
  int nvox;                /* number of voxels/vertices in search area */
  double totsize;          /* total volume (mm^3) or area (mm^2) in search*/
  double fwhm;             /* fwhm in mm */
  int nsmooth;             /* number of smooth steps, surf only */
  double ithr_lo, ithr_hi; /* intensity threshold range */
  int n_ithr;              /* Number ithreshs bet lo and hi*/
  char ithr_sign[50];      /* abs, pos, neg*/
  double sthr_lo, sthr_hi; /* volume threshold range */
  int n_sthr;              /* Number sthreshs bet lo and hi*/
  int i, v;
  double dummy;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: could not open %s for reading\n", fname);
    return (NULL);
  }

  fgets(tmpstr, 1000, fp);
  if (tmpstr[0] != '#') {
    printf("ERROR: %s is not formatted correctly (#)\n", fname);
    return (NULL);
  }

  sscanf(tmpstr, "%*s %s", tag);
  if (strcmp(tag, "CHT")) {
    printf("ERROR: %s is not formatted correctly (CHT)\n", fname);
    return (NULL);
  }

  while (1) {
    // Grab a line
    fgets(tmpstr, 1000, fp);

    // Test whether still in the TAG section
    if (tmpstr[0] != '#') break;

    // Scan the tag
    sscanf(tmpstr, "%*s %s", tag);
    // printf("%s \n",tag);

    if (!strcmp(tag, "nsim")) sscanf(tmpstr, "%*s %*s %d", &nsim);
    if (!strcmp(tag, "seed")) sscanf(tmpstr, "%*s %*s %ld", &seed);
    if (!strcmp(tag, "nvox")) sscanf(tmpstr, "%*s %*s %d", &nvox);
    if (!strcmp(tag, "totsize")) sscanf(tmpstr, "%*s %*s %lf", &totsize);
    if (!strcmp(tag, "fwhm")) sscanf(tmpstr, "%*s %*s %lf", &fwhm);
    if (!strcmp(tag, "nsmooth")) sscanf(tmpstr, "%*s %*s %d", &nsmooth);
    if (!strcmp(tag, "n_ithr")) sscanf(tmpstr, "%*s %*s %d", &n_ithr);
    if (!strcmp(tag, "ithr_lo")) sscanf(tmpstr, "%*s %*s %lf", &ithr_lo);
    if (!strcmp(tag, "ithr_hi")) sscanf(tmpstr, "%*s %*s %lf", &ithr_hi);
    if (!strcmp(tag, "ithr_sign")) sscanf(tmpstr, "%*s %*s %s", ithr_sign);
    if (!strcmp(tag, "n_sthr")) sscanf(tmpstr, "%*s %*s %d", &n_sthr);
    if (!strcmp(tag, "sthr_lo")) sscanf(tmpstr, "%*s %*s %lf", &sthr_lo);
    if (!strcmp(tag, "sthr_hi")) sscanf(tmpstr, "%*s %*s %lf", &sthr_hi);
  }

  cht = CHTalloc(n_ithr, ithr_lo, ithr_hi, n_sthr, sthr_lo, sthr_hi);
  cht->nsim = nsim;
  cht->seed = seed;
  cht->nvox = nvox;
  cht->totsize = totsize;
  cht->fwhm = fwhm;
  cht->nsmooth = nsmooth;
  strcpy(cht->ithr_sign, ithr_sign);
  if (!strcasecmp(ithr_sign, "pos"))
    cht->ithr_signid = +1;
  else if (!strcasecmp(ithr_sign, "abs"))
    cht->ithr_signid = 0;
  else if (!strcasecmp(ithr_sign, "neg"))
    cht->ithr_signid = -1;
  else {
    printf("ERROR: ithr_sign = %s, unrecognized.\n", ithr_sign);
    printf("       ithr_sign must be pos, neg, or abs.\n");
    return (NULL);
  }

  // The first data line has already been skipped by the fgets above

  // Read in the hit table
  for (i = 0; i < cht->n_ithr; i++) {
    // Skip first column as it is the ithr for that row
    fscanf(fp, "%lf ", &dummy);
    for (v = 0; v < cht->n_sthr; v++) fscanf(fp, "%d ", &(cht->hits[i][v]));
  }

  fclose(fp);

  return (cht);
}

/*-----------------------------------------------------------------
  CHTcompare() - compares two CHTs. Returns 1 if they are different
  and 0 if they are the same. If an item in the targ is less than 0,
  then it's value is replaced with the value of the source. This
  does not trigger a 1 return and allows for values to be filled in.
  Does not compare seeds.
  ----------------------------------------------------------------*/
int CHTcompare(CHT *src, CHT *targ)
{
  if (targ->nvox < 0)
    targ->nvox = src->nvox;
  else if (targ->nvox != src->nvox)
    return (1);

  if (targ->totsize < 0)
    targ->totsize = src->totsize;
  else if (targ->totsize != src->totsize)
    return (1);

  if (targ->fwhm < 0)
    targ->fwhm = src->fwhm;
  else if (targ->fwhm != src->fwhm)
    return (1);

  if (targ->nsmooth < 0)
    targ->nsmooth = src->nsmooth;
  else if (targ->nsmooth != src->nsmooth)
    return (1);

  if (targ->ithr_lo < 0)
    targ->ithr_lo = src->ithr_lo;
  else if (targ->ithr_lo != src->ithr_lo)
    return (1);

  if (targ->ithr_hi < 0)
    targ->ithr_hi = src->ithr_hi;
  else if (targ->ithr_hi != src->ithr_hi)
    return (1);

  if (targ->n_ithr < 0)
    targ->n_ithr = src->n_ithr;
  else if (targ->n_ithr != src->n_ithr)
    return (1);

  if (strlen(targ->ithr_sign) == 0)
    CHTsetSignString(targ, src->ithr_sign);
  else if (strcmp(targ->ithr_sign, src->ithr_sign))
    return (1);

  if (targ->sthr_lo < 0)
    targ->sthr_lo = src->sthr_lo;
  else if (targ->sthr_lo != src->sthr_lo)
    return (1);

  if (targ->sthr_hi < 0)
    targ->sthr_hi = src->sthr_hi;
  else if (targ->sthr_hi != src->sthr_hi)
    return (1);

  if (targ->n_sthr < 0)
    targ->n_sthr = src->n_sthr;
  else if (targ->n_sthr != src->n_sthr)
    return (1);

  return (0);
}

/*--------------------------------------------------------------
  CHTsetSignString() - sets the sign string. If the string is
  unrecognized, returns 1, otherwise 0. If the string is NULL,
  abs is used.
  --------------------------------------------------------------*/
int CHTsetSignString(CHT *cht, char *ithr_sign)
{
  int ithr_signid;

  ithr_signid = CHTsignId(ithr_sign);
  if (ithr_signid == -100) return (1);

  cht->ithr_signid = ithr_signid;

  if (ithr_sign == NULL)
    strcpy(cht->ithr_sign, "abs");
  else
    strcpy(cht->ithr_sign, ithr_sign);

  return (0);
}

/*--------------------------------------------------------------
  CHTsignId() - converts the sign string to a numeric code. The
  code is set in the CHT structure and returned.
  --------------------------------------------------------------*/
int CHTsignId(char *ithr_sign)
{
  if (ithr_sign == NULL) return (0);  // abs
  if (!strcasecmp(ithr_sign, "pos"))
    return (+1);
  else if (!strcasecmp(ithr_sign, "abs"))
    return (0);
  else if (!strcasecmp(ithr_sign, "neg"))
    return (-1);

  printf("ERROR: ithr_sign = %s, unrecognized.\n", ithr_sign);
  printf("       ithr_sign must be pos, neg, or abs.\n");
  return (-100);
}

/*-------###########################################--------*/
/*-------###########################################--------*/
/*-------###########################################--------*/

/*------------------------------------------------------------
  CSDalloc() - just allocs structure and sets some initial vars
  ------------------------------------------------------------*/
CLUSTER_SIM_DATA *CSDalloc(void)
{
  CLUSTER_SIM_DATA *csd;

  csd = (CLUSTER_SIM_DATA *)calloc(sizeof(CLUSTER_SIM_DATA), 1);
  csd->mergedflag = 0;  // not a merged data
  memset(csd->simtype, 0, strlen(csd->simtype));
  memset(csd->anattype, 0, strlen(csd->anattype));
  memset(csd->subject, 0, strlen(csd->subject));
  memset(csd->hemi, 0, strlen(csd->hemi));
  memset(csd->contrast, 0, strlen(csd->contrast));
  csd->seed = -1;
  csd->thresh = -1;
  csd->threshsign = 0;
  csd->nullfwhm = -1;
  csd->varfwhm = -1;
  csd->searchspace = -1;
  csd->nreps = -1;
  // Always do this now (4/9/10)
  csd->FixGroupSubjectArea = 1;
  // if(getenv("FIX_VERTEX_AREA") == NULL) csd->FixGroupSubjectArea = 0;
  // else                                  csd->FixGroupSubjectArea = 1;

  return (csd);
}

/*--------------------------------------------------------------
  CSDread() - reads a cluster simulation data file. The format
  of this file is currently defined by mri_glmfit.
  --------------------------------------------------------------*/
CLUSTER_SIM_DATA *CSDread(char *csdfile)
{
  FILE *fp;
  CLUSTER_SIM_DATA *csd;
  char tag[1000], tmpstr[1000], c;
  int r, nthrep, nrepstmp;
  double d;

  fp = fopen(csdfile, "r");
  if (fp == NULL) {
    printf("ERROR: CSDread(): could not open %s\n", csdfile);
    return (NULL);
  }

  csd = (CLUSTER_SIM_DATA *)calloc(sizeof(CLUSTER_SIM_DATA), 1);
  csd->mergedflag = 0;  // not a merged data

  // Go through each input line
  nthrep = 0;
  while (1) {
    r = fscanf(fp, "%s", tag);
    if (r == EOF) break;

    if (!strcmp(tag, "#")) {
      // ----------- Header ------------
      fscanf(fp, "%s", tag);
      if (!strcmp(tag, "simtype"))
        fscanf(fp, "%s", csd->simtype);
      else if (!strcmp(tag, "anattype")) {
        fscanf(fp, "%s", csd->anattype);
        if (!strcmp(csd->anattype, "surface")) {
          fscanf(fp, "%s", csd->subject);
          fscanf(fp, "%s", csd->hemi);
        }
      }
      else if (!strcmp(tag, "thresh"))
        fscanf(fp, "%lf", &(csd->thresh));
      else if (!strcmp(tag, "threshsign"))
        fscanf(fp, "%lf", &(csd->threshsign));
      else if (!strcmp(tag, "seed"))
        fscanf(fp, "%ld", &(csd->seed));
      else if (!strcmp(tag, "contrast"))
        fscanf(fp, "%s", csd->contrast);
      else if (!strcmp(tag, "nullfwhm"))
        fscanf(fp, "%lf", &csd->nullfwhm);
      else if (!strcmp(tag, "varfwhm"))
        fscanf(fp, "%lf", &csd->varfwhm);
      else if (!strcmp(tag, "searchspace"))
        fscanf(fp, "%lf", &csd->searchspace);
      else if (!strcmp(tag, "nreps")) {
        fscanf(fp, "%d", &(nrepstmp));
        if (!strcmp(csd->anattype, "surface") && nrepstmp != -1) {
          printf("ERROR: this CSD file (%s) is an older version and\n", csdfile);
          printf("is incompatible with post-Dec08 versions of FreeSurfer\n");
          // Post-Dec08 versions set nreps=-1 and save the "nrepetitions"
          // tag to store the number of repetitions for surface data.
          exit(1);
        }
        if (!strcmp(csd->anattype, "volume")) {
          csd->nreps = nrepstmp;
          CSDallocData(csd);
        }
      }
      else if (!strcmp(tag, "nrepetitions")) {
        fscanf(fp, "%d", &(csd->nreps));
        if (!strcmp(csd->anattype, "surface")) CSDallocData(csd);
      }
      else if (!strcmp(tag, "FixGroupSubjectArea"))
        fscanf(fp, "%d", &(csd->FixGroupSubjectArea));
      else
        fgets(tmpstr, 1000, fp);  // not an interesting line, so get past it
    }
    else {
      // ----------- Data ------------
      // printf("%s \n",tag);
      fscanf(fp, "%d", &(csd->nClusters[nthrep]));
      fscanf(fp, "%lf", &(csd->MaxClusterSize[nthrep]));
      fscanf(fp, "%lf", &(csd->MaxSig[nthrep]));
      // Check for max stat on this line
      c = fgetc(fp);
      while (c == ' ') c = fgetc(fp);
      ungetc(c, fp);
      if (c != '\n') {
        fscanf(fp, "%lf", &d);
        csd->MaxStat[nthrep] = d;
        // printf("d = %g\n",d);
      }
      // exit(1);
      nthrep++;
    }
  }
  return (csd);
}

/*--------------------------------------------------------------
  CSDreadMerge() - reads in a CSD file and merges it with
  another CSD. If the input csd is NULL, then it is the
  same as CSDread(). The purpose is to be able to do somthing
  like this: csd = CSDreadMerge(csdfile, csd);
  --------------------------------------------------------------*/
CLUSTER_SIM_DATA *CSDreadMerge(char *csdfile, CSD *csd)
{
  CLUSTER_SIM_DATA *csd2 = NULL, *csdmerged = NULL;

  if (csd == NULL) {
    csd = CSDread(csdfile);
    return (csd);
  }

  csd2 = CSDread(csdfile);
  if (csd2 == NULL) return (NULL);

  csdmerged = CSDmerge(csd, csd2);
  if (csdmerged == NULL) {
    CSDfreeData(csd2);
    return (NULL);
  }

  CSDcopy(csdmerged, csd);

  CSDfreeData(csd2);
  CSDfreeData(csdmerged);
  return (csd);
}

/*--------------------------------------------------------------
  CSDallocData() - allocates the arrays for a CSD (not the
  structure iteself).
  --------------------------------------------------------------*/
int CSDallocData(CLUSTER_SIM_DATA *csd)
{
  csd->nClusters = (int *)calloc(csd->nreps, sizeof(int));
  csd->MaxClusterSize = (double *)calloc(csd->nreps, sizeof(double));
  csd->MaxClusterSizeVtx = (double *)calloc(csd->nreps, sizeof(double));
  csd->MaxClusterWeightVtx = (double *)calloc(csd->nreps, sizeof(double));
  csd->MaxClusterWeightArea = (double *)calloc(csd->nreps, sizeof(double));
  csd->MaxSig = (double *)calloc(csd->nreps, sizeof(double));
  csd->MaxStat = (double *)calloc(csd->nreps, sizeof(double));
  return (0);
}

/*--------------------------------------------------------------
  CSDfreeData() - frees the data arrays and sets their pointers
  to NULL. Does not try to free the structure.
  --------------------------------------------------------------*/
int CSDfreeData(CLUSTER_SIM_DATA *csd)
{
  if (csd->nClusters) {
    free(csd->nClusters);
    csd->nClusters = NULL;
  }
  if (csd->MaxClusterSize) {
    free(csd->MaxClusterSize);
    csd->MaxClusterSize = NULL;
  }
  if (csd->MaxClusterSizeVtx) {
    free(csd->MaxClusterSizeVtx);
    csd->MaxClusterSizeVtx = NULL;
  }
  if (csd->MaxClusterWeightVtx) {
    free(csd->MaxClusterWeightVtx);
    csd->MaxClusterWeightVtx = NULL;
  }
  if (csd->MaxClusterWeightArea) {
    free(csd->MaxClusterWeightArea);
    csd->MaxClusterWeightArea = NULL;
  }
  if (csd->MaxSig) {
    free(csd->MaxSig);
    csd->MaxSig = NULL;
  }
  if (csd->MaxStat) {
    free(csd->MaxStat);
    csd->MaxStat = NULL;
  }
  if (csd->mcs_pdf) HISTOfree(&csd->mcs_pdf);
  if (csd->mcs_cdf) HISTOfree(&csd->mcs_cdf);
  if (csd->ms_pdf) HISTOfree(&csd->ms_pdf);
  if (csd->ms_cdf) HISTOfree(&csd->ms_cdf);
  if (csd->grf_cdf) free(csd->grf_cdf);
  csd->grf_cdf = NULL;

  return (0);
}

/*--------------------------------------------------------------
  CSDcopy() - copies csd into csdcopy. If csdcopy is NULL, it is
  allocated. If csdcopy is non-null, the data are freed and
  then re-allocated.
  --------------------------------------------------------------*/
CSD *CSDcopy(CSD *csd, CSD *csdcopy)
{
  int nthrep;

  if (csdcopy == NULL)
    csdcopy = (CLUSTER_SIM_DATA *)calloc(sizeof(CLUSTER_SIM_DATA), 1);
  else
    CSDfreeData(csdcopy);

  strcpy(csdcopy->simtype, csd->simtype);
  strcpy(csdcopy->anattype, csd->anattype);
  strcpy(csdcopy->subject, csd->subject);
  strcpy(csdcopy->hemi, csd->hemi);
  csdcopy->thresh = csd->thresh;
  csdcopy->threshsign = csd->threshsign;
  csdcopy->searchspace = csd->searchspace;
  csdcopy->nullfwhm = csd->nullfwhm;
  csdcopy->varfwhm = csd->varfwhm;
  csdcopy->FixGroupSubjectArea = csd->FixGroupSubjectArea;

  csdcopy->nreps = csd->nreps;
  CSDallocData(csdcopy);

  for (nthrep = 0; nthrep < csd->nreps; nthrep++) {
    csdcopy->nClusters[nthrep] = csd->nClusters[nthrep];
    csdcopy->MaxClusterSize[nthrep] = csd->MaxClusterSize[nthrep];
    csdcopy->MaxClusterSizeVtx[nthrep] = csd->MaxClusterSizeVtx[nthrep];
    csdcopy->MaxClusterWeightVtx[nthrep] = csd->MaxClusterWeightVtx[nthrep];
    csdcopy->MaxClusterWeightArea[nthrep] = csd->MaxClusterWeightArea[nthrep];
    csdcopy->MaxSig[nthrep] = csd->MaxSig[nthrep];
    csdcopy->MaxStat[nthrep] = csd->MaxStat[nthrep];
  }
  return (csdcopy);
}

/*--------------------------------------------------------------
  CSDmerge() - merge two CSDs into one. Requires that:
  (1) simtypes be the same
  (2) anattypes be the same
  (3) contrasts be the same
  (4) thresholds be the same
  (5) threshold signs be the same
  (6) nullfwhm be the same
  (7) varfwhm be the same
  (8) searchspace be the same
  (9) seeds be different
  (10) FixGroupSubjectArea be the same (for surfaces)
  The seed from the first is copied into the merge, and the
  mergeflag is set to 1.
  --------------------------------------------------------------*/
CSD *CSDmerge(CSD *csd1, CSD *csd2)
{
  int nthrep1, nthrep2, nthrep;
  CSD *csd;

  if (strcmp(csd1->simtype, csd2->simtype)) {
    printf("ERROR: CSDmerge: CSDs have different sim types\n");
    return (NULL);
  }
  if (strcmp(csd1->anattype, csd2->anattype)) {
    printf("ERROR: CSDmerge: CSDs have different anat types\n");
    return (NULL);
  }
  if (!strcmp(csd1->anattype, "surface")) {
    if (csd1->FixGroupSubjectArea != csd2->FixGroupSubjectArea) {
      printf("ERROR: CSDmerge: CSDs have different FixGroupSubjectArea flags\n");
      return (NULL);
    }
  }
  if (strcmp(csd1->contrast, csd2->contrast)) {
    printf("ERROR: CSDmerge: CSDs have different contrasts\n");
    return (NULL);
  }
  if (csd1->thresh != csd2->thresh) {
    printf("ERROR: CSDmerge: CSDs have different thresholds\n");
    return (NULL);
  }
  if (csd1->threshsign != csd2->threshsign) {
    printf("ERROR: CSDmerge: CSDs have different threshold signs\n");
    return (NULL);
  }
  if (fabs(csd1->searchspace - csd2->searchspace) > .01) {
    // These can sometimes be different when run on different machines
    printf("ERROR: CSDmerge: CSDs have different search spaces\n");
    printf("%g %g   %g\n,", csd1->searchspace, csd2->searchspace, csd1->searchspace - csd2->searchspace);
    return (NULL);
  }
  if (csd1->nullfwhm != csd2->nullfwhm) {
    printf("ERROR: CSDmerge: CSDs have different null fwhm\n");
    return (NULL);
  }
  if (csd1->varfwhm != csd2->varfwhm) {
    printf("ERROR: CSDmerge: CSDs have different variance fwhm\n");
    return (NULL);
  }
  if (csd1->seed == csd2->seed) {
    printf("ERROR: CSDmerge: CSDs have same seed\n");
    return (NULL);
  }

  csd = (CLUSTER_SIM_DATA *)calloc(sizeof(CLUSTER_SIM_DATA), 1);
  strcpy(csd->simtype, csd1->simtype);
  strcpy(csd->anattype, csd1->anattype);
  strcpy(csd->contrast, csd1->contrast);
  strcpy(csd->subject, csd1->subject);
  strcpy(csd->hemi, csd1->hemi);
  csd->thresh = csd1->thresh;
  csd->threshsign = csd1->threshsign;
  csd->searchspace = csd1->searchspace;
  csd->nullfwhm = csd1->nullfwhm;
  csd->varfwhm = csd1->varfwhm;
  csd->seed = csd1->seed;
  csd->mergedflag = 1;
  csd->FixGroupSubjectArea = csd1->FixGroupSubjectArea;

  csd->nreps = csd1->nreps + csd2->nreps;
  CSDallocData(csd);

  nthrep = 0;
  for (nthrep1 = 0; nthrep1 < csd1->nreps; nthrep1++) {
    csd->nClusters[nthrep] = csd1->nClusters[nthrep1];
    csd->MaxClusterSize[nthrep] = csd1->MaxClusterSize[nthrep1];
    csd->MaxClusterSizeVtx[nthrep] = csd1->MaxClusterSizeVtx[nthrep1];
    csd->MaxClusterWeightVtx[nthrep] = csd1->MaxClusterWeightVtx[nthrep1];
    csd->MaxClusterWeightArea[nthrep] = csd1->MaxClusterWeightArea[nthrep1];
    csd->MaxSig[nthrep] = csd1->MaxSig[nthrep1];
    csd->MaxStat[nthrep] = csd1->MaxStat[nthrep1];
    nthrep++;
  }
  for (nthrep2 = 0; nthrep2 < csd2->nreps; nthrep2++) {
    csd->nClusters[nthrep] = csd2->nClusters[nthrep2];
    csd->MaxClusterSize[nthrep] = csd2->MaxClusterSize[nthrep2];
    csd->MaxClusterSizeVtx[nthrep] = csd2->MaxClusterSizeVtx[nthrep1];
    csd->MaxClusterWeightVtx[nthrep] = csd2->MaxClusterWeightVtx[nthrep1];
    csd->MaxClusterWeightArea[nthrep] = csd2->MaxClusterWeightArea[nthrep1];
    csd->MaxSig[nthrep] = csd2->MaxSig[nthrep2];
    csd->MaxStat[nthrep] = csd2->MaxStat[nthrep2];
    nthrep++;
  }

  return (csd);
}

/*--------------------------------------------------------------
  CSDprintHeader() - print CSD header to the given stream.
  --------------------------------------------------------------*/
int CSDprintHeader(FILE *fp, CLUSTER_SIM_DATA *csd)
{
  fprintf(fp, "# simtype %s\n", csd->simtype);
  if (!strcmp(csd->anattype, "surface")) {
    fprintf(fp, "# anattype %s  %s %s\n", csd->anattype, csd->subject, csd->hemi);
    fprintf(fp, "# FixGroupSubjectArea %d\n", csd->FixGroupSubjectArea);
  }
  else
    fprintf(fp, "# anattype %s \n", csd->anattype);
  fprintf(fp, "# merged      %d\n", csd->mergedflag);
  fprintf(fp, "# contrast    %s\n", csd->contrast);
  fprintf(fp, "# seed        %ld\n", csd->seed);
  fprintf(fp, "# thresh      %lf\n", csd->thresh);
  fprintf(fp, "# threshsign  %lf\n", csd->threshsign);
  fprintf(fp, "# searchspace %lf\n", csd->searchspace);
  fprintf(fp, "# nullfwhm    %lf\n", csd->nullfwhm);
  fprintf(fp, "# varfwhm     %lf\n", csd->varfwhm);

  // Both volume and surface will have nrepetitions, and it will
  // always be accurate.
  fprintf(fp, "# nrepetitions %d\n", csd->nreps);
  fprintf(fp, "# NOTE: nreps and nrepetitions are both valid for volume data.\n");
  fprintf(fp, "# NOTE: nreps is invalid (-1) for surface data to assure.\n");
  fprintf(fp, "# NOTE:   backwards INcompatibility.\n");
  fprintf(fp, "# %s\n", getVersion().c_str());

  if (!strcmp(csd->anattype, "surface")) {
    // In Dec 2008, a bug was found in sclustSurfaceArea() which
    // caused the cluster surface area of group surfaces to be 20-25%
    // too large. That function is called by both mri_glmfit and
    // mri_surfcluster. The modifications below were made to assure
    // that CSD files generated by the new/fixed mri_glmfit are not
    // used by an old/unfixed version of mri_surfcluster. This
    // is done by setting nreps to 0, which will cause old versions
    // of the reader to die. New CSDs will instead have a tag
    // called nrepetitions which the new reader will use instead
    // of nreps.
    // fprintf(fp,"# nreps       %d\n",csd->nreps);
    fprintf(fp, "# nreps       -1\n");
    // This is just an explicit indication that it has been fixed.
    // Also assures that the new surfcluster.c is linked to.
    fprintf(fp, "# FixSurfClusterArea %d\n", FixSurfClusterArea);
    // Version of surfcluster.c
    fprintf(fp, "# %s\n", getVersion().c_str());
  }
  else {
    // Volume
    // Use the original for volume search space
    fprintf(fp, "# nreps       %d\n", csd->nreps);
  }

  return (0);
}

/*--------------------------------------------------------------
  CSDprint() - prints a CSD to the given stream.
  --------------------------------------------------------------*/
int CSDprint(FILE *fp, CLUSTER_SIM_DATA *csd)
{
  int nthrep;
  CSDprintHeader(fp, csd);
  fprintf(fp, "# LoopNo nClusters MaxClustSize        MaxSig    MaxStat\n");
  for (nthrep = 0; nthrep < csd->nreps; nthrep++) {
    fprintf(fp,
            "%7d       %3d      %g          %g     %g\n",
            nthrep,
            csd->nClusters[nthrep],
            csd->MaxClusterSize[nthrep],
            csd->MaxSig[nthrep],
            csd->MaxStat[nthrep]);
  }
  return (0);
}
/*--------------------------------------------------------------
  CSDprintWeight() - prints a CSD to the given stream (with nvtx and weight)
  --------------------------------------------------------------*/
int CSDprintWeight(FILE *fp, CLUSTER_SIM_DATA *csd)
{
  int nthrep;
  CSDprintHeader(fp, csd);
  fprintf(fp, "# LoopNo nClusters MaxClustSize         MaxSig      MaxStat ");
  fprintf(fp, "MaxClustSizeVtx    MaxClustWeightVtx\n");
  for (nthrep = 0; nthrep < csd->nreps; nthrep++) {
    fprintf(fp,
            "%7d       %3d      %g          %g     %g         %d           %g\n",
            nthrep,
            csd->nClusters[nthrep],
            csd->MaxClusterSize[nthrep],
            csd->MaxSig[nthrep],
            csd->MaxStat[nthrep],
            (int)csd->MaxClusterSizeVtx[nthrep],
            csd->MaxClusterWeightVtx[nthrep]);
  }
  return (0);
}

/*--------------------------------------------------------------
  CSDpvalClustSize() - computes the emperical pvalue for a given
  cluster size, as well as the confidence interval. ciPct is the conf
  interval given in percent (eg, 90 for 90%). This means that there
  will be a 90% chance that the "true" pvalue for the cluster will lie
  between pvalLow and pvalHi (based on binomial distribution). If
  no item from the simulation is larger than ClusterSize, then
  it is assumed that 1 item is so that things dont break.
  --------------------------------------------------------------*/
double CSDpvalClustSize(CLUSTER_SIM_DATA *csd, double ClusterSize, double ciPct, double *pvalLow, double *pvalHi)
{
  int nthrep, nover, k, nlow, nhi;
  double pval, psum, pcilow, pcihi;

  // First, count the number of MaxClusters whose size is greater than
  // the one under test
  nover = 0;
  for (nthrep = 0; nthrep < csd->nreps; nthrep++){
    if(csd->MaxClusterSize[nthrep] >  ClusterSize) nover++;
  }

  // Use greter-than-or-equal-to in addition to greather-than.  This
  // came up in a 1D application where the number of voxels in the
  // cluster can be quite small making the inequality important. Using
  // just > is too liberal. Using just >= is to conservative. This
  // effectively averages them together; seems to work. Probably not
  // important for 2D and 3D apps.
  char *UseGTE = getenv("FS_CSDPVALCLUSTSIZE_GTE");
  if(UseGTE != NULL){
    for (nthrep = 0; nthrep < csd->nreps; nthrep++){
      if(csd->MaxClusterSize[nthrep] >= ClusterSize) nover++;
    }
    nover /= 2.0;
  }

  // If none is over, then set nover = 1 so that things don't break
  if (nover == 0) nover = 1;

  // Compute the nomial pvalue
  pval = (double)nover / csd->nreps;

  // Ranges for confidence interval
  pcihi = ciPct / 100;
  pcilow = 1 - pcihi;

  // Loop thru all possible outcomes (k), computing the CDF of the
  // binomial distribution, until the upper conf interval is reached.
  nlow = -1;
  k = 0;
  psum = sc_ran_binomial_pdf(k, pval, csd->nreps);
  while (psum < pcihi && k <= csd->nreps) {
    // printf("%3d %lf\n",k,psum);
    if (nlow < 0 && psum > pcilow) nlow = k;
    k++;
    psum += sc_ran_binomial_pdf(k, pval, csd->nreps);
  }
  nhi = k;
  // Compute the pvalues at the lower and upper confidence intervals
  *pvalLow = (double)nlow / csd->nreps;
  *pvalHi = (double)nhi / csd->nreps;

  // printf("csize=%lf  p=%lf  ci=%lf  nLow=%d pLow=%lf nHi=%d pHi=%lf\n",
  //   ClusterSize,pval,ciPct,nlow,*pvalLow,nhi,*pvalHi);

  return (pval);
}

/*
  --------------------------------------------------------------------------
  CSDpvalMaxSig() - computes the emperical probability of getting a MaxSig
  greater than the given value (taking into account the sign).
  --------------------------------------------------------------------------
*/
double CSDpvalMaxSig(double val, CSD *csd)
{
  int n, nover;
  double pval;

  nover = 0;  // number of times maxsig exceeds the given value
  for (n = 0; n < csd->nreps; n++) {
    if (csd->threshsign == 0 && (fabs(csd->MaxSig[n]) > fabs(val)))
      nover++;
    else if (csd->threshsign > +0.5 && (csd->MaxSig[n] > val))
      nover++;
    else if (csd->threshsign < -0.5 && (csd->MaxSig[n] < val))
      nover++;
  }
  pval = (double)nover / csd->nreps;
  return (pval);
}

/*-------------------------------------------------------------------
  CSDpvalMaxSigMap() - computes the voxel-wise sig value of each voxel
  based on the CSD. The input and output are -log10(p)*sign, where
  sign is the sign of the input value. Bonf is for an additional
  Bonferroni correction (eg, 2 for across hemisphere or 3 for across
  hemis and subcortical).
  ------------------------------------------------------------------*/
MRI *CSDpvalMaxSigMap(MRI *sig, CSD *csd, MRI *mask, MRI *vwsig, double *maxmaxsig, int Bonf)
{
  int c, r, s, f, nhits, nvox;
  double m, val, voxsig, pval, maxvoxsig;

  if (vwsig == NULL) vwsig = MRIclone(sig, NULL);

  nvox = 0;
  nhits = 0;
  maxvoxsig = 0;
  for (s = 0; s < sig->depth; s++) {
    for (r = 0; r < sig->height; r++) {
      for (c = 0; c < sig->width; c++) {
        if (mask) {
          m = MRIgetVoxVal(mask, c, r, s, 0);
          if (m < 0.5) continue;
        }
        nvox++;
        for (f = 0; f < sig->nframes; f++) {
          val = MRIgetVoxVal(sig, c, r, s, f);
          if (fabs(val) > 0.0) {
            pval = CSDpvalMaxSig(val, csd);
            if (Bonf > 0) pval = 1 - pow((1 - pval), Bonf);
            voxsig = -SIGN(val) * log10(pval);
          }
          else
            voxsig = 0;
          if (fabs(voxsig) > 0) nhits++;
          MRIsetVoxVal(vwsig, c, r, s, f, voxsig);
          if (f == 0) {
            if (csd->threshsign == 0 && fabs(voxsig) > fabs(maxvoxsig)) maxvoxsig = voxsig;
            if (csd->threshsign > +0.5 && voxsig > maxvoxsig) maxvoxsig = voxsig;
            if (csd->threshsign < -0.5 && voxsig < maxvoxsig) maxvoxsig = voxsig;
            // printf("val=%g, voxsig %g max %g sign %g\n",val,voxsig,maxvoxsig,csd->threshsign);
          }
        }
      }
    }
  }
  printf("CSDpvalMaxSigMap(): found %d/%d above 0, max=%g\n", nhits, nvox, maxvoxsig);
  *maxmaxsig = maxvoxsig;
  return (vwsig);
}

/*-----------------------------------------------------------------------
  CSDcheckSimType() - checks simulation type string to make sure it
  is one that is recognized. Returns 0 if ok, 1 otherwise.
  -----------------------------------------------------------------------*/
int CSDcheckSimType(char *simtype)
{
  if (!strcmp(simtype, "perm")) return (0);
  if (!strcmp(simtype, "mc-full")) return (0);
  if (!strcmp(simtype, "mc-z")) return (0);
  if (!strcmp(simtype, "mc-t")) return (0);
  return (1);
}

/*-------------------------------------------------------------------
  CSDpdf() - computes pdf and cdf of Maximum Cluster Size and
  Maximum Sig. if nbins < 1, then nbins = sqrt(csd->nreps)
  -------------------------------------------------------------------*/
int CSDpdf(CSD *csd, int nbins)
{
  int n, dim;
  float min, max, ClusterSize;
  RFS *rfs;

  if (nbins < 1) {
    nbins = sqrt(csd->nreps);
    if (nbins == 0) nbins = 1;
  }

  if (!strcmp(csd->anattype, "surface"))
    dim = 2;  // surface
  else
    dim = 3;  // volume
  rfs = RFspecInit(0, NULL);
  rfs->name = strcpyalloc("gaussian");
  rfs->params[0] = 0;
  rfs->params[1] = 1;

  // Maximum Cluster Size ------------------------------------
  min = csd->MaxClusterSize[0];
  max = csd->MaxClusterSize[0];
  for (n = 0; n < csd->nreps; n++) {
    if (min > csd->MaxClusterSize[n]) min = csd->MaxClusterSize[n];
    if (max < csd->MaxClusterSize[n]) max = csd->MaxClusterSize[n];
  }
  csd->mcs_pdf = HISTObins(nbins, min, max);
  HISTOcount(csd->mcs_pdf, csd->MaxClusterSize, csd->nreps);

  csd->grf_cdf = (double *)calloc(nbins, sizeof(double));

  for (n = 0; n < csd->mcs_pdf->nbins; n++) {
    // printf("%d %g %g\n",n,csd->mcs_pdf->bins[n],csd->mcs_pdf->counts[n]);
    csd->mcs_pdf->counts[n] /= csd->nreps;
  }

  csd->mcs_cdf = HISTOcopy(csd->mcs_pdf, NULL);
  for (n = 1; n < csd->mcs_pdf->nbins; n++) {
    csd->mcs_cdf->counts[n] = csd->mcs_cdf->counts[n - 1] + csd->mcs_pdf->counts[n];
  }

  // Compute clusterwise sig with GRF. Do not subtract log10(2.0) one-sidedness.
  // The GRF will be bogus because needs to be in number of vertices/voxels
  for (n = 0; n < csd->mcs_pdf->nbins; n++) {
    ClusterSize = csd->mcs_pdf->bins[n];
    // csd->grf_cdf[n] =
    // RFprobZClusterSigThresh(ClusterSize, csd->thresh, csd->nullfwhm,
    //			      csd->searchspace, dim);
    csd->grf_cdf[n] = 1.0;  // bogus
  }

  // Maximum Sig ------------------------------------
  min = csd->MaxSig[0];
  max = csd->MaxSig[0];
  for (n = 0; n < csd->nreps; n++) {
    if (min > csd->MaxSig[n]) min = csd->MaxSig[n];
    if (max < csd->MaxSig[n]) max = csd->MaxSig[n];
  }
  csd->ms_pdf = HISTObins(nbins, min, max);
  HISTOcount(csd->ms_pdf, csd->MaxSig, csd->nreps);

  for (n = 0; n < csd->ms_pdf->nbins; n++) {
    // printf("%d %g %g\n",n,csd->ms_pdf->bins[n],csd->ms_pdf->counts[n]);
    csd->ms_pdf->counts[n] /= csd->nreps;
  }

  csd->ms_cdf = HISTOcopy(csd->ms_pdf, NULL);
  if (csd->threshsign >= 0) {
    for (n = 1; n < csd->ms_pdf->nbins; n++)
      csd->ms_cdf->counts[n] = csd->ms_cdf->counts[n - 1] + csd->ms_pdf->counts[n];
  }
  else {
    // For neg sign, integrate from pos to neg because the tail is neg
    for (n = csd->ms_pdf->nbins - 2; n >= 0; n--)
      csd->ms_cdf->counts[n] = csd->ms_cdf->counts[n + 1] + csd->ms_pdf->counts[n];
  }

  return (0);
}

/*------------------------------------------------------------------------*/
int CSDprintPDF(FILE *fp, CSD *csd)
{
  int nthbin;

  if (csd->mcs_pdf == NULL) {
    printf("ERROR: CSDprintPDF: csd pdf is NULL\n");
    return (1);
  }
  fprintf(fp, "# CSD PDF/CDF\n");
  CSDprintHeader(fp, csd);
  fprintf(fp, "# nbins %d\n", csd->mcs_pdf->nbins);
  fprintf(fp,
          "# BinNo  MaxClustBin MaxClustPDF MaxClustCDF GRFCDF"
          "    MaxSigBin MaxSigPDF MaxSigCDF\n");
  for (nthbin = 0; nthbin < csd->mcs_pdf->nbins; nthbin++) {
    fprintf(fp,
            "%4d      %f   %f     %f  %f     %f  %f  %f \n",
            nthbin,
            csd->mcs_pdf->bins[nthbin],
            csd->mcs_pdf->counts[nthbin],
            1.0 - csd->mcs_cdf->counts[nthbin],
            csd->grf_cdf[nthbin],
            csd->ms_pdf->bins[nthbin],
            csd->ms_pdf->counts[nthbin],
            1.0 - csd->ms_cdf->counts[nthbin]);
  }
  return (0);
}

/*------------------------------------------------------------------------*/
int CSDwritePDF(char *fname, CSD *csd)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if (!fp) {
    printf("ERROR: could not open %s\n", fname);
    return (1);
  }
  CSDprintPDF(fp, csd);
  return (0);
}
/*------------------------------------------------------------------------*/
int CSDwrite(char *fname, CSD *csd)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if (!fp) {
    printf("ERROR: could not open %s\n", fname);
    return (1);
  }
  CSDprint(fp, csd);
  return (0);
}
