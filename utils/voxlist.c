/**
 * @file  voxlist.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.21 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mri.h"
#include "diag.h"
#include "error.h"
#include "voxlist.h"
#include "macros.h"
#include "proto.h"

extern const char* Progname;

VOXEL_LIST  *
VLSTcreateInRegion(MRI *mri, float low_val, float hi_val ,
                   VOXEL_LIST *vl, int skip, int border_only,
                   MRI_REGION *box)
{
  int   x, y, z, nvox, i, width, height, depth ;
  double  val ;

  skip++ ;  /* next voxel + amount to skip */
  width = box->x+box->dx ;
  height = box->y+box->dy ;
  depth = box->z+box->dz ;
  for (nvox = 0, x = box->x ; x < width ; x+=skip)
  {
    for (y = box->y ; y < height ; y+=skip)
    {
      for (z = box->z ; z < depth ; z+=skip)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (val > 0)
          DiagBreak() ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (val >= low_val && val <= hi_val)
          nvox++ ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("allocating %d voxel indices...\n", nvox) ;
  if (vl == NULL)
    vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
  else if (vl->nvox < nvox)
  {
    free(vl->xi) ;
    free(vl->yi) ;
    free(vl->zi) ;
    vl->xi = vl->yi = vl->zi = NULL ;
  }
  if (vl->xi == NULL)
  {
    vl->xd = (float *)calloc(nvox, sizeof(float)) ;
    vl->yd = (float *)calloc(nvox, sizeof(float)) ;
    vl->zd = (float *)calloc(nvox, sizeof(float)) ;
    vl->vsrc = (float *)calloc(nvox, sizeof(float)) ;
    vl->vdst = (float *)calloc(nvox, sizeof(float)) ;
    vl->xi = (int *)calloc(nvox, sizeof(int)) ;
    vl->yi = (int *)calloc(nvox, sizeof(int)) ;
    vl->zi = (int *)calloc(nvox, sizeof(int)) ;
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
              Progname, nvox) ;
  vl->nvox = nvox ;
  for (nvox = 0, x = box->x ; x < width ; x+=skip)
  {
    for (y = box->y ; y < height ; y+=skip)
    {
      for (z = box->z ; z < depth ; z+=skip)
      {
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (val >= low_val && val <= hi_val)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          if ((border_only == 0) ||
              (MRIneighbors(mri, x, y, z, 0) > 0))
          {
            i = nvox++ ;
            vl->xi[i] = x ;
            vl->yi[i] = y ;
            vl->zi[i] = z ;
            vl->vsrc[i] = val ;
          }
        }
      }
    }
  }
  vl->mri = mri ;
  return(vl) ;
}
VOXEL_LIST *
VLSTcreate(MRI *mri,
           float low_val,
           float hi_val,
           VOXEL_LIST *vl,
           int skip,
           int border_only)
{
  int   x, y, z, nvox, i ;
  double  val ;

  skip++ ;  /* next voxel + amount to skip */
  for (nvox = x = 0 ; x < mri->width ; x+=skip)
  {
    for (y = 0 ; y < mri->height ; y+=skip)
    {
      for (z = 0 ; z < mri->depth ; z+=skip)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (val > 0)
          DiagBreak() ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (border_only && 
            (MRIneighborsInRange(mri, x, y, z, 0, low_val, hi_val) < 6))
          continue ;  // must be nbring at least one voxel not in [low_val, hi_val]
        if (val >= low_val && val <= hi_val)
          nvox++ ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("allocating %d voxel indices...\n", nvox) ;
  if (vl == NULL)
    vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
  else if (vl->nvox < nvox)
  {
    free(vl->xi) ;
    free(vl->yi) ;
    free(vl->zi) ;
    vl->xi = vl->yi = vl->zi = NULL ;
  }
  vl->nvox = nvox ;
  if (vl->xi == NULL)
  {
    if (nvox == 0)
      nvox = 1 ;
    vl->xd = (float *)calloc(nvox, sizeof(float)) ;
    vl->yd = (float *)calloc(nvox, sizeof(float)) ;
    vl->zd = (float *)calloc(nvox, sizeof(float)) ;
    vl->vsrc = (float *)calloc(nvox, sizeof(float)) ;
    vl->vdst = (float *)calloc(nvox, sizeof(float)) ;
    vl->xi = (int *)calloc(nvox, sizeof(int)) ;
    vl->yi = (int *)calloc(nvox, sizeof(int)) ;
    vl->zi = (int *)calloc(nvox, sizeof(int)) ;
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
              Progname, nvox) ;
  for (nvox = x = 0 ; x < mri->width ; x+=skip)
  {
    for (y = 0 ; y < mri->height ; y+=skip)
    {
      for (z = 0 ; z < mri->depth ; z+=skip)
      {
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (val >= low_val && val <= hi_val)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          if ((border_only == 0) ||
              (MRIneighborsInRange(mri, x, y, z, 0, low_val, hi_val) < 6))
          {
            i = nvox++ ;
            if (i == Gdiag_no)
              DiagBreak() ;
            if (x == Gx && y == Gy && z == Gz)
            {
              printf("voxel (%d, %d, %d) = %2.1f added to voxlist at %d\n",
                     x, y, z, val, i) ;
            }
            vl->xi[i] = x ;
            vl->yi[i] = y ;
            vl->zi[i] = z ;
            vl->vsrc[i] = val ;
          }
        }
      }
    }
  }
  vl->mri = mri ;
  return(vl) ;
}

typedef struct
{
  float vsrc, vdst, xd, yd, zd ;
  int   xi, yi, zi ;
} SORT_VOXEL ;
static int
compare_sort(const void *psv1, const void *psv2)
{
  SORT_VOXEL  *sv1, *sv2 ;

  sv1 = (SORT_VOXEL *)psv1 ;
  sv2 = (SORT_VOXEL *)psv2 ;

  if (sv1->vsrc < sv2->vsrc)
    return(1) ;
  else if (sv1->vsrc > sv2->vsrc)
    return(-1) ;

  return(0) ;
}

VOXEL_LIST *
VLSTsort(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst)
{
  SORT_VOXEL  *sort_voxels ;
  int         n ;

  if (vl_dst == NULL)
    vl_dst = VLSTcopy(vl_src, NULL, 0, vl_src->nvox) ;

  sort_voxels = (SORT_VOXEL *)calloc(vl_src->nvox, sizeof(SORT_VOXEL)) ;
  if (!sort_voxels)
    ErrorExit(ERROR_NOMEMORY,"MRIorderIndices: could not allocate sort table");

  
  for (n = 0 ; n < vl_dst->nvox ; n++)
  {
    sort_voxels[n].xi = vl_dst->xi[n] ;
    sort_voxels[n].yi = vl_dst->yi[n] ;
    sort_voxels[n].zi = vl_dst->zi[n] ;
    sort_voxels[n].vsrc = vl_dst->vsrc[n] ;
    sort_voxels[n].vdst = vl_dst->vdst[n] ;
    sort_voxels[n].xd = vl_dst->xd[n] ;
    sort_voxels[n].yd = vl_dst->yd[n] ;
    sort_voxels[n].zd = vl_dst->zd[n] ;
  }
  qsort(sort_voxels, vl_dst->nvox, sizeof(SORT_VOXEL), compare_sort) ;

  for (n = 0 ; n < vl_dst->nvox ; n++)
  {
    vl_dst->xi[n] = sort_voxels[n].xi  ;
    vl_dst->yi[n] = sort_voxels[n].yi  ;
    vl_dst->zi[n] = sort_voxels[n].zi  ;
    vl_dst->vsrc[n] = sort_voxels[n].vsrc  ;
    vl_dst->vdst[n] = sort_voxels[n].vdst  ;
    vl_dst->xd[n] = sort_voxels[n].xd  ;
    vl_dst->yd[n] = sort_voxels[n].yd  ;
    vl_dst->zd[n] = sort_voxels[n].zd  ;
  }

  free(sort_voxels) ;
  return(vl_dst) ;
}

int
VLSTfree(VOXEL_LIST **pvl)
{
  VOXEL_LIST *vl = *pvl ;
  *pvl = NULL ;

  if (!vl)
    return(ERROR_BADPARM) ;
  free(vl->xi) ;
  free(vl->yi) ;
  free(vl->zi) ;
  free(vl) ;
  return(NO_ERROR) ;
}
MRI *
VLSTcreateMri(VOXEL_LIST *vl, int val)
{
  int   i ;
  MRI *mri ;

  mri = MRIalloc(vl->mri->width, vl->mri->height, vl->mri->depth, MRI_UCHAR) ;
  MRIcopyHeader(vl->mri, mri) ;
  for (i = 0 ; i < vl->nvox ; i++)
  {
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
  }
  return(mri) ;
}

MRI *
VLSTaddToMri(VOXEL_LIST *vl, MRI *mri, int val)
{
  int   i ;

  if (mri == NULL)
  {
    mri = MRIalloc(vl->mri->width,
                   vl->mri->height,
                   vl->mri->depth,
                   MRI_UCHAR) ;
    MRIcopyHeader(vl->mri, mri) ;
  }
  for (i = 0 ; i < vl->nvox ; i++)
  {
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
  }
  return(mri) ;
}
MRI *
VLSTtoMri(VOXEL_LIST *vl, MRI *mri)
{
  int   i ;
  double  val ;

  if (mri == NULL)
    mri = MRIclone(vl->mri, NULL) ;

  for (i = 0 ; i < vl->nvox ; i++)
  {
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0) ;
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val) ;
  }
  return(mri) ;
}

VOXEL_LIST *
VLSTdilate(VOXEL_LIST *vl, int mode, MRI *mri_exclude)
{
  MRI         *mri_current, *mri_new ;
  int         i, xi, yi, zi, xk, yk, zk, nvox ;
  VOXEL_LIST  *vl_exp = NULL;

  // create volume to keep track of what voxels are in current list
  mri_current = VLSTcreateMri(vl, 1) ;
  mri_new = MRIclone(mri_current, NULL) ;

  // count how many vox will be in new list
  for (nvox = i = 0 ; i < vl->nvox ; i++)
  {
    for (xk = -1 ; xk <= 1 ; xk++)
    {
      xi = vl->mri->xi[vl->xi[i]+xk] ;
      for (yk = -1 ; yk <= 1 ; yk++)
      {
        yi = vl->mri->yi[vl->yi[i]+yk] ;
        for (zk = -1 ; zk <= 1 ; zk++)
        {
          zi = vl->mri->zi[vl->zi[i]+zk] ;
          if (nint(MRIgetVoxVal(mri_current, xi, yi, zi,0)) == 0 &&
              nint(MRIgetVoxVal(mri_new, xi, yi, zi,0)) == 0 &&
              (mri_exclude == NULL ||
               FZERO(MRIgetVoxVal(mri_exclude, xi, yi, zi, 0))))
          {
            MRIvox(mri_new, xi, yi, zi) = 1 ;
            nvox++ ;
          }
        }
      }
    }
  }

  if (nvox == 0)
    return(NULL) ;

  if (mode == VL_DILATE_ADD)
  {
    vl_exp = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
    vl_exp->nvox = vl->nvox + nvox ;
    vl_exp->mri = vl->mri ;
    vl_exp->xi = (int *)calloc(vl->nvox+nvox, sizeof(int)) ;
    vl_exp->yi = (int *)calloc(vl->nvox+nvox, sizeof(int)) ;
    vl_exp->zi = (int *)calloc(vl->nvox+nvox, sizeof(int)) ;
    if (!vl_exp || !vl_exp->xi || !vl_exp->yi || !vl_exp->zi)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
                Progname, nvox) ;
    for (i = 0 ; i < vl->nvox ; i++)
    {
      vl_exp->xi[i] = vl->xi[i] ;
      vl_exp->yi[i] = vl->yi[i] ;
      vl_exp->zi[i] = vl->zi[i] ;
    }
    nvox = vl->nvox ;
  }
  else if (mode == VL_DILATE_REPLACE)
  {
    vl_exp = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
    vl_exp->nvox = nvox ;
    vl_exp->mri = vl->mri ;
    vl_exp->xi = (int *)calloc(nvox, sizeof(int)) ;
    vl_exp->yi = (int *)calloc(nvox, sizeof(int)) ;
    vl_exp->zi = (int *)calloc(nvox, sizeof(int)) ;
    if (!vl_exp || !vl_exp->xi || !vl_exp->yi || !vl_exp->zi)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
                Progname, nvox) ;
    nvox = 0 ;
  }
  for (i = 0 ; i < vl->nvox ; i++)
  {
    for (xk = -1 ; xk <= 1 ; xk++)
    {
      xi = vl->mri->xi[vl->xi[i]+xk] ;
      for (yk = -1 ; yk <= 1 ; yk++)
      {
        yi = vl->mri->yi[vl->yi[i]+yk] ;
        for (zk = -1 ; zk <= 1 ; zk++)
        {
          zi = vl->mri->zi[vl->zi[i]+zk] ;
          if (MRIvox(mri_new, xi, yi, zi) == 1) // add it
          {
            MRIvox(mri_new, xi, yi, zi) = 0 ;  // only add it once
            vl_exp->xi[nvox] = xi ;
            vl_exp->yi[nvox] = yi ;
            vl_exp->zi[nvox] = zi ;
            nvox++ ;
          }
        }
      }
    }
  }

  MRIfree(&mri_current) ;
  MRIfree(&mri_new) ;

  return(vl_exp) ;
}



void VLSTcomputeStats(VOXEL_LIST *vl)
{
  double mean =0;
  double std = 0;
  int i;
  double val;

  if (vl->mri == NULL || vl->nvox <= 0)
  {
    vl->mean = 0;
    vl->std = 0;
    return;
  }

  for (i = 0 ; i < vl->nvox ; i++)
  {
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    mean += val;
    std += val*val;
  }

  mean /= (double)vl->nvox;
  std /= (double)vl->nvox;

  std = sqrt(std - mean*mean);

  vl->mean = mean;
  vl->std = std;

  return;
}
int
VLSTtransform(VOXEL_LIST *vl, MATRIX *m, MRI *mri, int sample_type)
{
  double   val, xd, yd, zd ;
  int    i ;

  VLSTtransformCoords(vl, m, 0) ;
  for (i = 0 ; i < vl->nvox ; i++)
  {
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0) ;
    vl->vsrc[i] = val ;
    xd = vl->xd[i] ; yd = vl->yd[i] ; zd = vl->zd[i] ;
    if (xd < 0)
      xd = 0 ;
    else if (xd >= mri->width-1)
      xd = mri->width-1 ;
    if (yd < 0)
      yd = 0 ;
    else if (yd >= mri->height-1)
      yd = mri->height-1 ;
    if (zd < 0)
      zd = 0 ;
    else if (zd >= mri->depth-1)
      zd = mri->depth-1 ;
    vl->xd[i] = xd ; vl->yd[i] = yd ; vl->zd[i] = zd ;
    MRIsampleVolumeType(mri, xd, yd, zd, &val, sample_type) ;
    vl->vdst[i] = val ;
  }
  return(NO_ERROR) ;
}
int
VLSTtransformCoords(VOXEL_LIST *vl, MATRIX *m, int skip)
{
  double   xd, yd, zd, val ;
  int    i, x, y, z ;
  static VECTOR *v1 = NULL, *v2 ;

  skip++ ;   // skip=0 means do every one
  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    v2 = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v1, 4, 1) = 1.0 ;
    *MATRIX_RELT(v2, 4, 1) = 1.0 ;
  }

  for (i = 0 ; i < vl->nvox ; i+=skip)
  {
    x = vl->xi[i] ; y = vl->yi[i] ; z = vl->zi[i] ;
    if (vl->mri)
    {
      val = MRIgetVoxVal(vl->mri, x, y, z, 0) ;
      vl->vsrc[i] = val ;
    }
    V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
    MatrixMultiply(m, v1, v2) ;
    xd = V3_X(v2) ; yd = V3_Y(v2) ; zd = V3_Z(v2) ;
    vl->xd[i] = xd ; vl->yd[i] = yd ; vl->zd[i] = zd ;
  }
  return(NO_ERROR) ;
}

VOXEL_LIST  *
VLSTcreateFromDifference(MRI *mri1, MRI *mri2, VOXEL_LIST *vl, int target_label)
{
  int   x, y, z, nvox, i ;
  double  val1, val2 ;

  for (nvox = x = 0 ; x < mri1->width ; x++)
  {
    for (y = 0 ; y < mri1->height ; y++)
    {
      for (z = 0 ; z < mri1->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val1 = MRIgetVoxVal(mri1, x, y, z, 0) ;
        val2 = MRIgetVoxVal(mri2, x, y, z, 0) ;
        if (val1 > 0)
          DiagBreak() ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if ((target_label >= 0) && (!FEQUAL(target_label,val1) && !FEQUAL(target_label,val2)))
          continue ;
        if (!FEQUAL(val1, val2))
          nvox++ ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("allocating %d voxel indices...\n", nvox) ;
  if (vl == NULL)
    vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
  else if (vl->nvox < nvox)
  {
    free(vl->xi) ;
    free(vl->yi) ;
    free(vl->zi) ;
    vl->xi = vl->yi = vl->zi = NULL ;
  }
  if (vl->xi == NULL)
  {
    vl->xd = (float *)calloc(nvox, sizeof(float)) ;
    vl->yd = (float *)calloc(nvox, sizeof(float)) ;
    vl->zd = (float *)calloc(nvox, sizeof(float)) ;
    vl->vsrc = (float *)calloc(nvox, sizeof(float)) ;
    vl->vdst = (float *)calloc(nvox, sizeof(float)) ;
    vl->xi = (int *)calloc(nvox, sizeof(int)) ;
    vl->yi = (int *)calloc(nvox, sizeof(int)) ;
    vl->zi = (int *)calloc(nvox, sizeof(int)) ;
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n",
              Progname, nvox) ;
  vl->nvox = nvox ;
  for (nvox = x = 0 ; x < mri1->width ; x++)
  {
    for (y = 0 ; y < mri1->height ; y++)
    {
      for (z = 0 ; z < mri1->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val1 = MRIgetVoxVal(mri1, x, y, z, 0) ;
        val2 = MRIgetVoxVal(mri2, x, y, z, 0) ;
        if (val1 > 0)
          DiagBreak() ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if ((target_label >= 0) && (!FEQUAL(target_label,val1) && !FEQUAL(target_label,val2)))
          continue ;
        if (!FEQUAL(val1, val2))
        {
          i = nvox++ ;
          vl->xi[i] = x ;
          vl->yi[i] = y ;
          vl->zi[i] = z ;
          vl->vsrc[i] = val1 ;
          vl->vdst[i] = val2 ;
        }
      }
    }
  }

  vl->mri = mri1 ;
  vl->mri2 = mri2 ;
  return(vl) ;
}
VOXEL_LIST *
VLSTalloc(int nvox)
{
  VOXEL_LIST *vl ;

  vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST)) ;
  if (vl == NULL)
    ErrorExit(ERROR_NOMEMORY, "VLSTalloc(%d): could not allocate VL struct", nvox) ;
  vl->xd = (float *)calloc(nvox, sizeof(float)) ;
  vl->yd = (float *)calloc(nvox, sizeof(float)) ;
  vl->zd = (float *)calloc(nvox, sizeof(float)) ;
  vl->vsrc = (float *)calloc(nvox, sizeof(float)) ;
  vl->vdst = (float *)calloc(nvox, sizeof(float)) ;
  vl->xi = (int *)calloc(nvox, sizeof(int)) ;
  vl->yi = (int *)calloc(nvox, sizeof(int)) ;
  vl->zi = (int *)calloc(nvox, sizeof(int)) ;
  vl->nvox = nvox ;
  if (vl->xi == NULL ||
      vl->yi == NULL ||
      vl->zi == NULL ||
      vl->xd == NULL ||
      vl->yd == NULL ||
      vl->zd == NULL ||
      vl->vsrc == NULL ||
      vl->vdst == NULL)
    ErrorExit(ERROR_NOMEMORY, "VLSTalloc(%d): could not allocate VL struct", nvox) ;
  return(vl) ;
}
VOXEL_LIST *
VLSTcopy(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst, int start_index, int num)
{
  int  i ;
  if (vl_dst == NULL)
    vl_dst = VLSTalloc(num) ;

  for (i = 0 ; i < num ; i++)
  {
    vl_dst->xi[i] = vl_src->xi[i+start_index] ;
    vl_dst->yi[i] = vl_src->yi[i+start_index] ;
    vl_dst->zi[i] = vl_src->zi[i+start_index] ;

    vl_dst->xd[i] = vl_src->xd[i+start_index] ;
    vl_dst->yd[i] = vl_src->yd[i+start_index] ;
    vl_dst->zd[i] = vl_src->zd[i+start_index] ;

    vl_dst->vsrc[i] = vl_src->vsrc[i+start_index] ;
    vl_dst->vdst[i] = vl_src->vdst[i+start_index] ;
  }

  return(vl_dst) ;
}
