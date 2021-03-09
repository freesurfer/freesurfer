/**
 * @brief utilities to create and process lists of voxels
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

#include "diag.h"
#include "error.h"
#include "mrisurf.h"
#include "macros.h"
#include "mri.h"
#include "proto.h"
#include "voxlist.h"

extern const char *Progname;

VOXEL_LIST *VLSTcreateInRegion(
    MRI *mri, float low_val, float hi_val, VOXEL_LIST *vl, int skip, int border_only, MRI_REGION *box)
{
  int x, y, z, f, nvox, i, width, height, depth;
  double val;

  skip++; /* next voxel + amount to skip */
  width = box->x + box->dx;
  height = box->y + box->dy;
  depth = box->z + box->dz;
  for (nvox = 0, x = box->x; x < width; x += skip) {
    for (y = box->y; y < height; y += skip) {
      for (z = box->z; z < depth; z += skip) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (val > 0) DiagBreak();
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        if (val >= low_val && val <= hi_val) nvox++;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("allocating %d voxel indices...\n", nvox);
  if (vl == NULL)
    vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST));
  else if (vl->max_vox < nvox) {
    free(vl->xi);
    free(vl->yi);
    free(vl->zi);
    free(vl->fi);
    vl->xi = vl->yi = vl->zi = vl->fi = NULL;
  }
  if (vl->xi == NULL) {
    vl->xd = (float *)calloc(nvox, sizeof(float));
    vl->yd = (float *)calloc(nvox, sizeof(float));
    vl->zd = (float *)calloc(nvox, sizeof(float));
    vl->vsrc = (float *)calloc(nvox, sizeof(float));
    vl->vdst = (float *)calloc(nvox, sizeof(float));
    vl->xi = (int *)calloc(nvox, sizeof(int));
    vl->yi = (int *)calloc(nvox, sizeof(int));
    vl->zi = (int *)calloc(nvox, sizeof(int));
    vl->fi = (int *)calloc(nvox, sizeof(int));
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi || !vl->fi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n", Progname, nvox);
  vl->nvox = nvox;

  for (nvox = 0, f = 0; f < mri->nframes; f++)
    for (x = box->x; x < width; x += skip) {
      for (y = box->y; y < height; y += skip) {
        for (z = box->z; z < depth; z += skip) {
          val = MRIgetVoxVal(mri, x, y, z, f);
          if (val >= low_val && val <= hi_val) {
            if (x == Gx && y == Gy && z == Gz) DiagBreak();
            if ((border_only == 0) || (MRIneighbors(mri, x, y, z, f) > 0)) {
              i = nvox++;
              vl->xi[i] = x;
              vl->yi[i] = y;
              vl->zi[i] = z;
              vl->fi[i] = f;
              vl->vsrc[i] = val;
            }
          }
        }
      }
    }
  vl->mri = mri;
  return (vl);
}
VOXEL_LIST *VLSTcreate(MRI *mri, float low_val, float hi_val, VOXEL_LIST *vl, int skip, int border_only)
{
  int x, y, z, nvox, i, f;
  double val;

  skip++; /* next voxel + amount to skip */
  for (nvox = f = 0; f < mri->nframes; f += skip)
    for (z = 0; z < mri->depth; z += skip) {
      for (y = 0; y < mri->height; y += skip) {
        for (x = 0; x < mri->width; x += skip) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          val = MRIgetVoxVal(mri, x, y, z, f);
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          if (border_only && (MRIneighborsInRange(mri, x, y, z, f, low_val, hi_val) < 6))
            continue;  // must be nbring at least one voxel not in [low_val, hi_val]
          if (val >= low_val && val <= hi_val) nvox++;
        }
      }
    }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("allocating %d voxel indices...\n", nvox);
  if (vl == NULL)
    vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST));
  else if (vl->max_vox < nvox) {
    free(vl->xi);
    free(vl->yi);
    free(vl->zi);
    free(vl->fi);
    vl->xi = vl->yi = vl->zi = vl->fi = NULL;
  }
  vl->nvox = nvox;
  if (vl->xi == NULL) {
    if (nvox == 0) nvox = 1;
    vl->xd = (float *)calloc(nvox, sizeof(float));
    vl->yd = (float *)calloc(nvox, sizeof(float));
    vl->zd = (float *)calloc(nvox, sizeof(float));
    vl->vsrc = (float *)calloc(nvox, sizeof(float));
    vl->vdst = (float *)calloc(nvox, sizeof(float));
    vl->xi = (int *)calloc(nvox, sizeof(int));
    vl->yi = (int *)calloc(nvox, sizeof(int));
    vl->zi = (int *)calloc(nvox, sizeof(int));
    vl->fi = (int *)calloc(nvox, sizeof(int));
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi || !vl->fi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n", Progname, nvox);
  for (nvox = f = 0; f < mri->nframes; f += skip)
    for (z = 0; z < mri->depth; z += skip) {
      for (y = 0; y < mri->height; y += skip) {
        for (x = 0; x < mri->width; x += skip) {
          val = MRIgetVoxVal(mri, x, y, z, f);
          if (val >= low_val && val <= hi_val) {
            if (x == Gx && y == Gy && z == Gz) DiagBreak();
            if ((border_only == 0) || (MRIneighborsInRange(mri, x, y, z, f, low_val, hi_val) < 6)) {
              i = nvox++;
              if (i == Gdiag_no) DiagBreak();
              if (x == Gx && y == Gy && z == Gz) {
                printf("voxel (%d, %d, %d) = %2.1f added to voxlist at %d\n", x, y, z, val, i);
              }
              vl->xi[i] = x;
              vl->yi[i] = y;
              vl->zi[i] = z;
              vl->fi[i] = f;
              vl->vsrc[i] = val;
            }
          }
        }
      }
    }
  vl->mri = mri;
  return (vl);
}

typedef struct
{
  float vsrc, vdst, xd, yd, zd;
  int xi, yi, zi;
} SORT_VOXEL;

static int compare_sort(const void *psv1, const void *psv2)
{
  SORT_VOXEL *sv1, *sv2;

  sv1 = (SORT_VOXEL *)psv1;
  sv2 = (SORT_VOXEL *)psv2;

  if (sv1->vsrc < sv2->vsrc) {
    return (1);
  }
  else if (sv1->vsrc > sv2->vsrc) {
    return (-1);
  }
  else {
    if (sv1->xi < sv2->xi) {
      return (1);
    }
    else if (sv1->xi > sv2->xi) {
      return (-1);
    }
    else {
      if (sv1->yi < sv2->yi) {
        return (1);
      }
      else if (sv1->yi > sv2->yi) {
        return (-1);
      }
      else {
        if (sv1->zi < sv2->zi) {
          return (1);
        }
        else if (sv1->zi > sv2->zi) {
          return (-1);
        }
      }
    }
  }

  return (0);
}

VOXEL_LIST *VLSTsort(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst)
{
  SORT_VOXEL *sort_voxels;
  int n;

  if (vl_dst == NULL) vl_dst = VLSTcopy(vl_src, NULL, 0, vl_src->nvox);

  sort_voxels = (SORT_VOXEL *)calloc(vl_src->nvox, sizeof(SORT_VOXEL));
  if (!sort_voxels) ErrorExit(ERROR_NOMEMORY, "MRIorderIndices: could not allocate sort table");

  for (n = 0; n < vl_dst->nvox; n++) {
    sort_voxels[n].xi = vl_dst->xi[n];
    sort_voxels[n].yi = vl_dst->yi[n];
    sort_voxels[n].zi = vl_dst->zi[n];
    sort_voxels[n].vsrc = vl_dst->vsrc[n];
    sort_voxels[n].vdst = vl_dst->vdst[n];
    sort_voxels[n].xd = vl_dst->xd[n];
    sort_voxels[n].yd = vl_dst->yd[n];
    sort_voxels[n].zd = vl_dst->zd[n];
  }
  qsort(sort_voxels, vl_dst->nvox, sizeof(SORT_VOXEL), compare_sort);

  for (n = 0; n < vl_dst->nvox; n++) {
    vl_dst->xi[n] = sort_voxels[n].xi;
    vl_dst->yi[n] = sort_voxels[n].yi;
    vl_dst->zi[n] = sort_voxels[n].zi;
    vl_dst->vsrc[n] = sort_voxels[n].vsrc;
    vl_dst->vdst[n] = sort_voxels[n].vdst;
    vl_dst->xd[n] = sort_voxels[n].xd;
    vl_dst->yd[n] = sort_voxels[n].yd;
    vl_dst->zd[n] = sort_voxels[n].zd;
  }

  free(sort_voxels);
  return (vl_dst);
}

int VLSTfree(VOXEL_LIST **pvl)
{
  VOXEL_LIST *vl = *pvl;
  *pvl = NULL;

  if (!vl) return (ERROR_BADPARM);
  free(vl->xi);
  free(vl->yi);
  free(vl->zi);
  free(vl->xd);
  free(vl->yd);
  free(vl->zd);

  free(vl->vsrc);
  free(vl->vdst);
  free(vl->fi);

  if (vl->t) free(vl->t);
  if (vl->mx) free(vl->mx);
  if (vl->my) free(vl->my);
  if (vl->mz) free(vl->mz);
  free(vl);
  return (NO_ERROR);
}
MRI *VLSTcreateMri(VOXEL_LIST *vl, int val)
{
  int i;
  MRI *mri;

  mri = MRIalloc(vl->mri->width, vl->mri->height, vl->mri->depth, MRI_UCHAR);
  MRIcopyHeader(vl->mri, mri);
  for (i = 0; i < vl->nvox; i++) {
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val);
  }
  return (mri);
}

MRI *VLSTaddToMri(VOXEL_LIST *vl, MRI *mri, int val)
{
  int i;

  if (mri == NULL) {
    mri = MRIalloc(vl->mri->width, vl->mri->height, vl->mri->depth, MRI_UCHAR);
    MRIcopyHeader(vl->mri, mri);
  }
  for (i = 0; i < vl->nvox; i++) {
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val);
  }
  return (mri);
}
MRI *VLSTtoMri(VOXEL_LIST *vl, MRI *mri)
{
  int i;
  double val;

  if (mri == NULL) {
    if (vl->mri == NULL) ErrorReturn(NULL, (ERROR_BADPARM, "VLSTtoMri: mri and vl->mri are both NULL"));
    mri = MRIclone(vl->mri, NULL);
  }

  for (i = 0; i < vl->nvox; i++) {
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val);
  }
  return (mri);
}

MRI *VLSTvsrcToMri(VOXEL_LIST *vl, MRI *mri)
{
  int i;
  double val;

  if (mri == NULL) mri = MRIclone(vl->mri, NULL);

  for (i = 0; i < vl->nvox; i++) {
    val = vl->vsrc[i];
    MRIsetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0, val);
  }
  return (mri);
}

VOXEL_LIST *VLSTdilate(VOXEL_LIST *vl, int mode, MRI *mri_exclude)
{
  MRI *mri_current, *mri_new;
  int i, xi, yi, zi, xk, yk, zk, nvox;
  VOXEL_LIST *vl_exp = NULL;

  // create volume to keep track of what voxels are in current list
  mri_current = VLSTcreateMri(vl, 1);
  mri_new = MRIclone(mri_current, NULL);

  // count how many vox will be in new list
  for (nvox = i = 0; i < vl->nvox; i++) {
    for (xk = -1; xk <= 1; xk++) {
      xi = vl->mri->xi[vl->xi[i] + xk];
      for (yk = -1; yk <= 1; yk++) {
        yi = vl->mri->yi[vl->yi[i] + yk];
        for (zk = -1; zk <= 1; zk++) {
          zi = vl->mri->zi[vl->zi[i] + zk];
          if (nint(MRIgetVoxVal(mri_current, xi, yi, zi, 0)) == 0 && nint(MRIgetVoxVal(mri_new, xi, yi, zi, 0)) == 0 &&
              (mri_exclude == NULL || FZERO(MRIgetVoxVal(mri_exclude, xi, yi, zi, 0)))) {
            MRIvox(mri_new, xi, yi, zi) = 1;
            nvox++;
          }
        }
      }
    }
  }

  if (nvox == 0) return (NULL);

  if (mode == VL_DILATE_ADD) {
    vl_exp = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST));
    vl_exp->nvox = vl->nvox + nvox;
    vl_exp->max_vox = vl->nvox + nvox;
    vl_exp->mri = vl->mri;
    vl_exp->xi = (int *)calloc(vl->nvox + nvox, sizeof(int));
    vl_exp->yi = (int *)calloc(vl->nvox + nvox, sizeof(int));
    vl_exp->zi = (int *)calloc(vl->nvox + nvox, sizeof(int));
    if (!vl_exp || !vl_exp->xi || !vl_exp->yi || !vl_exp->zi)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n", Progname, nvox);
    for (i = 0; i < vl->nvox; i++) {
      vl_exp->xi[i] = vl->xi[i];
      vl_exp->yi[i] = vl->yi[i];
      vl_exp->zi[i] = vl->zi[i];
    }
    nvox = vl->nvox;
  }
  else if (mode == VL_DILATE_REPLACE) {
    vl_exp = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST));
    vl_exp->nvox = vl_exp->max_vox = nvox;
    vl_exp->mri = vl->mri;
    vl_exp->xi = (int *)calloc(nvox, sizeof(int));
    vl_exp->yi = (int *)calloc(nvox, sizeof(int));
    vl_exp->zi = (int *)calloc(nvox, sizeof(int));
    if (!vl_exp || !vl_exp->xi || !vl_exp->yi || !vl_exp->zi)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n", Progname, nvox);
    nvox = 0;
  }
  for (i = 0; i < vl->nvox; i++) {
    for (xk = -1; xk <= 1; xk++) {
      xi = vl->mri->xi[vl->xi[i] + xk];
      for (yk = -1; yk <= 1; yk++) {
        yi = vl->mri->yi[vl->yi[i] + yk];
        for (zk = -1; zk <= 1; zk++) {
          zi = vl->mri->zi[vl->zi[i] + zk];
          if (MRIvox(mri_new, xi, yi, zi) == 1)  // add it
          {
            MRIvox(mri_new, xi, yi, zi) = 0;  // only add it once
            vl_exp->xi[nvox] = xi;
            vl_exp->yi[nvox] = yi;
            vl_exp->zi[nvox] = zi;
            nvox++;
          }
        }
      }
    }
  }

  MRIfree(&mri_current);
  MRIfree(&mri_new);

  return (vl_exp);
}

void VLSTcomputeStats(VOXEL_LIST *vl)
{
  double mean = 0;
  double std = 0;
  int i;
  double val;

  if (vl->mri == NULL || vl->nvox <= 0) {
    vl->mean = 0;
    vl->std = 0;
    return;
  }

  for (i = 0; i < vl->nvox; i++) {
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    mean += val;
    std += val * val;
  }

  mean /= (double)vl->nvox;
  std /= (double)vl->nvox;

  std = sqrt(std - mean * mean);

  vl->mean = mean;
  vl->std = std;

  return;
}
int VLSTtransform(VOXEL_LIST *vl, MATRIX *m, MRI *mri, int sample_type)
{
  double val, xd, yd, zd;
  int i;

  VLSTtransformCoords(vl, m, 0);
  for (i = 0; i < vl->nvox; i++) {
    val = MRIgetVoxVal(vl->mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    vl->vsrc[i] = val;
    xd = vl->xd[i];
    yd = vl->yd[i];
    zd = vl->zd[i];
    if (xd < 0)
      xd = 0;
    else if (xd >= mri->width - 1)
      xd = mri->width - 1;
    if (yd < 0)
      yd = 0;
    else if (yd >= mri->height - 1)
      yd = mri->height - 1;
    if (zd < 0)
      zd = 0;
    else if (zd >= mri->depth - 1)
      zd = mri->depth - 1;
    vl->xd[i] = xd;
    vl->yd[i] = yd;
    vl->zd[i] = zd;
    MRIsampleVolumeType(mri, xd, yd, zd, &val, sample_type);
    vl->vdst[i] = val;
  }
  return (NO_ERROR);
}
int VLSTtransformCoords(VOXEL_LIST *vl, MATRIX *m, int skip)
{
  double xd, yd, zd, val;
  int i, x, y, z;
  static VECTOR *v1 = NULL, *v2;

  skip++;  // skip=0 means do every one
  if (v1 == NULL) {
    v1 = VectorAlloc(4, MATRIX_REAL);
    v2 = VectorAlloc(4, MATRIX_REAL);
    *MATRIX_RELT(v1, 4, 1) = 1.0;
    *MATRIX_RELT(v2, 4, 1) = 1.0;
  }

  for (i = 0; i < vl->nvox; i += skip) {
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    if (vl->mri) {
      val = MRIgetVoxVal(vl->mri, x, y, z, 0);
      vl->vsrc[i] = val;
    }
    V3_X(v1) = x;
    V3_Y(v1) = y;
    V3_Z(v1) = z;
    MatrixMultiply(m, v1, v2);
    xd = V3_X(v2);
    yd = V3_Y(v2);
    zd = V3_Z(v2);
    vl->xd[i] = xd;
    vl->yd[i] = yd;
    vl->zd[i] = zd;
  }
  return (NO_ERROR);
}

VOXEL_LIST *VLSTcreateFromDifference(MRI *mri1, MRI *mri2, VOXEL_LIST *vl, int target_label)
{
  int x, y, z, f, nvox, i;
  double val1, val2;

  for (nvox = f = 0; f < mri1->nframes; f++)
    for (nvox = x = 0; x < mri1->width; x++) {
      for (y = 0; y < mri1->height; y++) {
        for (z = 0; z < mri1->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          val1 = MRIgetVoxVal(mri1, x, y, z, f);
          val2 = MRIgetVoxVal(mri2, x, y, z, f);
          if (val1 > 0) DiagBreak();
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          if ((target_label >= 0) && (!FEQUAL(target_label, val1) && !FEQUAL(target_label, val2))) continue;
          if (!FEQUAL(val1, val2)) nvox++;
        }
      }
    }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("allocating %d voxel indices...\n", nvox);
  if (vl == NULL)
    vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST));
  else if (vl->max_vox < nvox) {
    free(vl->fi);
    free(vl->xi);
    free(vl->yi);
    free(vl->zi);
    vl->xi = vl->yi = vl->zi = vl->fi = NULL;
  }
  if (vl->xi == NULL) {
    vl->xd = (float *)calloc(nvox, sizeof(float));
    vl->yd = (float *)calloc(nvox, sizeof(float));
    vl->zd = (float *)calloc(nvox, sizeof(float));
    vl->vsrc = (float *)calloc(nvox, sizeof(float));
    vl->vdst = (float *)calloc(nvox, sizeof(float));
    vl->xi = (int *)calloc(nvox, sizeof(int));
    vl->yi = (int *)calloc(nvox, sizeof(int));
    vl->zi = (int *)calloc(nvox, sizeof(int));
    vl->fi = (int *)calloc(nvox, sizeof(int));
    vl->max_vox = nvox;
  }
  if (!vl || !vl->xi || !vl->yi || !vl->zi)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d voxel list\n", Progname, nvox);
  vl->nvox = nvox;
  for (nvox = f = 0; f < mri1->nframes; f++)
    for (x = 0; x < mri1->width; x++) {
      for (y = 0; y < mri1->height; y++) {
        for (z = 0; z < mri1->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          val1 = MRIgetVoxVal(mri1, x, y, z, f);
          val2 = MRIgetVoxVal(mri2, x, y, z, f);
          if (val1 > 0) DiagBreak();
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          if ((target_label >= 0) && (!FEQUAL(target_label, val1) && !FEQUAL(target_label, val2))) continue;
          if (!FEQUAL(val1, val2)) {
            i = nvox++;
            vl->xi[i] = x;
            vl->yi[i] = y;
            vl->zi[i] = z;
            vl->fi[i] = f;
            vl->vsrc[i] = val1;
            vl->vdst[i] = val2;
          }
        }
      }
    }

  vl->mri = mri1;
  vl->mri2 = mri2;
  return (vl);
}
VOXEL_LIST *VLSTalloc(int nvox)
{
  VOXEL_LIST *vl;

  vl = (VOXEL_LIST *)calloc(1, sizeof(VOXEL_LIST));
  if (vl == NULL) ErrorExit(ERROR_NOMEMORY, "VLSTalloc(%d): could not allocate VL struct", nvox);
  vl->xd = (float *)calloc(nvox, sizeof(float));
  vl->yd = (float *)calloc(nvox, sizeof(float));
  vl->zd = (float *)calloc(nvox, sizeof(float));
  vl->vsrc = (float *)calloc(nvox, sizeof(float));
  vl->vdst = (float *)calloc(nvox, sizeof(float));
  vl->xi = (int *)calloc(nvox, sizeof(int));
  vl->yi = (int *)calloc(nvox, sizeof(int));
  vl->zi = (int *)calloc(nvox, sizeof(int));
  vl->fi = (int *)calloc(nvox, sizeof(int));
  vl->nvox = nvox;
  vl->max_vox = nvox;
  if (vl->xi == NULL || vl->yi == NULL || vl->zi == NULL || vl->fi == NULL || vl->xd == NULL || vl->yd == NULL ||
      vl->zd == NULL || vl->vsrc == NULL || vl->vdst == NULL)
    ErrorExit(ERROR_NOMEMORY, "VLSTalloc(%d): could not allocate VL struct", nvox);
  return (vl);
}
VOXEL_LIST *VLSTcopy(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst, int start_src_index, int num)
{
  int i;

  if (vl_dst == NULL) vl_dst = VLSTalloc(num);

  if (vl_dst->max_vox < num)
    ErrorExit(ERROR_NOMEMORY,
              "VLSTcopy: destination voxlist not big enough (%d, but must be %d)",
              vl_dst->max_vox,
              start_src_index + num);

  vl_dst->type = vl_src->type;

  if (vl_src->mx) vl_dst->mx = (float *)calloc(num, sizeof(float));
  if (vl_src->my) vl_dst->my = (float *)calloc(num, sizeof(float));
  if (vl_src->mz) vl_dst->mz = (float *)calloc(num, sizeof(float));
  if (vl_src->t) vl_dst->t = (float *)calloc(num, sizeof(float));

  if (vl_src->mx && (vl_dst->t == NULL || vl_dst->mx == NULL || vl_dst->my == NULL || vl_dst->mz == NULL))
    ErrorExit(ERROR_NOMEMORY, "VLSTcopy: could not allocate %d-len slope arrays", num);

  for (i = 0; i < num; i++) {
    vl_dst->xi[i] = vl_src->xi[i + start_src_index];
    vl_dst->yi[i] = vl_src->yi[i + start_src_index];
    vl_dst->zi[i] = vl_src->zi[i + start_src_index];
    vl_dst->fi[i] = vl_src->fi[i + start_src_index];

    vl_dst->xd[i] = vl_src->xd[i + start_src_index];
    vl_dst->yd[i] = vl_src->yd[i + start_src_index];
    vl_dst->zd[i] = vl_src->zd[i + start_src_index];

    vl_dst->vsrc[i] = vl_src->vsrc[i + start_src_index];
    vl_dst->vdst[i] = vl_src->vdst[i + start_src_index];
  }

  if (vl_dst->nvox < num) vl_dst->nvox = num;
  return (vl_dst);
}
VOXEL_LIST *VLSTcopyInto(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst, int start_dst_index, int num)
{
  int i;

  if (vl_dst == NULL) vl_dst = VLSTalloc(start_dst_index + num);

  if (vl_dst->max_vox < start_dst_index + num)
    ErrorExit(ERROR_NOMEMORY,
              "VLSTcopy: destination voxlist not big enough (%d, but must be %d)",
              vl_dst->max_vox,
              start_dst_index + num);

  vl_dst->type = vl_src->type;

  if (vl_src->mx) vl_dst->mx = (float *)calloc(num, sizeof(float));
  if (vl_src->my) vl_dst->my = (float *)calloc(num, sizeof(float));
  if (vl_src->mz) vl_dst->mz = (float *)calloc(num, sizeof(float));
  if (vl_src->t) vl_dst->t = (float *)calloc(num, sizeof(float));

  if (vl_src->mx && (vl_dst->t == NULL || vl_dst->mx == NULL || vl_dst->my == NULL || vl_dst->mz == NULL))
    ErrorExit(ERROR_NOMEMORY, "VLSTcopy: could not allocate %d-len slope arrays", num);

  for (i = 0; i < num; i++) {
    vl_dst->xi[start_dst_index + i] = vl_src->xi[i];
    vl_dst->yi[start_dst_index + i] = vl_src->yi[i];
    vl_dst->zi[start_dst_index + i] = vl_src->zi[i];
    vl_dst->fi[start_dst_index + i] = vl_src->fi[i];

    vl_dst->xd[start_dst_index + i] = vl_src->xd[i];
    vl_dst->yd[start_dst_index + i] = vl_src->yd[i];
    vl_dst->zd[start_dst_index + i] = vl_src->zd[i];

    vl_dst->vsrc[start_dst_index + i] = vl_src->vsrc[i];
    vl_dst->vdst[start_dst_index + i] = vl_src->vdst[i];
  }

  if (vl_dst->nvox < start_dst_index + num) vl_dst->nvox = start_dst_index + num;
  return (vl_dst);
}
int VLSTaddUnique(VOXEL_LIST *vl, int x, int y, int z, float xd, float yd, float zd)
{
  if (VLSTinList(vl, x, y, z)) return (0);
  VLSTadd(vl, x, y, z, xd, yd, zd);
  return (1);
}
int
VLSTaddWithValue(VOXEL_LIST *vl, int x, int y, int z, float xd, float yd, float zd, float vsrc, float vdst) 
{
  if (VLSTadd(vl,x,y,z,xd,yd,zd) == NO_ERROR)
  {
    int ind ;
    ind = vl->nvox-1 ;
    vl->vdst[ind] = vdst ;
    vl->vsrc[ind] = vsrc ;
  }
  else
    return(Gerror) ;
  return(NO_ERROR) ;
}

int VLSTadd(VOXEL_LIST *vl, int x, int y, int z, float xd, float yd, float zd)
{
  if (vl->nvox >= vl->max_vox)
    ErrorReturn(ERROR_NOMEMORY, (ERROR_NOMEMORY, "VLSTadd: too many voxels (%d)", vl->max_vox));
  vl->xi[vl->nvox] = x;
  vl->yi[vl->nvox] = y;
  vl->zi[vl->nvox] = z;
  vl->xd[vl->nvox] = xd;
  vl->yd[vl->nvox] = yd;
  vl->zd[vl->nvox] = zd;
  vl->nvox++;
  return (NO_ERROR);
}

int VLSTinList(VOXEL_LIST *vl, int x, int y, int z)
{
  int n;

  for (n = 0; n < vl->nvox; n++)
    if (vl->xi[n] == x && vl->yi[n] == y && vl->zi[n] == z) return (1);
  return (0);
}

VOXEL_LIST *VLSTsplineFit(VOXEL_LIST *vl, int num_control)
{
  VOXEL_LIST *vl_spline;
  int k, i, km1, kp1;
  float len, total_len, dx, dy, dz, tx, ty, tz;

  if (num_control > vl->nvox)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "VLSTsplineFit: input vl has %d points, not enough for %d control points",
                 vl->nvox,
                 num_control));

  vl_spline = VLSTalloc(num_control);
  vl_spline->type = VOXLIST_SPLINE;
  vl_spline->mx = (float *)calloc(num_control, sizeof(float));
  vl_spline->my = (float *)calloc(num_control, sizeof(float));
  vl_spline->mz = (float *)calloc(num_control, sizeof(float));
  vl_spline->t = (float *)calloc(num_control, sizeof(float));
  if (vl_spline->t == NULL || vl_spline->mx == NULL || vl_spline->my == NULL || vl_spline->mz == NULL)
    ErrorExit(ERROR_NOMEMORY, "VLSTsplineFit: could not allocate %d-len slope arrays", num_control);

  for (k = 0; k < num_control; k++) {
    i = nint((float)k * (float)(vl->nvox - 1) / (float)(num_control - 1));

    vl_spline->xi[k] = vl->xi[i];
    vl_spline->yi[k] = vl->yi[i];
    vl_spline->zi[k] = vl->zi[i];
    vl_spline->fi[k] = vl->fi[i];
    vl_spline->vsrc[k] = vl->vsrc[i];
    vl_spline->vdst[k] = vl->vdst[i];
    vl_spline->xd[k] = vl->xd[i];
    vl_spline->yd[k] = vl->yd[i];
    vl_spline->zd[k] = vl->zd[i];
  }

  // compute total length of control point line segments
  for (total_len = k = 0; k < num_control - 1; k++) {
    dx = vl_spline->xd[k + 1] - vl_spline->xd[k];
    dy = vl_spline->yd[k + 1] - vl_spline->yd[k];
    dz = vl_spline->zd[k + 1] - vl_spline->zd[k];
    len = sqrt(dx * dx + dy * dy + dz * dz);
    total_len += len;
  }

  // compute parameterization
  vl_spline->t[0] = 0;
  vl_spline->t[num_control - 1] = 1;
  for (k = 1; k < num_control - 1; k++) {
    dx = vl_spline->xd[k] - vl_spline->xd[k - 1];
    dy = vl_spline->yd[k] - vl_spline->yd[k - 1];
    dz = vl_spline->zd[k] - vl_spline->zd[k - 1];
    len = sqrt(dx * dx + dy * dy + dz * dz);
    vl_spline->t[k] = len / total_len;
  }

  // compute slopes
  for (len = k = 0; k < num_control; k++) {
    if (k == 0)
      km1 = 0;
    else
      km1 = k - 1;
    if (k == num_control - 1)
      kp1 = num_control - 1;
    else
      kp1 = k + 1;

    dx = vl_spline->xd[kp1] - vl_spline->xd[km1];
    dy = vl_spline->yd[kp1] - vl_spline->yd[km1];
    dz = vl_spline->zd[kp1] - vl_spline->zd[km1];
    tx = vl_spline->t[kp1] - vl_spline->t[km1];
    ty = vl_spline->t[kp1] - vl_spline->t[km1];
    tz = vl_spline->t[kp1] - vl_spline->t[km1];
    tx = ty = tz = 2;

    vl_spline->mx[k] = dx / tx;
    vl_spline->my[k] = dy / ty;
    vl_spline->mz[k] = dz / tz;
  }

  return (vl_spline);
}
int VLSTwriteLabel(VOXEL_LIST *vl, const char *fname, MRI_SURFACE *mris, MRI *mri)
{
  LABEL *area = VLSTtoLabel(vl, mris, mri);
  LabelWrite(area, fname);
  LabelFree(&area);
  return (NO_ERROR);
}
LABEL *VLSTtoLabel(VOXEL_LIST *vl, MRI_SURFACE *mris, MRI *mri)
{
  int n;
  LABEL *area = LabelAlloc(vl->nvox, NULL, "");
  double xs, ys, zs;

  if (mris) {
    for (n = 0; n < vl->nvox; n++) {
      MRISsurfaceRASFromVoxel(mris, mri, (double)vl->xd[n], (double)vl->yd[n], (double)vl->zd[n], &xs, &ys, &zs);
      area->lv[n].x = xs;
      area->lv[n].y = ys;
      area->lv[n].z = zs;
      area->lv[n].stat = vl->vsrc[n];
    }
  }
  else  // use scanner coords
  {
    for (n = 0; n < vl->nvox; n++) {
      MRIvoxelToWorld(mri, (double)vl->xd[n], (double)vl->yd[n], (double)vl->zd[n], &xs, &ys, &zs);
      area->lv[n].x = xs;
      area->lv[n].y = ys;
      area->lv[n].z = zs;
      area->lv[n].stat = vl->vsrc[n];
    }
    sprintf(area->space, "coords=scanner");
  }
  area->n_points = vl->nvox;

  return (area);
}

VOXEL_LIST *VLSTinterpolate(VOXEL_LIST *vl, float spacing)
{
  int k, nvox, km1, kp1;
  VOXEL_LIST *vl_interp;
  float x_k, y_k, z_k, x_kp1, y_kp1, z_kp1, x, y, z, t, dx, dy, dz, mx_k, my_k, mz_k, mx_kp1, my_kp1, mz_kp1, len;

  // compute slopes, and # of points in interpolated spline
  if (vl->mx == 0) {
    int km1, kp1;
    double dx, dy, dz;

    vl->mx = (float *)calloc(vl->nvox, sizeof(float));
    vl->my = (float *)calloc(vl->nvox, sizeof(float));
    vl->mz = (float *)calloc(vl->nvox, sizeof(float));
    vl->t = (float *)calloc(vl->nvox, sizeof(float));
    if (vl->t == NULL || vl->mx == NULL || vl->my == NULL || vl->mz == NULL)
      ErrorExit(ERROR_BADPARM, "VLSTinterpolate: cannot allocate slopes");
    for (k = 0; k < vl->nvox; k++) {
      if (k == 0)
        km1 = 0;
      else
        km1 = k - 1;
      if (k == vl->nvox - 1)
        kp1 = vl->nvox - 1;
      else
        kp1 = k + 1;
      dx = vl->xd[kp1] - vl->xd[km1];
      dy = vl->yd[kp1] - vl->yd[km1];
      dz = vl->zd[kp1] - vl->zd[km1];
      vl->mx[k] = dx / 2;
      vl->my[k] = dy / 2;
      vl->mz[k] = dz / 2;
    }
  }

  for (k = nvox = 0; k < vl->nvox - 1; k++) {
    if (k == 0)
      km1 = 0;
    else
      km1 = k - 1;
    if (k == vl->nvox - 1)
      kp1 = vl->nvox - 1;
    else
      kp1 = k + 1;

    dx = vl->xd[kp1] - vl->xd[km1];
    dy = vl->yd[kp1] - vl->yd[km1];
    dz = vl->zd[kp1] - vl->zd[km1];

    vl->mx[k] = dx / 2;
    vl->my[k] = dy / 2;
    vl->mz[k] = dz / 2;

    x_k = vl->xd[k];
    y_k = vl->yd[k];
    z_k = vl->zd[k];
    x_kp1 = vl->xd[k + 1];
    y_kp1 = vl->yd[k + 1];
    z_kp1 = vl->zd[k + 1];
    dx = x_kp1 - x_k;
    dy = y_kp1 - y_k;
    dz = z_kp1 - z_k;
    len = sqrt(dx * dx + dy * dy + dz * dz);
    if (FZERO(len)) continue;
    nvox += ceil(len / spacing) + 1;
  }

  vl_interp = VLSTalloc(nvox + 1);
  vl_interp->nvox = 0;
  for (k = 0; k < vl->nvox - 1; k++) {
    x_k = vl->xd[k];
    y_k = vl->yd[k];
    z_k = vl->zd[k];
    x_kp1 = vl->xd[k + 1];
    y_kp1 = vl->yd[k + 1];
    z_kp1 = vl->zd[k + 1];
    dx = x_kp1 - x_k;
    dy = y_kp1 - y_k;
    dz = z_kp1 - z_k;
    len = sqrt(dx * dx + dy * dy + dz * dz);
    if (FZERO(len)) continue;
    mx_k = vl->mx[k];
    my_k = vl->my[k];
    mz_k = vl->mz[k];
    mx_kp1 = vl->mx[k + 1];
    my_kp1 = vl->my[k + 1];
    mz_kp1 = vl->mz[k + 1];
    for (t = 0; t <= 1.0; t += spacing / len) {
      x = h00(t) * x_k + h10(t) * mx_k + h01(t) * x_kp1 + h11(t) * mx_kp1;
      y = h00(t) * y_k + h10(t) * my_k + h01(t) * y_kp1 + h11(t) * my_kp1;
      z = h00(t) * z_k + h10(t) * mz_k + h01(t) * z_kp1 + h11(t) * mz_kp1;
      VLSTadd(vl_interp, nint(x), nint(y), nint(z), x, y, z);
      if (nint(x) == Gx && nint(y) == Gy && nint(z) == Gz) DiagBreak();
    }
  }

  if (vl_interp->xi[vl_interp->nvox - 1] != vl->xi[vl->nvox - 1] ||
      vl_interp->yi[vl_interp->nvox - 1] != vl->yi[vl->nvox - 1] ||
      vl_interp->zi[vl_interp->nvox - 1] != vl->zi[vl->nvox - 1])
    VLSTadd(vl_interp,
            vl->xi[vl->nvox - 1],
            vl->yi[vl->nvox - 1],
            vl->zi[vl->nvox - 1],
            vl->xd[vl->nvox - 1],
            vl->yd[vl->nvox - 1],
            vl->zd[vl->nvox - 1]);
  return (vl_interp);
}

MRI *VLSTwriteOrderToMRI(VOXEL_LIST *vl, MRI *mri)
{
  int n;

  if (mri == NULL) {
    if (vl->mri == NULL)
      ErrorExit(ERROR_BADPARM, "VLSTwriteOrderToMRI: mri must be specified as parameter or in vl->mri");

    mri = MRIclone(vl->mri, NULL);
  }

  for (n = 0; n < vl->nvox; n++) MRIsetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0, n + 1);

  return (mri);
}
VOXEL_LIST *VLSTfromMRI(MRI *mri, int vno)
{
  VOXEL_LIST *vl;
  int n;
  double xd, yd, zd;

  vl = VLSTalloc(mri->height);
  vl->nvox = 0;
  for (n = 0; n < vl->max_vox; n++) {
    xd = MRIgetVoxVal(mri, vno, n, 0, 0);
    yd = MRIgetVoxVal(mri, vno, n, 1, 0);
    zd = MRIgetVoxVal(mri, vno, n, 2, 0);
    VLSTadd(vl, nint(xd), nint(yd), nint(zd), xd, yd, zd);
  }

  return (vl);
}
int VLSTinterpolateSplineIntoVolume(VOXEL_LIST *vl, MRI *mri, double spacing, VOXEL_LIST *vl_total, float val)
{
  VOXEL_LIST *vl_interp;

  vl_interp = VLSTinterpolate(vl, spacing);
  if (vl_interp == NULL) return (Gerror);
  VLSTinterpolateIntoVolume(vl_interp, mri, val);
  VLSTcopyInto(vl_interp, vl_total, vl_total->nvox, vl_interp->nvox);
  VLSTfree(&vl_interp);
  return (NO_ERROR);
}
int VLSTinterpolateIntoVolume(VOXEL_LIST *vl, MRI *mri, float val)
{
  int n;

  for (n = 0; n < vl->nvox; n++) MRIinterpolateIntoVolume(mri, vl->xd[n], vl->yd[n], vl->zd[n], val);

  return (NO_ERROR);
}

double VLSTcomputeEntropy(VOXEL_LIST *vl, MRI *mri, int num)
{
  double val, entropy, total;
  int n;

  for (total = 0.0, n = 0; n < vl->nvox; n++) {
    val = MRIgetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0);
    total += val;
  }

  for (entropy = 0.0, n = 0; n < vl->nvox; n++) {
    val = MRIgetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0);
    if (FZERO(val)) continue;
    MRIsetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0, 0);
    val /= total;
    entropy -= val * log(val);
  }

  return (entropy);
}

double VLSTcomputeSplineMean(VOXEL_LIST *vl_spline, MRI *mri, double step_size)
{
  VOXEL_LIST *vl = VLSTinterpolate(vl_spline, step_size);
  double mean, val;
  int n;

  for (mean = 0.0, n = 0; n < vl->nvox; n++) {
    val = MRIgetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0);
    if (FZERO(val)) continue;
    mean += val;
  }

  if (vl->nvox > 1) mean /= vl->nvox;
  VLSTfree(&vl);
  return (mean);
}

float VLSTcomputeSplineMedian(VOXEL_LIST *vl_spline, MRI *mri, double step_size)
{
  VOXEL_LIST *vl = VLSTinterpolate(vl_spline, step_size);
  float *vals, median;
  int n;

  vals = (float *)calloc(vl->nvox, sizeof(float));
  for (n = 0; n < vl->nvox; n++) {
    vals[n] = MRIgetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0);
    if (fabs(vals[n]) > 1e20) DiagBreak();
  }

  qsort(vals, vl->nvox, sizeof(float), compare_floats);
  median = vals[vl->nvox / 2];
  VLSTfree(&vl);
  free(vals);
  return (median);
}

double VLSTcomputeSplineSegmentMean(VOXEL_LIST *vl_spline, MRI *mri, double step_size, int start, int stop)
{
  VOXEL_LIST *vl = VLSTinterpolate(vl_spline, step_size);
  double mean, val;
  int n;

  if (stop >= vl->nvox - 1) stop = vl->nvox - 1;
  if (start > stop) start = stop;
  for (mean = 0.0, n = start; n <= stop; n++) {
    val = MRIgetVoxVal(mri, vl->xi[n], vl->yi[n], vl->zi[n], 0);
    if (FZERO(val)) continue;
    mean += val;
  }

  if (stop - start + 1) mean /= (float)(stop - start + 1);
  VLSTfree(&vl);
  return (mean);
}

double VLSTrmsDistance(VOXEL_LIST *vl1, VOXEL_LIST *vl2, double max_dist, MRI **pmri_dist)
{
  MRI *mri_dist;
  double dist, rms;
  int i;

  if (*pmri_dist == NULL) {
    MRI *mri_tmp;
    if (vl1->mri == NULL) ErrorExit(ERROR_BADPARM, "VLSTrmsDistance: vl1->mri must be set prior to calling");

    mri_tmp = VLSTtoMri(vl1, NULL);
    MRIbinarize(mri_tmp, mri_tmp, 1, 0, 1);
    *pmri_dist = mri_dist = MRIdistanceTransform(mri_tmp, NULL, 1, max_dist, DTRANS_MODE_SIGNED, NULL);
    MRIfree(&mri_tmp);
  }
  else
    mri_dist = *pmri_dist;  // already allocated and computed (hopefully!)

  for (rms = 0.0, i = 0; i < vl2->nvox; i++) {
    dist = MRIgetVoxVal(mri_dist, vl2->xi[i], vl2->yi[i], vl2->zi[i], 0);
    if (dist < 0) dist = 0;
    rms += dist * dist;
  }
  rms = sqrt(rms / vl2->nvox);
  return (rms);
}

double VLSThausdorffDistance(VOXEL_LIST *vl1, VOXEL_LIST *vl2, double max_dist, MRI **pmri_dist)
{
  MRI *mri_dist;
  double dist, hdist;
  int i;

  if (*pmri_dist == NULL) {
    MRI *mri_tmp;
    if (vl1->mri == NULL) ErrorExit(ERROR_BADPARM, "VLSTrmsDistance: vl1->mri must be set prior to calling");

    mri_tmp = VLSTtoMri(vl1, NULL);
    MRIbinarize(mri_tmp, mri_tmp, 1, 0, 1);
    *pmri_dist = mri_dist = MRIdistanceTransform(mri_tmp, NULL, 1, max_dist, DTRANS_MODE_SIGNED, NULL);
    MRIfree(&mri_tmp);
  }
  else
    mri_dist = *pmri_dist;  // already allocated and computed (hopefully!)

  hdist = 0;
  for (i = 0; i < vl2->nvox; i++) {
    dist = MRIgetVoxVal(mri_dist, vl2->xi[i], vl2->yi[i], vl2->zi[i], 0);
    if (dist < 0) dist = 0;
    if (dist > hdist) hdist = dist;
  }
  return (hdist);
}
double VLSTmean(VOXEL_LIST *vl, MRI *mri, double *pvar)
{
  double mean, var, val;
  int x, y, z, f, n;

  if (mri == NULL) mri = vl->mri;
  for (mean = var = 0.0, n = 0; n < vl->nvox; n++) {
    x = vl->xi[n];
    y = vl->yi[n];
    z = vl->zi[n];
    f = vl->fi[n];
    val = MRIgetVoxVal(mri, x, y, z, f);
    var += val * val;
    mean += val;
  }
  mean /= vl->nvox;
  var = (var - vl->nvox * mean * mean) / (vl->nvox - 1);
  if (pvar) *pvar = var;
  return (mean);
}

int VLSTsample(VOXEL_LIST *vl, MRI *mri)
{
  double val;
  int x, y, z, f, n;

  if (mri == NULL) mri = vl->mri;
  for (n = 0; n < vl->nvox; n++) {
    x = vl->xi[n];
    y = vl->yi[n];
    z = vl->zi[n];
    f = vl->fi[n];
    if (f < mri->nframes)
      val = MRIgetVoxVal(mri, x, y, z, f);
    else
      val = MRIgetVoxVal(mri, x, y, z, f);
    vl->vsrc[n] = val;
  }
  return (NO_ERROR);
}
int VLSTsampleFloat(VOXEL_LIST *vl, MRI *mri)
{
  double val;
  double x, y, z;
  int n, f;

  if (mri == NULL) mri = vl->mri;
  for (n = 0; n < vl->nvox; n++) {
    x = vl->xd[n];
    y = vl->yd[n];
    z = vl->zd[n];
    f = vl->fi[n];
    MRIsampleVolumeFrame(mri, x, y, z, f, &val);
    vl->vsrc[n] = val;
  }
  return (NO_ERROR);
}

int VLmostCommonLabel(VOXEL_LIST *vl)
{
  int max_label, *label_counts, n, label, max_label_count;

  max_label = 0;
  for (n = 0; n < vl->nvox; n++) {
    label = nint(vl->vsrc[n]);
    if (label > max_label) max_label = label;
  }
  label_counts = (int *)calloc(max_label + 1, sizeof(int));
  for (n = 0; n <= max_label; n++) label_counts[n] = 0;
  for (n = 0; n < vl->nvox; n++) {
    label = nint(vl->vsrc[n]);
    label_counts[label]++;
  }
  for (max_label_count = max_label = n = 0; n <= max_label; n++)
    if (label_counts[n] > max_label_count) {
      max_label_count = label_counts[n];
      max_label = n;
    }
  return (max_label);
}
