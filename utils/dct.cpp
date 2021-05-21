/**
 * @brief utilities for computing and applying a discrete cosine transform.
 *
 * utilities for computing and applying a discrete cosine transform.
 */
/*
 * Original Author: Bruce Fischl
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
#include <stdlib.h>
#include <string.h>

#include "dct.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "voxlist.h"

#define FOV 256.0

DCT *DCTalloc(int ncoef, MRI *mri_source)
{
  DCT *dct;
  dct = (DCT *)calloc(1, sizeof(DCT));
  if (dct == NULL) ErrorExit(ERROR_NOMEMORY, "DCTalloc(%d): allocation failed", ncoef);

  dct->b = 1.0;
  dct->ncoef = ncoef;
  dct->res = .25 * MIN(MIN(mri_source->xsize, mri_source->ysize), mri_source->zsize);
  dct->mri_source = mri_source;
  dct->v_xk = VectorAlloc(ncoef, MATRIX_REAL);
  dct->v_yk = VectorAlloc(ncoef, MATRIX_REAL);
  dct->v_zk = VectorAlloc(ncoef, MATRIX_REAL);
  if (dct->v_xk == NULL || dct->v_yk == NULL || dct->v_zk == NULL)
    ErrorExit(ERROR_NOMEMORY, "DCTalloc(%d): allocation failed", ncoef);
  dct->x_inv = (double *)calloc(ceil(mri_source->width / dct->res), sizeof(double));
  dct->y_inv = (double *)calloc(ceil(mri_source->height / dct->res), sizeof(double));
  dct->z_inv = (double *)calloc(ceil(mri_source->depth / dct->res), sizeof(double));
  dct->x = (double *)calloc(mri_source->width, sizeof(double));
  dct->y = (double *)calloc(mri_source->height, sizeof(double));
  dct->z = (double *)calloc(mri_source->depth, sizeof(double));
  if (dct->x_inv == NULL || dct->y_inv == NULL || dct->z_inv == NULL || dct->x == NULL || dct->y == NULL ||
      dct->z == NULL)
    ErrorExit(ERROR_NOMEMORY, "DCTalloc(%d): 2nd allocation failed", ncoef);
  DCTcreateMatrix(dct, mri_source, 0);
  DCTupdate(dct);
  return (dct);
}
int DCTfree(DCT **pdct)
{
  DCT *dct;

  dct = *pdct;
  *pdct = NULL;

  if (dct->m_x_basis) MatrixFree(&dct->m_x_basis);  // free old one
  if (dct->m_y_basis) MatrixFree(&dct->m_y_basis);  // free old one
  if (dct->m_z_basis) MatrixFree(&dct->m_z_basis);  // free old one
  if (dct->v_xk) VectorFree(&dct->v_xk);
  if (dct->v_yk) VectorFree(&dct->v_yk);
  if (dct->v_zk) VectorFree(&dct->v_zk);
  free(dct->x_inv);
  free(dct->y_inv);
  free(dct->z_inv);
  free(dct->x);
  free(dct->y);
  free(dct->z);
  free(dct);
  return (NO_ERROR);
}
int DCTcreateMatrix(DCT *dct, MRI *mri, int skip)
{
  int Nx, Ny, Nz, vox, k, N, i;
  double dk;
  MATRIX *m;

  if (dct->m_x_basis) MatrixFree(&dct->m_x_basis);  // free old one
  if (dct->m_y_basis) MatrixFree(&dct->m_y_basis);  // free old one
  if (dct->m_z_basis) MatrixFree(&dct->m_z_basis);  // free old one

  Nx = (int)((mri->width + skip) / (double)(skip + 1));
  Ny = (int)((mri->height + skip) / (double)(skip + 1));
  Nz = (int)((mri->depth + skip) / (double)(skip + 1));
  //  Nx = Ny = Nz = FOV ;
  dct->m_x_basis = MatrixAlloc(Nx, dct->ncoef, MATRIX_REAL);
  dct->m_y_basis = MatrixAlloc(Ny, dct->ncoef, MATRIX_REAL);
  dct->m_z_basis = MatrixAlloc(Nz, dct->ncoef, MATRIX_REAL);
  if (dct->m_x_basis == NULL || dct->m_y_basis == NULL || dct->m_z_basis == NULL)
    ErrorExit(ERROR_NOMEMORY, "DCTcreateMatrix(%d, %d, %d, %d,): allocation failed", dct->ncoef, Nx, Ny, Nz);

  for (i = 0; i < 3; i++)  // for each coordinate direction
  {
    switch (i) {
      default:
      case 0:
        N = Nx;
        m = dct->m_x_basis;
        break;
      case 1:
        N = Ny;
        m = dct->m_y_basis;
        break;
      case 2:
        N = Nz;
        m = dct->m_z_basis;
        break;
    }
    for (vox = 0; vox < N; vox++) {
      dk = sqrt(1.0 / N);
      *MATRIX_RELT(m, vox + 1, 1) = dk;
      for (k = 1; k < dct->ncoef; k++) {
        dk = sqrt(2.0 / N) * cos(M_PI * (2 * vox + 1) * k / (2.0 * N));
        ;
        *MATRIX_RELT(m, vox + 1, k + 1) = dk;
      }
    }
  }

  return (NO_ERROR);
}

DCT *DCTcopy(DCT *dct_src, DCT *dct_dst)
{
  if (dct_dst == NULL) dct_dst = DCTalloc(dct_src->ncoef, dct_src->mri_source);

  MatrixCopyRealRegion(dct_src->v_xk, dct_dst->v_xk, 1, 1, MIN(dct_src->v_xk->rows, dct_dst->v_xk->rows), 1, 1, 1);
  MatrixCopyRealRegion(dct_src->v_yk, dct_dst->v_yk, 1, 1, MIN(dct_src->v_yk->rows, dct_dst->v_yk->rows), 1, 1, 1);
  MatrixCopyRealRegion(dct_src->v_zk, dct_dst->v_zk, 1, 1, MIN(dct_src->v_zk->rows, dct_dst->v_zk->rows), 1, 1, 1);
  return (dct_dst);
}

MRI *DCTapply(DCT *dct, MRI *mri_src, MRI *mri_target, MRI *mri_dst, int sample_type)
{
  int x, y, z;
  double xd, yd, zd;
  double val;
  MATRIX *m_vox2vox;
  VECTOR *v1, *v2;

  if (mri_target == NULL) mri_target = mri_src;  // assume output geometry is same as input
  if (mri_dst == NULL) mri_dst = MRIclone(mri_target, NULL);
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_target, mri_src);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0;

  DCTupdate(dct);
  for (x = 0; x < mri_dst->width; x++) {
    V3_X(v1) = x;
    for (y = 0; y < mri_dst->height; y++) {
      V3_Y(v1) = y;
      for (z = 0; z < mri_dst->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        V3_Z(v1) = z;
        MatrixMultiply(m_vox2vox, v1, v2);
        DCTinverseTransformPoint(dct, V3_X(v2), V3_Y(v2), V3_Z(v2), &xd, &yd, &zd);
        MRIsampleVolumeType(mri_src, xd, yd, zd, &val, sample_type);
        MRIsetVoxVal(mri_dst, x, y, z, 0, val);
      }
    }
  }

  MatrixFree(&m_vox2vox);
  VectorFree(&v1);
  VectorFree(&v2);
  return (mri_dst);
}
MRI *DCTapplyInverse(DCT *dct, MRI *mri_src, MRI *mri_dst, int sample_type)
{
  int x, y, z;
  double xs, ys, zs;
  double val;

  if (mri_dst == NULL) mri_dst = MRIclone(mri_src, NULL);

  for (x = 0; x < mri_dst->width; x++) {
    for (y = 0; y < mri_dst->height; y++) {
      for (z = 0; z < mri_dst->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        DCTtransformPoint(dct, x, y, z, &xs, &ys, &zs);
        MRIsampleVolumeType(mri_src, xs, ys, zs, &val, sample_type);
        MRIsetVoxVal(mri_dst, x, y, z, 0, val);
      }
    }
  }

  return (mri_dst);
}

int DCTtransformVoxlist(DCT *dct, VOXEL_LIST *vl)
{
  int i;
  double x, y, z, xd, yd, zd;

  for (i = 0; i < vl->nvox; i++) {
    x = vl->xi[i];
    y = vl->yi[i];
    z = vl->zi[i];
    DCTtransformPoint(dct, x, y, z, &xd, &yd, &zd);
    vl->xd[i] = xd;
    vl->yd[i] = yd;
    vl->zd[i] = zd;
  }

  return (NO_ERROR);
}

int DCTinverseTransformVoxlist(DCT *dct, VOXEL_LIST *vl)
{
  int i, x, y, z;
  MATRIX *m_vox2vox;
  static VECTOR *v1 = NULL, *v2;
  double xd, yd, zd;

  if (DIAG_VERBOSE_ON) {
    MRIwrite(vl->mri, "t2.mgz");
    MRIwrite(dct->mri_source, "s2.mgz");
  }
  if (v1 == NULL) {
    v1 = VectorAlloc(4, MATRIX_REAL);
    v2 = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(v1, 4) = 1.0;
    VECTOR_ELT(v1, 4) = 1.0;
  }
  m_vox2vox = MRIgetVoxelToVoxelXform(vl->mri, dct->mri_source);
  for (i = 0; i < vl->nvox; i++) {
    V3_X(v1) = vl->xi[i];
    V3_Y(v1) = vl->yi[i];
    V3_Z(v1) = vl->zi[i];
    MatrixMultiply(m_vox2vox, v1, v2);
    x = nint(V3_X(v2));
    y = nint(V3_Y(v2));
    z = nint(V3_Z(v2));
    if (x < 0)
      x = 0;
    else if (x >= dct->mri_source->width)
      x = dct->mri_source->width - 1;
    if (y < 0)
      y = 0;
    else if (y >= dct->mri_source->height)
      y = dct->mri_source->height - 1;
    if (z < 0)
      z = 0;
    else if (z >= dct->mri_source->depth)
      z = dct->mri_source->depth - 1;

#if 0    
    vl->xd[i] = dct->x_inv[x] ;
    vl->yd[i] = dct->y_inv[y] ;
    vl->zd[i] = dct->z_inv[z] ;
#else
    DCTinverseTransformPoint(dct, x, y, z, &xd, &yd, &zd);
    vl->xd[i] = xd;
    vl->yd[i] = yd;
    zd = vl->zd[i];
#endif
  }

  MatrixFree(&m_vox2vox);
  return (NO_ERROR);
}
int DCTtransformPoint(DCT *dct, int x, int y, int z, double *px, double *py, double *pz)
{
  double xd, yd, zd;

  xd = MatrixRowDotProduct(dct->m_x_basis, x + 1, dct->v_xk);
  yd = MatrixRowDotProduct(dct->m_y_basis, y + 1, dct->v_yk);
  zd = MatrixRowDotProduct(dct->m_z_basis, z + 1, dct->v_zk);
  *px = x + xd;
  *py = y + yd;
  *pz = z + zd;
  return (NO_ERROR);
}

int DCTinverseTransformPoint(DCT *dct, double x, double y, double z, double *px, double *py, double *pz)
{
  double xd, yd, zd;
  int xi, yi, zi;

  xi = nint(x / dct->res);
  yi = nint(y / dct->res);
  zi = nint(z / dct->res);
  if (xi < 0) xi = 0;
  if (yi < 0) yi = 0;
  if (zi < 0) zi = 0;
  if (xi >= ceil(dct->mri_source->width / dct->res)) xi = ceil(dct->mri_source->width / dct->res) - 1;
  if (yi >= ceil(dct->mri_source->height / dct->res)) yi = ceil(dct->mri_source->height / dct->res) - 1;
  if (zi >= ceil(dct->mri_source->depth / dct->res)) zi = ceil(dct->mri_source->depth / dct->res) - 1;

  xd = dct->x_inv[xi];
  yd = dct->y_inv[yi];
  zd = dct->z_inv[zi];
  *px = xd;
  *py = yd;
  *pz = zd;
  return (NO_ERROR);
}

int DCTdump(DCT *dct, FILE *fp)
{
  int k;
  fprintf(fp, "DCT (%d):\n", dct->ncoef);
  for (k = 1; k <= dct->ncoef; k++)
    fprintf(
        fp, "\t(%2.3f, %2.3f, %2.3f)\n", VECTOR_ELT(dct->v_xk, k), VECTOR_ELT(dct->v_yk, k), VECTOR_ELT(dct->v_zk, k));
  fflush(fp);
  return (NO_ERROR);
}
static int soap_bubble(double *in_vals, double *ctrl, double *out_vals, int N, double max_change_allowed);
int DCTupdate(DCT *dct)
{
  double *x_wts, *y_wts, *z_wts, xd, yd, zd, jcd, jfd, *wts, *inv, *fwd, jd;
  int jc, jf, x, y, z, N, i, j, Ninv;

  x_wts = (double *)calloc(ceil(dct->mri_source->width / dct->res), sizeof(double));
  y_wts = (double *)calloc(ceil(dct->mri_source->height / dct->res), sizeof(double));
  z_wts = (double *)calloc(ceil(dct->mri_source->depth / dct->res), sizeof(double));
  if (x_wts == NULL || y_wts == NULL || z_wts == NULL)
    ErrorExit(ERROR_NOMEMORY, "DCTupdate(%d): allocation failed", dct->ncoef);

  for (i = 0; i < 3; i++) {
    switch (i) {
      default:
      case 0:
        wts = x_wts;
        N = dct->mri_source->width;
        inv = dct->x_inv;
        fwd = dct->x;
        break;
      case 1:
        wts = y_wts;
        N = dct->mri_source->height;
        inv = dct->y_inv;
        fwd = dct->y;
        break;
      case 2:
        wts = z_wts;
        N = dct->mri_source->depth;
        inv = dct->z_inv;
        fwd = dct->z;
        break;
    }
    Ninv = ceil(N / dct->res);
    memset(inv, 0, Ninv * sizeof(inv[0]));
    for (j = 0; j < N; j++) {
      x = y = z = 0;
      switch (i) {
        case 0:
          x = j;
          break;
        case 1:
          y = j;
          break;
        case 2:
          z = j;
          break;
      }
      DCTtransformPoint(dct, x, y, z, &xd, &yd, &zd);
      switch (i) {
        default:
        case 0:
          jd = xd / dct->res;
          break;
        case 1:
          jd = yd / dct->res;
          break;
        case 2:
          jd = zd / dct->res;
          break;
      }
      fwd[j] = jd;
      jf = (int)floor(jd);
      jc = jf + 1;
      jfd = jd - jf;
      jcd = 1 - jfd;
      if (jc >= 0 && jc < Ninv) {
        wts[jc] += jfd;
        inv[jc] += jfd * j;
      }
      if (jf >= 0 && jf < Ninv) {
        wts[jf] += jcd;
        inv[jf] += jcd * j;
      }
    }
    for (j = 0; j < Ninv; j++)
      if (wts[j] > 0) inv[j] /= wts[j];

    for (j = 0; j < Ninv; j++)
      if (wts[j] < 0) DiagBreak();

    soap_bubble(inv, wts, inv, Ninv, 0.01);
  }

  free(x_wts);
  free(y_wts);
  free(z_wts);
  return (NO_ERROR);
}
static int soap_bubble(double *in_vals, double *ctrl, double *out_vals, int N, double max_change_allowed)
{
  double *tmp_vals, max_change, change, out_val, min_val, max_val;
  int iter, i, j, num;

  tmp_vals = (double *)calloc(N, sizeof(double));
  memmove(out_vals, in_vals, N * sizeof(double));

  // fix first and last point to min and max index +- 1
  min_val = 1e8;
  max_val = -min_val;
  for (i = 0; i < N; i++)
    if (ctrl[i] > 0) {
      if (in_vals[i] > max_val) max_val = in_vals[i];
      if (in_vals[i] < min_val) min_val = in_vals[i];
    }
  if (DZERO(ctrl[0])) {
    out_vals[0] = min_val;
    ctrl[0] = 1;
  }
  if (DZERO(ctrl[N - 1])) {
    out_vals[N - 1] = max_val;
    ctrl[N - 1] = 1;
  }

  iter = 0;
  do {
    max_change = 0;
    memmove(tmp_vals, out_vals, N * sizeof(double));

    for (i = 0; i < N; i++) {
      if (ctrl[i] > 0) {
        out_vals[i] = tmp_vals[i];
        continue;  // fixed point
      }
      for (out_val = 0.0, num = 0, j = MAX(i - 1, 0); j <= MIN(i + 1, N); j++, num++) out_val += tmp_vals[j];
      out_vals[i] = out_val / num;
      change = fabs(out_vals[i] - tmp_vals[i]);
      if (change > max_change) max_change = change;
    }
    iter++;
    if (iter > 1000) break;
    memmove(tmp_vals, out_vals, N * sizeof(double));
  } while (max_change > 0.001);

  free(tmp_vals);
  return (NO_ERROR);
}
