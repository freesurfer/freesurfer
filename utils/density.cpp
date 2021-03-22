/**
 * @brief utilities for computing joint densities from images
 *
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

#include <stdio.h>
#include <stdlib.h>

#include "density.h"
#include "diag.h"
#include "error.h"
#include "image.h"
#include "macros.h"  // DZERO
#include "mri.h"
#include "proto.h"  // nint
#include "utils.h"  // fgetl

DENSITY *DensityHistogramEstimate(MRI *mri1, MRI *mri2, int nbins, float sigma, int *valid1, int *valid2)
{
  DENSITY *pdf;
  int x, y, bin1, bin2, nvox, n1, n2, n;
  float val1, val2, scale1, scale2;

  pdf = (DENSITY *)calloc(1, sizeof(DENSITY));
  if (!pdf) ErrorExit(ERROR_NOMEMORY, "DensityHistogramEstimate(%d): could not allocate pdf", nbins);
  pdf->sigma = sigma;
  pdf->Ipdf = ImageAlloc(nbins, nbins, PFFLOAT, 1);
  if (pdf->Ipdf == NULL)
    ErrorExit(ERROR_NOMEMORY, "DensityHistogramEstimate(%d): could not allocate density image", nbins);

  MRIvalRange(mri1, &pdf->min_val1, &pdf->max_val1);
  MRIvalRange(mri2, &pdf->min_val2, &pdf->max_val2);

  pdf->max_val1 = MAX(255, pdf->max_val1);
  pdf->max_val2 = MAX(255, pdf->max_val2);
  scale1 = (float)(nbins - 1) / (pdf->max_val1 - pdf->min_val1);
  scale2 = (float)(nbins - 1) / (pdf->max_val2 - pdf->min_val2);

  strcpy(pdf->fname1, mri1->fname);
  strcpy(pdf->fname2, mri2->fname);
  n1 = ceil(pdf->max_val1 - pdf->min_val1 + 1);
  n2 = ceil(pdf->max_val2 - pdf->min_val2 + 1);
  pdf->valid1 = (int *)calloc(n1, sizeof(int));
  pdf->valid2 = (int *)calloc(n2, sizeof(int));
  if (!pdf->valid1 || !pdf->valid2)
    ErrorExit(ERROR_NOMEMORY, "DensityHistogramEstimate: could not allocate lookup tables (%d, %d)\n", n1, n2);
  if (valid1) {
    for (n = 0; n < n1; n++) pdf->valid1[n] = valid1[n];
  }
  else
    memset(pdf->valid1, 1, n1 * sizeof(int));

  if (valid2) {
    for (n = 0; n < n2; n++) pdf->valid2[n] = valid2[n];
  }
  else
    memset(pdf->valid2, 1, n2 * sizeof(int));

  for (nvox = 0, x = 0; x < mri1->width; x++) {
    for (y = 0; y < mri1->height; y++) {
      if (x == Gx && y == Gy) DiagBreak();
      val1 = MRIgetVoxVal(mri1, x, y, 0, 0);
      if (valid1 && (valid1[nint(val1 - pdf->min_val1)] == 0)) continue;
      val2 = MRIgetVoxVal(mri2, x, y, 0, 0);
      if (valid2 && (valid2[nint(val2 - pdf->min_val2)] == 0)) continue;
      nvox++;
      bin1 = (val1 - pdf->min_val1) * scale1;
      bin2 = (val2 - pdf->min_val2) * scale2;
      if (bin1 < 0 || bin1 >= nbins || bin2 < 0 || bin2 >= nbins) {
        printf("illegal bin!!!! (%d, %d)\n", bin1, bin2);
        DiagBreak();
      }
      *IMAGEFpix(pdf->Ipdf, bin1, bin2) += 1.0;
    }
  }

  for (x = 0; x < nbins; x++) {
    for (y = 0; y < nbins; y++) {
      *IMAGEFpix(pdf->Ipdf, x, y) /= nvox;
    }
  }

  pdf->dof = nvox;
  if (sigma > 0) {
    IMAGE *Ikernel, *Ismooth;

    if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) ImageWrite(pdf->Ipdf, "Ipdf.tif");
    Ikernel = ImageGaussian1d(sigma, 0);
    Ismooth = ImageConvolveGaussian(pdf->Ipdf, Ikernel, NULL, 0);
    ImageFree(&pdf->Ipdf);
    ImageFree(&Ikernel);
    pdf->Ipdf = Ismooth;
    if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) ImageWrite(pdf->Ipdf, "Ipdf_smooth.tif");
  }

  printf("histogram computed with %d dofs\n", nvox);
  return (pdf);
}
int DensityWrite(DENSITY *pdf, char *fname)
{
  FILE *fp;
  int i, j, n;
#if 0
  char buf[STRLEN] ;
  sprintf(buf, "%s.tif", fname) ;

  ImageWrite(pdf->Ipdf, buf) ;
#endif

  fp = fopen(fname, "w");
  if (fp == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "DensityWrite(%s): could not open file", fname));

  fprintf(fp, ":nbins %d\n", pdf->Ipdf->rows);
  fprintf(fp, ":vol1 %f %f\n", pdf->min_val1, pdf->max_val1);
  fprintf(fp, ":vol2 %f %f\n", pdf->min_val2, pdf->max_val2);
  fprintf(fp, ":sigma %f\n", pdf->sigma);
  fprintf(fp, ":dof %d\n", pdf->dof);
  fprintf(fp, ":mri1 %s\n", pdf->fname1);
  fprintf(fp, ":mri2 %s\n", pdf->fname2);
  fprintf(fp, ":valid1 ");
  n = ceil(pdf->max_val1 - pdf->min_val1 + 1);
  for (i = 0; i < n; i++) fprintf(fp, "%d ", pdf->valid1[i]);
  fprintf(fp, "\n");
  fprintf(fp, ":valid2 ");
  n = ceil(pdf->max_val2 - pdf->min_val2 + 1);
  for (i = 0; i < n; i++) fprintf(fp, "%d ", pdf->valid2[i]);
  fprintf(fp, "\n");

  for (i = 0; i < pdf->Ipdf->rows; i++) {
    for (j = 0; j < pdf->Ipdf->cols; j++) {
      fprintf(fp, "%f ", *IMAGEFpix(pdf->Ipdf, i, j));
    }
    fprintf(fp, "\n");
  }

  return (NO_ERROR);
}
DENSITY *DensityRead(char *fname)
{
  DENSITY *pdf;
  FILE *fp;
  int i, j, nbins, n;
  char *cp, line[MAX_LINE_LEN];

  pdf = (DENSITY *)calloc(1, sizeof(DENSITY));
  if (!pdf) ErrorExit(ERROR_NOMEMORY, "DensityRead(%s): could not density", fname);

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_BADPARM, "DensityRead(%s): could not open file", fname));

  cp = fgetl(line, MAX_LINE_LEN, fp);
  i = sscanf(cp, "%*s %d\n", &nbins);
  pdf->Ipdf = ImageAlloc(nbins, nbins, PFFLOAT, 1);
  if (pdf->Ipdf == NULL) ErrorExit(ERROR_NOMEMORY, "DensityRead(%d): could not allocate density image", nbins);
  if (fscanf(fp,
             ":vol1 %f %f\n:vol2 %f %f\n:sigma %f\n:dof %d\n:mri1 %s\n:mri2 %s\n:valid1 ",
             &pdf->min_val1,
             &pdf->max_val1,
             &pdf->min_val2,
             &pdf->max_val2,
             &pdf->sigma,
             &pdf->dof,
             pdf->fname1,
             pdf->fname2) != 8) {
    ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): could not read parameter(s)", fname);
  }

  n = ceil(pdf->max_val1 - pdf->min_val1 + 1);
  pdf->valid1 = (int *)calloc(n, sizeof(int));
  if (!pdf->valid1) ErrorExit(ERROR_NOMEMORY, "DensityHistogramEstimate: could not allocate lookup tables (%d)\n", n);
  for (i = 0; i < n; i++) {
    if (fscanf(fp, "%d ", &pdf->valid1[i]) != 1) {
      ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): could not read parameter(s)", fname);
    }
  }
  if (fscanf(fp, "\n:valid2 ") != 0) {
    ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): could not read expected line", fname);
  }
  n = ceil(pdf->max_val2 - pdf->min_val2 + 1);
  pdf->valid2 = (int *)calloc(n, sizeof(int));
  if (!pdf->valid2) ErrorExit(ERROR_NOMEMORY, "DensityHistogramEstimate: could not allocate lookup tables (%d)\n", n);
  for (i = 0; i < n; i++) {
    if (fscanf(fp, "%d ", &pdf->valid2[i]) != 1) {
      ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): could not read parameter(s)", fname);
    }
  }
  if (fscanf(fp, "\n") != 0) {
    ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): could not read expected line", fname);
  }

  pdf->min_p = 1.0;
  for (i = 0; i < pdf->Ipdf->rows; i++) {
    for (j = 0; j < pdf->Ipdf->cols; j++) {
      if (fscanf(fp, "%f ", IMAGEFpix(pdf->Ipdf, i, j)) != 1) {
        ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): could not read parameter(s)", fname);
      }
      if (!DZERO(*IMAGEFpix(pdf->Ipdf, i, j)) && (*IMAGEFpix(pdf->Ipdf, i, j) < pdf->min_p)) {
        pdf->min_p = *IMAGEFpix(pdf->Ipdf, i, j);
      }
    }
    if (fscanf(fp, "\n") != 0) {
      ErrorPrintf(ERROR_BAD_FILE, "DensityRead(%s): reached end of file", fname);
    }
  }

  return (pdf);
}
#define BIG_AND_NEGATIVE -1000000
double DensityLogLikelihood(DENSITY *pdf, float val1, float val2)
{
  int bin1, bin2, nbins;
  double p;

  nbins = pdf->Ipdf->rows;
  bin1 = (val1 - pdf->min_val1) * (float)(nbins - 1) / (pdf->max_val1 - pdf->min_val1);
  bin2 = (val2 - pdf->min_val2) * (float)(nbins - 1) / (pdf->max_val2 - pdf->min_val2);
#if 0
  if ((pdf->valid1[bin1] == 0) && (pdf->valid2[bin2] == 0))
    return(0) ;
#endif

  if (bin1 < 0 || bin1 >= nbins || bin2 < 0 || bin2 >= nbins) return (log(pdf->min_p / 10));
  p = *IMAGEFpix(pdf->Ipdf, bin1, bin2);
  if (DZERO(p)) return (log(pdf->min_p / 10));
  return (log(p));
}

MRI *DensityLikelihoodImage(MRI *mri1, MRI *mri2, MRI *mri_ll, MATRIX *m, DENSITY *pdf, MRI *mri_seg, int inverse)
{
  double val_src, val_dst, xs, ys, thick, ll;
  int x, y;
  VECTOR *v_src, *v_dst;

  thick = mri1->xsize;
  v_src = VectorAlloc(3, MATRIX_REAL);
  v_dst = VectorAlloc(3, MATRIX_REAL);

  if (!mri_ll) {
    mri_ll = MRIalloc(mri2->width, mri2->height, mri2->depth, MRI_FLOAT);
    MRIcopyHeader(mri2, mri_ll);
  }

  VECTOR_ELT(v_src, 3) = 1.0;
  VECTOR_ELT(v_dst, 3) = 1.0 / thick;

  for (x = 0; x < mri2->width; x++) {
    VECTOR_ELT(v_dst, 1) = (double)x;
    for (y = 0; y < mri2->height; y++) {
      if (x == Gx && y == Gy) DiagBreak();
      if (mri_seg && MRIvox(mri_seg, x, y, 0) == 0) continue;

      VECTOR_ELT(v_dst, 2) = (double)y;
      MatrixMultiply(m, v_dst, v_src);
      xs = VECTOR_ELT(v_src, 1);
      ys = VECTOR_ELT(v_src, 2);

      if (MRIindexNotInVolume(mri1, xs, ys, 0) == 0) {
        val_dst = MRIvox(mri2, x, y, 0);
        MRIsampleVolumeSlice(mri1, xs, ys, 0, &val_src, MRI_CORONAL);
      }
      else {
        val_src = val_dst = -1000000;
      }
      if (inverse)
        ll = -DensityLogLikelihood(pdf, val_src, val_dst); /* nissl,block */
      else
        ll = -DensityLogLikelihood(pdf, val_dst, val_src); /* nissl,block */
      MRIsetVoxVal(mri_ll, x, y, 0, 0, ll);
    }
  }

  VectorFree(&v_src);
  VectorFree(&v_dst);
  return (mri_ll);
}
