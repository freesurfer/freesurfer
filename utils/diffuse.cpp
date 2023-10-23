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

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "const.h"
#include "diag.h"
#include "error.h"
#include "image.h"
#include "kernel.h"
#include "macros.h"
#include "proto.h"

/* sobel x and y filter coefficients */
float sy[9] = {-0.25f, -0.5f, -0.25f, 0.0f, 0.0f, 0.0f, 0.25f, 0.5f, 0.25f};
float sx[9] = {-0.25f, 0.0f, 0.25f, -0.5f, 0.0f, 0.5f, -0.25f, 0.0f, 0.25f};

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#define g(I, k) exp((-1.0 / (k)) * fabs((double)(I)))
#define LAMBDA 0.25
int ImageDiffuse(IMAGE *inImage, IMAGE *outImage, double k, int niter, int which, double slope, KIMAGE *kimage)
{
  IMAGE *fSrcImage, *fDstImage;

  if (outImage->pixel_format != inImage->pixel_format) {
    fprintf(stderr, "ImageDiffuse: input and output image format must match.\n");
    exit(2);
  }

  if (!ImageCheckSize(inImage, outImage, 0, 0, 0)) {
    fprintf(stderr, "ImageDiffuse: input and output image sizes must match.\n");
    exit(2);
  }

  outImage->rows = inImage->rows;
  outImage->cols = inImage->cols;

  fSrcImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  fDstImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);

  ImageCopy(inImage, fSrcImage);

  switch (which) {
    case DIFFUSE_PERONA:
      ImageDiffusePerona(fSrcImage, fDstImage, k, niter, slope, kimage);
      break;
    case DIFFUSE_CURVATURE:
      ImageDiffuseCurvature(fSrcImage, fDstImage, k, niter, slope, kimage);
      break;
    case DIFFUSE_HV:
      ImageDiffuseHV(fSrcImage, fDstImage, k, niter, slope, kimage);
      break;
    default:
      fprintf(stderr, "ImageDiffuse: unknown diffusion type %d\n", which);
      break;
  }

  ImageCopy(fDstImage, outImage);
  ImageFree(&fSrcImage);
  ImageFree(&fDstImage);

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#define FOUR_CONNECTED 0
#if FOUR_CONNECTED
#define KERNEL_MUL 0.25f
#else
#define KERNEL_MUL 0.125f
#endif

int ImageDiffuseCurvature(IMAGE *inImage, IMAGE *outImage, double A, int niter, double slope, KIMAGE *kimage)
{
  int x, y, i, rows, cols, *xE, *yN, *xW, *yS, ys, yn, xe, xw, ci;
  float c[9], fvals[9], dst_val;
  FILE *fp;
  KIMAGE *ksrc = NULL;
  static IMAGE *tmpImage = NULL, *gradImage = NULL;

  rows = inImage->rows;
  cols = inImage->cols;

  if (kimage) ksrc = KernelImageClone(kimage);

  if ((outImage->pixel_format != inImage->pixel_format) || (outImage->pixel_format != PFFLOAT)) {
    fprintf(stderr, "ImageDiffuseCurv: input and output image format must both be float.\n");
    exit(2);
  }

  if ((outImage->rows != inImage->rows) || (outImage->cols != inImage->cols)) {
    fprintf(stderr, "ImageDiffuseCurv: input and output image sizes must match.\n");
    exit(2);
  }

  if (!ImageCheckSize(inImage, tmpImage, 0, 0, 0)) {
    if (tmpImage) ImageFree(&tmpImage);
    tmpImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  }
  if (!ImageCheckSize(inImage, gradImage, 0, 0, 0)) {
    if (gradImage) ImageFree(&gradImage);
    gradImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  }

  /* build index tables for border */
  xE = (int *)calloc((unsigned int)cols, sizeof(int));
  xW = (int *)calloc((unsigned int)cols, sizeof(int));
  yN = (int *)calloc((unsigned int)rows, sizeof(int));
  yS = (int *)calloc((unsigned int)rows, sizeof(int));

  xW[0] = 0;
  for (x = 1; x < cols; x++) {
    xW[x] = x - 1;
    xE[x - 1] = x;
  }
  xE[cols - 1] = cols - 1;
  yS[0] = 0;
  for (y = 1; y < rows; y++) {
    yN[y] = y - 1;
    yS[y - 1] = y;
  }
  yS[rows - 1] = rows - 1;

  ImageCopy(inImage, tmpImage);

#if 0
  if (0 && (Gdiag & DIAG_WRITE))
    fp = fopen("diffuse.dat", "w") ;
  else
#endif
  fp = NULL;

#if FOUR_CONNECTED
  fvals[1] = c[1] = 0.0f;
  fvals[3] = c[3] = 0.0f;
  fvals[6] = c[6] = 0.0f;
  fvals[8] = c[8] = 0.0f;
#endif

  for (i = 0; i < niter; i++) {
#if 0
    if (Gdiag & DIAG_WRITE)
      fprintf(stdout, "iteration %d\n", i) ;
#endif

    if (kimage) KernelImageCopy(kimage, ksrc);

    ImageCurvature(tmpImage, (float)A, gradImage);
    for (x = 0; x < cols; x++) {
      xe = xE[x];
      xw = xW[x];
      for (y = 0; y < rows; y++) {
        /*
          C1 |  C2  | C3
          ---------------
          C4 |  C0  | C5
          ---------------
          C6 |  C7  | C8
        */

        yn = yN[y];
        ys = yS[y];

        fvals[0] = *IMAGEFpix(tmpImage, x, y);
        fvals[2] = *IMAGEFpix(tmpImage, x, yn);
        fvals[4] = *IMAGEFpix(tmpImage, xw, y);
        fvals[5] = *IMAGEFpix(tmpImage, xe, y);
        fvals[7] = *IMAGEFpix(tmpImage, x, ys);

#if !FOUR_CONNECTED
        fvals[1] = *IMAGEFpix(tmpImage, xw, yn);
        fvals[3] = *IMAGEFpix(tmpImage, xe, yn);
        fvals[6] = *IMAGEFpix(tmpImage, xw, ys);
        fvals[8] = *IMAGEFpix(tmpImage, xe, ys);
#endif

        c[2] = KERNEL_MUL / *IMAGEFpix(gradImage, x, yn);
        c[4] = KERNEL_MUL / *IMAGEFpix(gradImage, xw, y);
        c[5] = KERNEL_MUL / *IMAGEFpix(gradImage, xe, y);
        c[7] = KERNEL_MUL / *IMAGEFpix(gradImage, x, ys);
#if !FOUR_CONNECTED
        c[1] = KERNEL_MUL / *IMAGEFpix(gradImage, xw, yn);
        c[3] = KERNEL_MUL / *IMAGEFpix(gradImage, xe, yn);
        c[6] = KERNEL_MUL / *IMAGEFpix(gradImage, xw, ys);
        c[8] = KERNEL_MUL / *IMAGEFpix(gradImage, xe, ys);
#endif

        c[0] = 1.0f;
        for (ci = 1; ci <= 8; ci++) c[0] -= c[ci];

        if (kimage) /* track evolution of kernel */
        {
          /* kernel stuff is row,col */
          KernelDiscount(kimage, y, x, c[0]);
          KernelUpdate(ksrc, kimage, y, x, yn, x, c[2]);
          KernelUpdate(ksrc, kimage, y, x, y, xw, c[4]);
          KernelUpdate(ksrc, kimage, y, x, y, xe, c[5]);
          KernelUpdate(ksrc, kimage, y, x, ys, x, c[7]);
#if !FOUR_CONNECTED
          KernelUpdate(ksrc, kimage, y, x, yn, xw, c[1]);
          KernelUpdate(ksrc, kimage, y, x, yn, xe, c[3]);
          KernelUpdate(ksrc, kimage, y, x, ys, xw, c[6]);
          KernelUpdate(ksrc, kimage, y, x, ys, xe, c[8]);
#endif
        }

        for (dst_val = 0.0f, ci = 0; ci <= 8; ci++) dst_val += fvals[ci] * c[ci];
        *IMAGEFpix(outImage, x, y) = dst_val;
      }
    }

    ImageCopy(outImage, tmpImage);
    A = A + A * slope;
    if (kimage) KernelImageNormalize(kimage);
  }

  free(xE);
  free(xW);
  free(yN);
  free(yS);

  if (fp) fclose(fp);

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#if 0
#define Chv(dx, dy, k) ((float)exp(-0.5f * SQR((float)(fabs(dx) - fabs(dy)) / k)))
#else
#define Chv(adx, ady, k) ((float)exp(-(float)SQR((adx - ady)) / k))
#endif

#if 0
#define KSIZE 9
static float xkernel[3*9] =
  {
    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
    -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f
  } ;

static float ykernel[9*3] =
  {
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f,
    1.0f, 0.0f, -1.0f
  } ;
#else

#define KSIZE 15
static float xkernel[3 * KSIZE];
static float ykernel[KSIZE * 3];

#endif

int ImageDiffuseHV(IMAGE *inImage, IMAGE *outImage, double k, int niter, double slope, KIMAGE *kimage)
{
  int x, y, i, rows, cols, *xE, *yN, *xW, *yS, ys, yn, xe, xw, ci;
  float c[9], fvals[9], dst_val;
  FILE *fp;
  KIMAGE *ksrc = NULL;
  static IMAGE *xImage = NULL, *yImage = NULL, *tmpImage;
  IMAGE xKernelImage, yKernelImage;

  yKernelImage.cols = yKernelImage.ocols = xKernelImage.cols = xKernelImage.ocols = 3;
  yKernelImage.cols = yKernelImage.ocols = xKernelImage.cols = xKernelImage.ocols = 3;
  yKernelImage.image = yKernelImage.firstpix = (ubyte *)ykernel;
  xKernelImage.image = xKernelImage.firstpix = (ubyte *)xkernel;
  if (tmpImage == NULL) /* initialize kernels */
  {
    int row, col;
    float *kptr;

    kptr = xkernel;
    for (col = 0; col < KSIZE; col++, kptr++) {
#if 0
      *kptr = 1.0f/(float)KSIZE ;
      *(kptr + KSIZE) = 0.0f ;
      *(kptr + 2*KSIZE) = -1.0f/(float)KSIZE ;
#else
      *kptr = 1.0f;
      *(kptr + KSIZE) = 0.0f;
      *(kptr + 2 * KSIZE) = -1.0f;
#endif
    }
    kptr = ykernel;
    for (row = 0; row < KSIZE; row++, kptr += 3) {
#if 0
      *kptr = 1.0f/(float)KSIZE ;
      *(kptr + 1) = 0.0f ;
      *(kptr + 2) = -1.0f/(float)KSIZE ;
#else
      *kptr = 1.0f;
      *(kptr + 1) = 0.0f;
      *(kptr + 2) = -1.0f;
#endif
    }
  }

  rows = inImage->rows;
  cols = inImage->cols;

  if (kimage) ksrc = KernelImageClone(kimage);

  if ((outImage->pixel_format != inImage->pixel_format) || (outImage->pixel_format != PFFLOAT)) {
    fprintf(stderr,
            "ImageDiffuseHV: input and output image format must both be "
            "in float format.\n");
    exit(2);
  }

  if ((outImage->rows != inImage->rows) || (outImage->cols != inImage->cols)) {
    fprintf(stderr, "ImageDiffuseHV: input and output image sizes must match.\n");
    exit(2);
  }

  if (!ImageCheckSize(inImage, tmpImage, 0, 0, 0)) {
    if (tmpImage) ImageFree(&tmpImage);
    tmpImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  }

  if (!ImageCheckSize(inImage, xImage, 0, 0, 0)) {
    if (xImage) ImageFree(&xImage);
    xImage = ImageAlloc(rows, cols, PFFLOAT, 1);
  }
  else {
    xImage->rows = rows;
    xImage->cols = cols;
  }
  if (!ImageCheckSize(inImage, yImage, 0, 0, 0)) {
    if (yImage) ImageFree(&yImage);
    yImage = ImageAlloc(rows, cols, PFFLOAT, 1);
  }
  else {
    yImage->rows = rows;
    yImage->cols = cols;
  }

  /* build index tables for border */
  xE = (int *)calloc((unsigned int)cols, sizeof(int));
  xW = (int *)calloc((unsigned int)cols, sizeof(int));
  yN = (int *)calloc((unsigned int)rows, sizeof(int));
  yS = (int *)calloc((unsigned int)rows, sizeof(int));

  xW[0] = 0;
  for (x = 1; x < cols; x++) {
    xW[x] = x - 1;
    xE[x - 1] = x;
  }
  xE[cols - 1] = cols - 1;
  yS[0] = 0;
  for (y = 1; y < rows; y++) {
    yN[y] = y - 1;
    yS[y - 1] = y;
  }
  yS[rows - 1] = rows - 1;

  ImageCopy(inImage, tmpImage);

#if 0
  if (0 && (Gdiag & DIAG_WRITE))
    fp = fopen("diffuse.dat", "w") ;
  else
#endif
  fp = NULL;

  for (i = 0; i < niter; i++) {
#if 0
    if (Gdiag & DIAG_WRITE)
      fprintf(stdout, "iteration %d\n", i) ;
#endif

    if (kimage) KernelImageCopy(kimage, ksrc);

#if 0
    ImageConvolve3x3(inImage, sx, xImage) ;
    ImageConvolve3x3(inImage, sy, yImage) ;
#else
    ImageCorrelate(&xKernelImage, inImage, 0, xImage);
    ImageAbs(xImage, xImage);
    ImageCorrelate(&yKernelImage, inImage, 0, yImage);
    ImageAbs(yImage, yImage);
#endif

    for (x = 0; x < cols; x++) {
      xe = xE[x];
      xw = xW[x];
      for (y = 0; y < rows; y++) {
        /*
          C1 |  C2  | C3
          ---------------
          C4 |  C0  | C5
          ---------------
          C6 |  C7  | C8
        */

        yn = yN[y];
        ys = yS[y];

        fvals[0] = *IMAGEFpix(tmpImage, x, y);
        fvals[2] = *IMAGEFpix(tmpImage, x, yn);
        fvals[4] = *IMAGEFpix(tmpImage, xw, y);
        fvals[5] = *IMAGEFpix(tmpImage, xe, y);
        fvals[7] = *IMAGEFpix(tmpImage, x, ys);
#if !FOUR_CONNECTED
        fvals[1] = *IMAGEFpix(tmpImage, xw, yn);
        fvals[3] = *IMAGEFpix(tmpImage, xe, yn);
        fvals[6] = *IMAGEFpix(tmpImage, xw, ys);
        fvals[8] = *IMAGEFpix(tmpImage, xe, ys);
#endif

        c[2] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, x, yn), *IMAGEFpix(yImage, x, yn), k);
        c[4] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, xw, y), *IMAGEFpix(yImage, xw, y), k);
        c[5] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, xe, y), *IMAGEFpix(yImage, xe, y), k);
        c[7] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, x, ys), *IMAGEFpix(yImage, x, ys), k);
#if !FOUR_CONNECTED
        c[1] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, xw, yn), *IMAGEFpix(yImage, xw, yn), k);
        c[3] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, xe, yn), *IMAGEFpix(yImage, xe, yn), k);
        c[6] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, xw, ys), *IMAGEFpix(yImage, xw, ys), k);
        c[8] = KERNEL_MUL * Chv(*IMAGEFpix(xImage, xe, ys), *IMAGEFpix(yImage, xe, ys), k);
#endif

        c[0] = 1.0f;
        for (ci = 1; ci <= 8; ci++) c[0] -= c[ci];

        if (kimage) /* track evolution of kernel */
        {
          /* kernel stuff is row,col */
          KernelDiscount(kimage, y, x, c[0]);
          KernelUpdate(ksrc, kimage, y, x, yn, x, c[2]);
          KernelUpdate(ksrc, kimage, y, x, y, xw, c[4]);
          KernelUpdate(ksrc, kimage, y, x, y, xe, c[5]);
          KernelUpdate(ksrc, kimage, y, x, ys, x, c[7]);
#if !FOUR_CONNECTED
          KernelUpdate(ksrc, kimage, y, x, yn, xw, c[1]);
          KernelUpdate(ksrc, kimage, y, x, yn, xe, c[3]);
          KernelUpdate(ksrc, kimage, y, x, ys, xw, c[6]);
          KernelUpdate(ksrc, kimage, y, x, ys, xe, c[8]);
#endif
        }

        for (dst_val = 0.0f, ci = 0; ci <= 8; ci++) dst_val += fvals[ci] * c[ci];
        *IMAGEFpix(outImage, x, y) = dst_val;
      }
    }

    ImageCopy(outImage, tmpImage);
    k = k + k * slope;
    if (kimage) KernelImageNormalize(kimage);
  }

  free(xE);
  free(xW);
  free(yN);
  free(yS);

  if (fp) fclose(fp);

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#if 1

#define C(grad, k) ((float)exp(-0.5f * SQR((float)fabs((grad)) / k)))

int ImageDiffusePerona(IMAGE *inImage, IMAGE *outImage, double k, int niter, double slope, KIMAGE *kimage)
{
  int x, y, i, rows, cols, *xE, *yN, *xW, *yS, ys, yn, xe, xw, ci;
  float c[9], fvals[9], dst_val;
  FILE *fp;
  KIMAGE *ksrc = NULL;
  static IMAGE *tmpImage = NULL, *gradImage = NULL;

  rows = inImage->rows;
  cols = inImage->cols;

  if (kimage) ksrc = KernelImageClone(kimage);

  if ((outImage->pixel_format != inImage->pixel_format) || (outImage->pixel_format != PFFLOAT)) {
    fprintf(stderr,
            "ImageDiffusePerona: input and output image format must both be "
            "in float format.\n");
    exit(2);
  }

  if ((outImage->rows != inImage->rows) || (outImage->cols != inImage->cols)) {
    fprintf(stderr, "ImageDiffusePerona: input and output image sizes must match.\n");
    exit(2);
  }

  if (!ImageCheckSize(inImage, tmpImage, 0, 0, 0)) {
    if (tmpImage) ImageFree(&tmpImage);
    tmpImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  }
  if (!ImageCheckSize(inImage, gradImage, 0, 0, 0)) {
    if (gradImage) ImageFree(&gradImage);
    gradImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  }

  /* build index tables for border */
  xE = (int *)calloc((unsigned int)cols, sizeof(int));
  xW = (int *)calloc((unsigned int)cols, sizeof(int));
  yN = (int *)calloc((unsigned int)rows, sizeof(int));
  yS = (int *)calloc((unsigned int)rows, sizeof(int));

  xW[0] = 0;
  for (x = 1; x < cols; x++) {
    xW[x] = x - 1;
    xE[x - 1] = x;
  }
  xE[cols - 1] = cols - 1;
  yS[0] = 0;
  for (y = 1; y < rows; y++) {
    yN[y] = y - 1;
    yS[y - 1] = y;
  }
  yS[rows - 1] = rows - 1;

  ImageCopy(inImage, tmpImage);

#if 0
  if (0 && (Gdiag & DIAG_WRITE))
    fp = fopen("diffuse.dat", "w") ;
  else
#endif
  fp = NULL;

#if FOUR_CONNECTED
  fvals[1] = c[1] = 0.0f;
  fvals[3] = c[3] = 0.0f;
  fvals[6] = c[6] = 0.0f;
  fvals[8] = c[8] = 0.0f;
#endif
  for (i = 0; i < niter; i++) {
#if 0
    if (Gdiag & DIAG_WRITE)
      fprintf(stdout, "iteration %d\n", i) ;
#endif

    if (kimage) KernelImageCopy(kimage, ksrc);

    ImageSobel(tmpImage, gradImage, NULL, NULL);

    for (x = 0; x < cols; x++) {
      xe = xE[x];
      xw = xW[x];
      for (y = 0; y < rows; y++) {
        /*
          C1 |  C2  | C3
          ---------------
          C4 |  C0  | C5
          ---------------
          C6 |  C7  | C8
        */

        yn = yN[y];
        ys = yS[y];

        fvals[0] = *IMAGEFpix(tmpImage, x, y);
        fvals[2] = *IMAGEFpix(tmpImage, x, yn);
        fvals[4] = *IMAGEFpix(tmpImage, xw, y);
        fvals[5] = *IMAGEFpix(tmpImage, xe, y);
        fvals[7] = *IMAGEFpix(tmpImage, x, ys);
#if !FOUR_CONNECTED
        fvals[1] = *IMAGEFpix(tmpImage, xw, yn);
        fvals[3] = *IMAGEFpix(tmpImage, xe, yn);
        fvals[6] = *IMAGEFpix(tmpImage, xw, ys);
        fvals[8] = *IMAGEFpix(tmpImage, xe, ys);
#endif

        c[2] = KERNEL_MUL * C(*IMAGEFpix(gradImage, x, yn), k);
        c[4] = KERNEL_MUL * C(*IMAGEFpix(gradImage, xw, y), k);
        c[5] = KERNEL_MUL * C(*IMAGEFpix(gradImage, xe, y), k);
        c[7] = KERNEL_MUL * C(*IMAGEFpix(gradImage, x, ys), k);
#if !FOUR_CONNECTED
        c[1] = KERNEL_MUL * C(*IMAGEFpix(gradImage, xw, yn), k);
        c[3] = KERNEL_MUL * C(*IMAGEFpix(gradImage, xe, yn), k);
        c[6] = KERNEL_MUL * C(*IMAGEFpix(gradImage, xw, ys), k);
        c[8] = KERNEL_MUL * C(*IMAGEFpix(gradImage, xe, ys), k);
#endif

        c[0] = 1.0f;
        for (ci = 1; ci <= 8; ci++) c[0] -= c[ci];

        if (kimage) /* track evolution of kernel */
        {
          /* kernel stuff is row,col */
          KernelDiscount(kimage, y, x, c[0]);
          KernelUpdate(ksrc, kimage, y, x, yn, x, c[2]);
          KernelUpdate(ksrc, kimage, y, x, y, xw, c[4]);
          KernelUpdate(ksrc, kimage, y, x, y, xe, c[5]);
          KernelUpdate(ksrc, kimage, y, x, ys, x, c[7]);
#if !FOUR_CONNECTED
          KernelUpdate(ksrc, kimage, y, x, yn, xw, c[1]);
          KernelUpdate(ksrc, kimage, y, x, yn, xe, c[3]);
          KernelUpdate(ksrc, kimage, y, x, ys, xw, c[6]);
          KernelUpdate(ksrc, kimage, y, x, ys, xe, c[8]);
#endif
        }

        for (dst_val = 0.0f, ci = 0; ci <= 8; ci++) dst_val += fvals[ci] * c[ci];
        *IMAGEFpix(outImage, x, y) = dst_val;
      }
    }

    ImageCopy(outImage, tmpImage);
    k = k + k * slope;
    if (kimage) KernelImageNormalize(kimage);
  }

  free(xE);
  free(xW);
  free(yN);
  free(yS);

  if (fp) fclose(fp);

  return (0);
}
#else

int ImageDiffusePerona(IMAGE *inImage, IMAGE *outImage, double k, int niter, double slope, KIMAGE *kimage)
{
  UCHAR *csrc, *cdst;
  float *fsrc, *fdst, fcval, feval, fwval, fnval, fsval, fdval;
  int x, y, cval, eval, wval, nval, sval, i, dval;
  double deltaN, deltaE, deltaW, deltaS, gnval, gsval, geval, gwval;
  FILE *fp;
  static IMAGE *tmpImage = NULL;

  if ((outImage->pixel_format != inImage->pixel_format) || (outImage->pixel_format != PFFLOAT)) {
    fprintf(stderr, "ImageDiffuse: input and output image format must both be float.\n");
    exit(2);
  }

  if ((outImage->rows != inImage->rows) || (outImage->cols != inImage->cols)) {
    fprintf(stderr, "ImageDiffuse: input and output image sizes must match.\n");
    exit(2);
  }

  if (!ImageCheckSize(inImage, tmpImage, 0, 0, 0)) {
    if (tmpImage) ImageFree(&tmpImage);
    tmpImage = ImageAlloc(inImage->rows, inImage->cols, PFFLOAT, 1);
  }

  ImageCopy(inImage, tmpImage);

#if 0
  if (0 && (Gdiag & DIAG_WRITE))
    fp = fopen("diffuse.dat", "w") ;
  else
#endif
  fp = NULL;

  for (i = 0; i < niter; i++) {
    for (x = 0; x < tmpImage->cols; x++) {
      for (y = 0; y < tmpImage->rows; y++) {
        fsrc = IMAGEFpix(tmpImage, x, y);
        fcval = *fsrc;
        if (x > 0)
          fwval = *(fsrc - 1);
        else
          fwval = fcval;

        if (x < tmpImage->cols - 1)
          feval = *(fsrc + 1);
        else
          feval = fcval;

        if (y > 0)
          fnval = *IMAGEFpix(tmpImage, x, y - 1);
        else
          fnval = fcval;

        if (y < tmpImage->rows - 1)
          fsval = *IMAGEFpix(tmpImage, x, y + 1);
        else
          fsval = fcval;

        deltaN = (double)(fnval - fcval);
        deltaS = (double)(fsval - fcval);
        deltaE = (double)(feval - fcval);
        deltaW = (double)(fwval - fcval);

        gnval = g(deltaN, k);
        gsval = g(deltaS, k);
        geval = g(deltaE, k);
        gwval = g(deltaW, k);
        fdval = (float)(LAMBDA * (gnval * deltaN + gsval * deltaS + gwval * deltaW + geval * deltaE));

        if (fp) {
          fprintf(fp, "(%d, %d) = %2.2f, del = %2.2f\n", x, y, fcval, fdval);
          fprintf(fp, "N: (%2.2f) delta = %2.3lf, g = %2.3lf\n", fnval, deltaN, gnval);
          fprintf(fp, "S: (%2.2f) delta = %2.3lf, g = %2.3lf\n", fsval, deltaS, gsval);
          fprintf(fp, "W: (%2.2f) delta = %2.3lf, g = %2.3lf\n", fwval, deltaW, gwval);
          fprintf(fp, "E: (%2.2f) delta = %2.3lf, g = %2.3lf\n", feval, deltaE, geval);
          fprintf(fp, "val =  %2.2f + %2.2f = %2.2f\n\n", fcval, fdval, fcval + fdval);
        }

        fcval = fcval + fdval;
        *IMAGEFpix(outImage, x, y) = fcval;
      }
    }

    ImageCopy(outImage, tmpImage);
  }

  if (fp) fclose(fp);

  return (0);
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int ImageCurvature(IMAGE *inImage, float A, IMAGE *gradImage)
{
  static IMAGE *xImage = NULL, *yImage = NULL;
  int x, y, rows, cols;
  float *xpix, *ypix, *gradpix, xval, yval, gval, Asq;

  rows = inImage->rows;
  cols = inImage->cols;
  if (!ImageCheckSize(inImage, xImage, 0, 0, 0)) {
    if (xImage) ImageFree(&xImage);
    xImage = ImageAlloc(rows, cols, PFFLOAT, 1);
  }
  else {
    xImage->rows = rows;
    xImage->cols = cols;
  }
  if (!ImageCheckSize(inImage, yImage, 0, 0, 0)) {
    if (yImage) ImageFree(&yImage);
    yImage = ImageAlloc(rows, cols, PFFLOAT, 1);
  }
  else {
    yImage->rows = rows;
    yImage->cols = cols;
  }

  ImageConvolve3x3(inImage, sx, xImage);
  ImageConvolve3x3(inImage, sy, yImage);

  xpix = IMAGEFpix(xImage, 0, 0);
  ypix = IMAGEFpix(yImage, 0, 0);
  gradpix = IMAGEFpix(gradImage, 0, 0);
  Asq = A * A;
  for (y = 0; y < rows; y++) {
    for (x = 0; x < cols; x++) {
      xval = *xpix++;
      yval = *ypix++;
      gval = (float)sqrt(1.0 + (double)(Asq * (xval * xval + yval * yval)));
      *gradpix++ = gval;
    }
  }

  return (NO_ERROR);
}
