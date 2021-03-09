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

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include "diag.h"
#include "image.h"
#include "kernel.h"
#include "macros.h"
#include "proto.h"

/*----------------------------------------------------------------------
                           MACROS AND CONSTANTS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           LOCAL PROTOTYPES
----------------------------------------------------------------------*/
/*----------------------------------------------------------------------
                           FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
KIMAGE *KernelImageAlloc(int rows, int cols, int krows, int kcols)
{
  KIMAGE *kimage;
  int kernel_bytes, row, col;

  kernel_bytes = rows * cols * sizeof(KERNEL);
  kimage = (KIMAGE *)calloc(1, sizeof(KIMAGE) + kernel_bytes);
  if (!kimage) {
    fprintf(stderr, "KernelImageAlloc(%d, %d, %d, %d) - kimage failed\n", rows, cols, krows, kcols);
    exit(2);
  }
  kimage->krows = krows;
  kimage->kcols = kcols;
  kimage->rows = rows;
  kimage->cols = cols;
  kimage->fname = NULL;

  for (row = 0; row < rows; row++) {
    for (col = 0; col < cols; col++) KernelInit(kimage, row, col);
  }

  return (kimage);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
KIMAGE *KernelImageClone(KIMAGE *kimage)
{
  KIMAGE *knew;

  knew = KernelImageAlloc(kimage->rows, kimage->cols, kimage->krows, kimage->kcols);
  return (knew);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelImageCopy(KIMAGE *ksrc, KIMAGE *kdst)
{
  int row, col;

  for (row = 0; row < ksrc->rows; row++) {
    for (col = 0; col < ksrc->cols; col++) KernelCopy(ksrc, kdst, row, col, row, col);
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelCopy(KIMAGE *ksrc, KIMAGE *kdst, int src_row, int src_col, int dst_row, int dst_col)
{
  int row, col, cols;
  KERNEL *src_kernel, *dst_kernel;
  float *src_w, *dst_w;

  // if (src_row == 8 && dst_row == 8) src_col = src_row = 8; /* remove warning for now */

  src_kernel = KIMAGEpix(ksrc, dst_row, dst_col);
  dst_kernel = KIMAGEpix(kdst, dst_row, dst_col);
  cols = src_kernel->cols;
  for (row = 0; row < src_kernel->rows; row++) {
    src_w = src_kernel->weights[row];
    dst_w = dst_kernel->weights[row];
    for (col = 0; col < cols; col++) *dst_w++ = *src_w++;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelInit(KIMAGE *kimage, int row, int col)
{
  KERNEL *kernel;
  int halfrows, halfcols, wrow;

  kernel = KIMAGEpix(kimage, row, col);
  kernel->rows = kimage->krows;
  kernel->cols = kimage->kcols;

  /*
    center the kernel if possible, otherwise move kernel center so that
    it is as centered as possible around row,col without leaving the domain
    of the image.
  */
  halfrows = (kernel->rows) / 2;
  halfcols = (kernel->cols) / 2;
#if 0
  if (row-halfrows < 0)                             /* off top of image */
    kernel->row0 = 0 ;
  else if (row+halfrows >= kimage->rows)          /* off bottom of image */
    kernel->row0 = kimage->rows - kernel->rows ;
  else                                              /* okay - center kernel */
#endif
  kernel->row0 = row - halfrows;

#if 0
  if (col-halfcols < 0)                             /* off left of image */
    kernel->col0 = 0 ;
  else if (col+halfcols >= kimage->cols)           /* off right of image */
    kernel->col0 = kimage->cols - kernel->cols ;
  else                                              /* okay - center kernel */
#endif
  kernel->col0 = col - halfcols;

  /* allocate and initialize weights */
  kernel->weights = (float **)calloc(kernel->rows, sizeof(float *));
  if (!kernel->weights) {
    fprintf(stderr, "KernelInit(%d, %d): could not allocate weights\n", row, col);
    exit(2);
  }

  for (wrow = 0; wrow < kernel->rows; wrow++) {
    kernel->weights[wrow] = (float *)calloc(kernel->cols, sizeof(float));
    if (!kernel->weights[wrow]) {
      fprintf(stderr, "KernelInit: could not allocate %dth row\n", wrow);
      exit(3);
    }
  }

  kernel->weights[row - kernel->row0][col - kernel->col0] = 1.0f;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelImageFree(KIMAGE *kimage)
{
  int row, col;

  for (row = 0; row < kimage->rows; row++) {
    for (col = 0; col < kimage->cols; col++) {
      KernelFree(KIMAGEpix(kimage, row, col));
    }
  }
  free(kimage);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              free all the contents of the kernel, but not the kernel
              itself
----------------------------------------------------------------------*/
void KernelFree(KERNEL *kernel)
{
  int wrow;

  for (wrow = 0; wrow < kernel->rows; wrow++) free(kernel->weights[wrow]);

  free(kernel->weights);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelImageDump(KIMAGE *kimage, FILE *fp)
{
  int row, col, krow, kcol;
  KERNEL *kernel;

  fprintf(fp, "image size: %d x %d, kernel size: %d x %d\n", kimage->rows, kimage->cols, kimage->krows, kimage->kcols);

  for (row = 0; row < kimage->rows; row++) {
    for (col = 0; col < kimage->cols; col++) {
      kernel = KIMAGEpix(kimage, row, col);
      fprintf(fp, "kernel at (%d, %d), row0, col0 = (%d, %d)\n", row, col, kernel->row0, kernel->col0);
      for (krow = 0; krow < kernel->rows; krow++) {
        for (kcol = 0; kcol < kernel->cols; kcol++) {
          if (fabs(kernel->weights[krow][kcol]) > 0.00001) {
            fprintf(fp,
                    "\t(%d, %d) --> (%d, %d) = %f\n",
                    krow,
                    kcol,
                    krow + kernel->row0,
                    kcol + kernel->col0,
                    kernel->weights[krow][kcol]);
          }
        }
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelDiscount(KIMAGE *kimage, int row, int col, float weight)
{
  KERNEL *kernel;
  float *w;
  int cols;

  DiagPrintf(DIAG_KERNELS, "KernelDiscount(%d, %d, %f)\n", row, col, weight);

  kernel = KIMAGEpix(kimage, row, col);
  cols = kernel->cols;
  for (row = 0; row < kernel->rows; row++) {
    w = kernel->weights[row];
    for (col = 0; col < cols; col++) *w++ *= weight;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelUpdate(KIMAGE *ksrc, KIMAGE *kdst, int dst_row, int dst_col, int src_row, int src_col, float weight)
{
  KERNEL *src_kernel, *dst_kernel;
  float *w;
  int row, col, src_rows, src_cols, dst_rows, dst_cols, dst_row0, dst_col0, src_row0, src_col0, drow, dcol;

  DiagPrintf(DIAG_KERNELS, "KernelUpdate(%d, %d, %d, %d, %f)\n", dst_row, dst_col, src_row, src_col, weight);

  src_kernel = KIMAGEpix(ksrc, src_row, src_col);
  dst_kernel = KIMAGEpix(kdst, dst_row, dst_col);

  src_rows = src_kernel->rows;
  src_cols = src_kernel->cols;
  src_row0 = src_kernel->row0;
  src_col0 = src_kernel->col0;

  dst_rows = dst_kernel->rows;
  dst_cols = dst_kernel->cols;
  dst_row0 = dst_kernel->row0;
  dst_col0 = dst_kernel->col0;

  /*
    go through each point in the source kernel. If it maps to a point that
    is within the domain of the destination kernel, then add its weight
    to the corresponding point in the destination kernel.
  */
  drow = src_row0 - dst_row0;
  dcol = src_col0 - dst_col0;
  for (row = 0; row < src_rows; row++) {
    w = src_kernel->weights[row];
    for (col = 0; col < src_cols; col++, w++) {
      if (src_kernel->weights[row][col] > 0.0001) weight = weight * 1.0f;

      dst_row = row + drow;
      dst_col = col + dcol;
      if ((dst_row < 0) || (dst_col < 0) || (dst_row >= dst_rows) || (dst_col >= dst_cols)) continue; /* out of range */
      dst_kernel->weights[dst_row][dst_col] += *w * weight;
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelImageNormalize(KIMAGE *kimage)
{
  int row, col;

  for (row = 0; row < kimage->rows; row++) {
    for (col = 0; col < kimage->cols; col++) KernelNormalize(kimage, row, col);
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelNormalize(KIMAGE *kimage, int row, int col)
{
  KERNEL *kernel;
  float *w, total;
  int krow, kcol, cols, krows, kcols;

  kernel = KIMAGEpix(kimage, row, col);
  cols = kernel->cols;

  /* set kernel locations outside of image to 0 */
  krows = -kernel->row0;
  for (krow = 0; krow < krows; krow++) /* erase 1st krows */
  {
    w = kernel->weights[krow];
    memset((char *)w, 0, kernel->cols * sizeof(float));
  }
  kcols = -kernel->col0;
  if (kcols > 0) /* erase 1st kcols */
  {
    for (krow = 0; krow < kernel->rows; krow++) {
      w = kernel->weights[krow];
      memset((char *)w, 0, kcols * sizeof(float));
    }
  }
  kcols = kernel->col0 + kernel->cols - kimage->cols;
  if (kcols > 0) /* erase last kcols */
  {
    for (krow = 0; krow < kernel->rows; krow++) {
      w = kernel->weights[krow];
      memset((char *)(w + kernel->cols - kcols), 0, kcols * sizeof(float));
    }
  }
  krows = kernel->row0 + kernel->rows - kimage->rows;
  kcols = kernel->cols;
  for (krow = kernel->rows - krows; krow < kernel->rows; krow++) { /* erase last krows */
    w = kernel->weights[krow];
    memset((char *)w, 0, kcols * sizeof(float));
  }

#if 0
  /* make kernel positive */
  for (kmin = 100000.0f, krow = 0 ; krow < kernel->rows ; krow++)
  {
    w = kernel->weights[krow] ;
    for (kcol = 0 ; kcol < cols ; kcol++, w++)
    {
      if (*w < kmin)
        kmin = *w ;
    }
  }
  for (krow = 0 ; krow < kernel->rows ; krow++)
  {
    w = kernel->weights[krow] ;
    for (kcol = 0 ; kcol < cols ; kcol++)
      *w++ -= kmin ;
  }
#endif

  /* now normalize total weights to be 1 */
  for (total = 0.0f, krow = 0; krow < kernel->rows; krow++) {
    w = kernel->weights[krow];
    for (kcol = 0; kcol < cols; kcol++) {
      total += (float)fabs(*w++);
    }
  }

  if (FZERO(total)) /* shouldn't happen */
  {
    fprintf(stderr, "KernelNormalize(%d, %d): zero total!\n", row, col);
    return;
  }

  for (krow = 0; krow < kernel->rows; krow++) {
    w = kernel->weights[krow];
    for (col = 0; col < cols; col++) *w++ /= total;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelImageConvolve(KIMAGE *kimage, IMAGE *src_image, IMAGE *dst_image)
{
  int src_row, src_rows, src_cols, start_col, start_row;
  int dst_row, dst_col, dst_rows, dst_cols;
  int krow, kcol, krows, kcols, row0, col0;
  KERNEL *kernel;
  float *src_pix, *dst_pix;
  float *w, total;

  src_rows = src_image->rows;
  src_cols = src_image->cols;
  dst_rows = dst_image->rows;
  dst_cols = dst_image->cols;

  dst_pix = IMAGEFpix(dst_image, 0, 0);
  for (dst_row = 0; dst_row < dst_rows; dst_row++) {
    for (dst_col = 0; dst_col < dst_cols; dst_col++) {
      /* calculate kernel value at this point */
      kernel = KIMAGEpix(kimage, dst_row, dst_col);
      row0 = kernel->row0;
      col0 = kernel->col0;

      /*
        if the kernel extends off of an edge of the image, we want to ignore
        those weights, making sure that the image and kernel are still properly
        aligned at the new starting point. If the kernel extends off of the
        end of the image, then clip the effective number of rows and cols in
        the kernel to reflect the overlap with the image.
      */
      krows = kimage->krows;
      kcols = kimage->kcols;
      if (col0 + kcols >= src_cols) /* kernel extends off bottom of image */
        kcols = src_cols - col0;

      if (row0 + krows >= src_rows) /* kernel extends off right of image */
        krows = src_rows - row0;

      if (col0 < 0) /* kernel extends off of top edge of image */
      {
        start_col = -col0;
        col0 = 0;
      }
      else
        start_col = 0;

      if (row0 < 0) /* kernel extends off of left edge of image */
      {
        start_row = -row0;
        row0 = 0;
      }
      else
        start_row = 0;

      total = 0.0f;

      src_row = row0;
      for (krow = start_row; krow < krows; krow++, src_row++) {
        src_pix = IMAGEFpix(src_image, col0, src_row);
        w = kernel->weights[krow] + start_col;

        for (kcol = start_col; kcol < kcols; kcol++) total += *w++ * *src_pix++;
      }
      *dst_pix++ = total;
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void KernelImageWrite(KIMAGE *kimage, char *fname, int argc, char *argv[])
{
  IMAGE *image;
  int i;
  char str[100];

  image = KernelImageToSeq(kimage);
  for (i = 0; i < argc; i++) {
    argv[0] = argv[i];
    update_header(image, 1, argv);
  }
  if (kimage->fname) {
    free(image->orig_name);
    image->orig_name = STRCPALLOC(kimage->fname);
  }
  sprintf(str, "%d %d", kimage->rows, kimage->cols);
  image->seq_name = STRCPALLOC(str);
  ImageWrite(image, fname);
  ImageFree(&image);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
KIMAGE *KernelImageRead(char *fname)
{
  KIMAGE *kimage;
  IMAGE *image;

  image = ImageRead(fname);
  if (!image) return (NULL);

  kimage = KernelImageFromSeq(image);

  ImageFree(&image);
  return (kimage);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *KernelImageToSeq(KIMAGE *kimage)
{
  IMAGE *image;
  int num_frame, row, col, rows, cols, krow, kcol, krows, kcols, pix_per_frame, row0, col0, drow, dcol, npix, col_dif,
      row_dif;
  KERNEL *kernel;
  float *fsrc, *fdst, *fbase;

  krows = kimage->krows;
  kcols = kimage->kcols;
  rows = kimage->rows;
  cols = kimage->cols;
  pix_per_frame = krows * kcols;
  num_frame = rows * cols;

  image = ImageAlloc(krows, kcols, PFFLOAT, num_frame);
  if (!image) {
    fprintf(stderr, "KernelImageToSeq: could not allocate sequence\n");
    return (NULL);
  }

  for (row = 0; row < rows; row++) {
    for (col = 0; col < cols; col++) {
      kernel = KIMAGEpix(kimage, row, col);
      row0 = kernel->row0;
      col0 = kernel->col0;
      fbase = IMAGEFpix(image, 0, 0) + ((row * cols) + col) * pix_per_frame;
      memset((char *)fbase, 0, pix_per_frame * sizeof(float));

      /*
            copy kernel into Image frame, changing coordinate systems so that
            the kernel is centered in the frame. This involves truncating some of
            the kernel.
      */
      for (krow = 0; krow < krows; krow++) {
        row_dif = (row0 + krows / 2) - row;
        col_dif = (col0 + kcols / 2) - col;

        /*
                calculate destination location. The point (row,col) should map
                to the center of the Image image.
        */
        drow = krow + row_dif;
        dcol = col_dif;
        if (dcol < 0) /* destination before start of image */
        {
          kcol = -dcol;
          dcol = 0;
        }
        else
          kcol = 0;

        npix = kcols - abs(col_dif);

        fsrc = kernel->weights[krow] + kcol;
        fdst = fbase + drow * kcols + dcol;
        if (drow >= 0 && (drow < krows) && (npix > 0)) memmove((char *)fdst, (char *)fsrc, npix * sizeof(float));
      }
    }
  }

  return (image);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
KIMAGE *KernelImageFromSeq(IMAGE *image)
{
  int rows, cols, krows, kcols, row, col, krow, kcol, pix_per_frame, row0, col0, drow, dcol, npix, col_dif, row_dif;
  KERNEL *kernel;
  float *fsrc, *fdst, *fbase;
  KIMAGE *kimage;

  rows = nint(sqrt(image->num_frame));
  cols = rows;

  if (image->seq_name) /* rows and cols may be stored in seq_name field */
    sscanf(image->seq_name, "rows=%d cols=%d", &rows, &cols);

  krows = image->rows;
  kcols = image->cols;
  pix_per_frame = krows * kcols;

  kimage = KernelImageAlloc(rows, cols, krows, kcols);
  if (!kimage) {
    fprintf(stderr, "KernelImageFromSeq: could not allocate kimage\n");
    return (NULL);
  }

  if (image->orig_name) kimage->fname = STRCPALLOC(image->orig_name);
  for (row = 0; row < rows; row++) {
    for (col = 0; col < cols; col++) {
      kernel = KIMAGEpix(kimage, row, col);
      row0 = kernel->row0;
      col0 = kernel->col0;
      fbase = IMAGEFpix(image, 0, 0) + ((row * cols) + col) * pix_per_frame;

      /*
            copy Image frame into kernel, changing coordinate systems so that
            the kernel is centered in the frame. This involves truncating some of
            the kernel.
      */
      for (krow = 0; krow < krows; krow++) {
        row_dif = (row0 + krows / 2) - row;
        col_dif = (col0 + kcols / 2) - col;

        /*
                calculate destination location. The point (row,col) should map
                to the center of the hips image.
        */
        drow = krow + row_dif;
        dcol = col_dif;
        if (dcol < 0) /* destination before start of image */
        {
          kcol = -dcol;
          dcol = 0;
        }
        else
          kcol = 0;

        npix = kcols - abs(col_dif);

        fdst = kernel->weights[krow] + kcol;
        fsrc = fbase + drow * kcols + dcol;
        if (drow >= 0 && (drow < krows) && (npix > 0)) memmove((char *)fdst, (char *)fsrc, npix * sizeof(float));
      }
    }
  }

  return (kimage);
}
