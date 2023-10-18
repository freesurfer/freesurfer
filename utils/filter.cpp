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

/*
 *       FILE NAME:   filter.c
 *
 *       DESCRIPTION: image processing filters
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        6/18/96
 *
 */

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hips.h"
#include "diag.h"
#include "error.h"
#include "filter.h"
#include "image.h"
#include "machine.h"
#include "macros.h"
#include "proto.h"
#include "utils.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define NITSHI_TAU 2.0f

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/
/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

#define DEBUG_FILTER 0
IMAGE *Ifilter = NULL;
#if DEBUG_FILTER
extern int Gx, Gy;
#endif

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageNitShiFilter(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy, int wsize, double sigma, IMAGE *Idst)
{
  static IMAGE *IE = NULL, *IF = NULL, *IG = NULL, *Iblur = NULL;
  IMAGE *Iin, *Iout;
  int rows, cols, x, y, whalf, xk, yk, xs, ys;
  float norm, total, *dpix, *spix, fmin, fmax, fval, Eval, Fval, Gval, sigma_sq;

#if DEBUG_FILTER
  if (Ifilter) ImageFree(&Ifilter);
  Ifilter = ImageAlloc(wsize, wsize, PFFLOAT, 1);
#endif

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (Isrc->pixel_format != PFFLOAT) {
    Iin = ImageAlloc(rows, cols, PFFLOAT, 1);
    ImageCopy(Isrc, Iin);
    ImageValRange(Isrc, &fmin, &fmax);
  }
  else
    Iin = Isrc;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else
    Iout = Idst;

  if (!ImageCheckSize(Iin, IE, 0, 0, 0)) {
    if (IE) {
      ImageFree(&IE);
      ImageFree(&IF);
      ImageFree(&IG);
    }
    IE = ImageAlloc(rows, cols, PFFLOAT, 1);
    IF = ImageAlloc(rows, cols, PFFLOAT, 1);
    IG = ImageAlloc(rows, cols, PFFLOAT, 1);
  }

  if (!Iblur) Iblur = ImageGaussian1d(NITSHI_TAU, 0);
  whalf = (wsize - 1) / 2;

  /* calculate IE, IF, and IG */
  ImageMul(Ix, Ix, IE);
  ImageMul(Ix, Iy, IF);
  ImageMul(Iy, Iy, IG);
  ImageConvolveGaussian(IE, Iblur, IE, 0);
  ImageConvolveGaussian(IF, Iblur, IF, 0);
  ImageConvolveGaussian(IG, Iblur, IG, 0);

#if 0
  ImageWrite(Ix, "Ix.hipl") ;
  ImageWrite(Iy, "Iy.hipl") ;
  ImageWrite(IE, "IE.hipl");
  ImageWrite(IF, "IF.hipl");
  ImageWrite(IG, "IG.hipl");
#endif

  /* now apply actual filter */
  sigma_sq = 2.0f * (float)(sigma * sigma);
  dpix = IMAGEFpix(Iout, 0, 0);
  for (y = 0; y < rows; y++) {
    for (x = 0; x < cols; x++) {
      norm = total = 0.0f;

      /* apply kernel */
      for (yk = -whalf; yk <= whalf; yk++) {
        ys = y + yk;
        if (ys < 0)
          ys = 0;
        else if (ys >= rows)
          ys = rows - 1;

        for (xk = -whalf; xk <= whalf; xk++) {
          xs = x + xk;
          if (xs < 0)
            xs = 0;
          else if (xs >= cols)
            xs = cols - 1;

          spix = IMAGEFpix(Iin, xs, ys);
          Eval = *IMAGEFpix(IE, xs, ys);
          Fval = *IMAGEFpix(IF, xs, ys);
          Gval = *IMAGEFpix(IG, xs, ys);

          fval = (float)exp(-(Eval * xk * xk + 2 * Fval * xk * yk + Gval * yk * yk) / sigma_sq);
#if DEBUG_FILTER
          if (x == Gx && y == Gy) *IMAGEFpix(Ifilter, xk + whalf, yk + whalf) = fval;
#endif
          total += fval * *spix;
          norm += fval;
        }
      }

#if DEBUG_FILTER
      if (x == Gx && y == Gy)
        for (yk = -whalf; yk <= whalf; yk++) {
          for (xk = -whalf; xk <= whalf; xk++) {
            *IMAGEFpix(Ifilter, xk + whalf, yk + whalf) /= norm;
          }
        }
#endif
      *dpix++ = total / norm;
    }
  }

  if (Iin != Isrc) ImageFree(&Iin);

  if (Iout != Idst) {
    /* should worry about using fmin and fmax to rescale byte images here... */
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int ImageBuildExponentialFilter(IMAGE *gradImage, int wsize, float k, IMAGE *offsetImage, IMAGE *filterSequence)
{
  int rows, cols, x, y, whalf, xc, yc, x0, y0, dx, dy, frame;
  float fpix, *g, norm, val, *filterPix;
  static float *gaussian = NULL;
  static int w = 0;

  rows = gradImage->rows;
  cols = gradImage->cols;

  whalf = (wsize - 1) / 2;

  if (wsize != w) {
    w = wsize;
    free(gaussian);
    gaussian = NULL;
  }

  if (!gaussian) /* allocate a gaussian bump */
  {
    float den;

    gaussian = (float *)calloc(wsize * wsize, sizeof(float));
    den = wsize * wsize + wsize + 1;
    norm = 0.0f;
    for (g = gaussian, y = 0; y < wsize; y++) {
      yc = y - whalf;
      for (x = 0; x < wsize; x++, g++) {
        xc = x - whalf;
        *g = (float)exp(-36.0 * sqrt((double)(xc * xc + yc * yc)) / (double)den);
        norm += *g;
      }
    }

    /* normalize gaussian */
    for (g = gaussian, y = 0; y < wsize; y++) {
      for (x = 0; x < wsize; x++, g++) *g /= norm;
    }
  }

  /*
    x and y are in window coordinates, while xc and yc are in image
    coordinates.
  */
  for (frame = y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, frame++) {
      if (offsetImage) {
        dx = nint(*IMAGEFpix(offsetImage, x0, y0));
        dy = nint(*IMAGEIseq_pix(offsetImage, x0, y0, 1));
      }
      else
        dx = dy = 0;

      norm = 0.0f;
      filterPix = IMAGEFseq_pix(filterSequence, 0, 0, frame);

      if (x0 == 5 && y0 == 10) x0 = 70;

      for (g = gaussian, y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0 + dy;
        if (yc < 0)
          yc = -yc;
        else if (yc >= rows)
          yc = rows - (yc - rows + 1);

        for (x = -whalf; x <= whalf; x++) {
          xc = x0 + x + dx;
          if (xc < 0)
            xc = -xc;
          else if (xc >= cols)
            xc = cols - (xc - cols + 1);

          fpix = *IMAGEFpix(gradImage, xc, yc);
          val = (float)exp((double)(-fpix * fpix / k)) /*  * *g++ */;
          norm += val;
          *filterPix++ = val;
        }
      }

      if (FZERO(norm)) continue;

      /* normalize kernel weights to sum to 1 */
      filterPix = IMAGEFseq_pix(filterSequence, 0, 0, frame);
      for (y = 0; y < wsize; y++) {
        for (x = 0; x < wsize; x++) *filterPix++ /= norm;
      }
    }
  }

  /*  ImageWrite(filterSequence, "filter.hipl") ;*/
  return (0);
}
#if 0
/*----------------------------------------------------------------------
            Parameters:
           Description:
----------------------------------------------------------------------*/
int ImageCalculateMomentOffset(IMAGE *gradImage, int wsize, float c,
                               IMAGE *offsetImage)
{
  return(0) ;
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int ImageSpaceVariantFilter(IMAGE *inImage, IMAGE *filterSequence, IMAGE *outImage) { return (0); }
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int ImageExponentialFilter(IMAGE *inImage, IMAGE *gradImage, int wsize, float k, IMAGE *offsetImage, IMAGE *outImage)
{
  int rows, cols, x, y, whalf, xc, yc, x0, y0, dx, dy;
  float fpix, *g, norm, val, *filterPix, *filter, *outPix;
  static float w, *gaussian;

  xc = yc = 0; /* eliminate compiler warning */
  filter = (float *)calloc(wsize * wsize, sizeof(float));
  if (!filter) ErrorReturn(ERROR_NO_MEMORY, (ERROR_NO_MEMORY, "ImageExponentialFilter: could not allocate filter"));

  rows = gradImage->rows;
  cols = gradImage->cols;

  whalf = (wsize - 1) / 2;

  if (wsize != w) {
    w = wsize;
    free(gaussian);
    gaussian = NULL;
  }

  if (!gaussian) /* allocate a gaussian bump */
  {
    float den;

    gaussian = (float *)calloc(wsize * wsize, sizeof(float));
    den = wsize * wsize + wsize + 1;
    norm = 0.0f;
    for (g = gaussian, y = 0; y < wsize; y++) {
      yc = y - whalf;
      for (x = 0; x < wsize; x++, g++) {
        xc = x - whalf;
        *g = (float)exp(-36.0 * sqrt((double)(xc * xc + yc * yc)) / (double)den);
        norm += *g;
      }
    }

    /* normalize gaussian */
    for (g = gaussian, y = 0; y < wsize; y++) {
      for (x = 0; x < wsize; x++, g++) *g /= norm;
    }
  }

  /*
    x and y are in window coordinates, while xc and yc are in image
    coordinates.
  */
  outPix = IMAGEFpix(outImage, 0, 0);
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++) {
      if (offsetImage) {
        dx = nint(*IMAGEFpix(offsetImage, x0, y0));
        dy = nint(*IMAGEFseq_pix(offsetImage, x0, y0, 1));
      }
      else
        dx = dy = 0;

      norm = 0.0f;
      filterPix = filter;

      for (g = gaussian, y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0 + dy;
        if (yc < 0)
          yc = 0;
        else if (yc >= rows)
          yc = rows - 1;
        for (x = -whalf; x <= whalf; x++) {
          xc = x0 + x + dx;
          if (xc < 0)
            xc = 0;
          else if (xc >= cols)
            xc = cols - 1;

          fpix = *IMAGEFpix(gradImage, xc, yc);
          val = (float)exp((double)(-fpix * fpix / k)) /*  * *g++ */;
          norm += val;
          *filterPix++ = val;
        }
      }

      if (FZERO(norm)) /* neigborhood is all zeros */
      {
        *outPix++ = 0.0f;
        continue;
      }

      /* normalize kernel weights to sum to 1 */
      filterPix = filter;
      for (y = 0; y < wsize; y++) {
        for (x = 0; x < wsize; x++) *filterPix++ /= norm;
      }

      /*
            now apply filter to this point in the image, taking possible
            offset into accound
      */
      filterPix = filter;
      val = 0.0f;
      for (y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0 + dy;
        if (yc < 0)
          yc = 0;
        else if (yc >= rows)
          yc = rows - 1;

        for (x = -whalf; x <= whalf; x++) {
          xc = x0 + x + dx;
          if (xc < 0)
            xc = 0;
          else if (xc >= cols)
            xc = cols - 1;

          fpix = *IMAGEFpix(inImage, xc, yc);
          val += fpix * *filterPix++;
        }
      }

      if (isnan(val)) {
        fprintf(stderr, "(%d, %d) = NaN!\n", xc, yc);
      }
      *outPix++ = val;
    }
  }

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int compare_sort_array(const void *pf1, const void *pf2);
int ImageMedianFilter(IMAGE *inImage, int wsize, IMAGE *offsetImage, IMAGE *outImage)
{
  static float *sort_array = NULL;
  static int sort_size = 0;
  int ecode, x0, y0, rows, cols, x, y, whalf, yc, dx, dy, frame, wsq, median_index;
  float *sptr, *outPix, min_val, max_val, *inPix;
  ubyte *in_image, *out_image;
  IMAGE *Iout, *Iin;
  int xc;

  rows = inImage->rows;
  cols = inImage->cols;
  if (!offsetImage) {
    /* h_median only takes byte formatted input and output images */
    if (inImage->pixel_format != PFBYTE) {
      Iin = ImageAlloc(rows, cols, PFBYTE, inImage->num_frame);
      ImageValRange(inImage, &min_val, &max_val);
      ImageScale(inImage, inImage, 0.0f, 255.0f);
      ImageCopy(inImage, Iin);
      ImageScale(inImage, inImage, min_val, max_val); /* restore old image */
    }
    else
      Iin = inImage;

    if (inImage->pixel_format != PFBYTE)
      Iout = ImageAlloc(rows, cols, PFBYTE, outImage->num_frame);
    else
      Iout = outImage;

    out_image = Iout->image;
    in_image = Iin->image;
    for (frame = 0; frame < inImage->num_frame; frame++) {
      ecode = h_median(Iin, Iout, wsize);
      if (ecode != HIPS_OK) ErrorReturn(-1, (ecode, "ImageMedian: h_median failed (%d)", ecode));

      Iout->firstpix += Iout->sizeimage;
      Iout->image += Iout->sizeimage;

      Iin->firstpix += Iin->sizeimage;
      Iin->image += Iin->sizeimage;
    }

    Iout->firstpix = Iout->image = out_image;
    Iin->firstpix = Iin->image = in_image;
    if (Iin != inImage) ImageFree(&Iin);
    if (outImage != Iout) {
      ImageCopy(Iout, outImage);
      ImageScale(outImage, outImage, min_val, max_val);
      ImageFree(&Iout);
    }
    return (NO_ERROR);
  }

  median_index = wsize * wsize / 2;
  wsq = wsize * wsize;
  whalf = (wsize - 1) / 2;

  /* create a static array for sorting pixels in */
  if (wsize > sort_size) {
    sort_size = wsize;
    if (sort_array) sort_array = NULL;
  }

  if (!sort_array) sort_array = (float *)calloc(wsq, sizeof(float));

  for (frame = 0; frame < inImage->num_frame; frame++) {
    outPix = IMAGEFseq_pix(outImage, 0, 0, frame);
    for (y0 = 0; y0 < rows; y0++) {
      for (x0 = 0; x0 < cols; x0++) {
        /*
                 x and y are in window coordinates, while xc and yc are in image
                 coordinates.
         */
        if (offsetImage) {
          dx = nint(*IMAGEFpix(offsetImage, x0, y0));
          dy = nint(*IMAGEFseq_pix(offsetImage, x0, y0, 1));
        }
        else
          dx = dy = 0;

        for (sptr = sort_array, y = -whalf; y <= whalf; y++) {
          /* reflect across the boundary */
          yc = y + y0 + dy;
#if 0
          if (yc < 0)
            yc = -yc ;
          else if (yc >= rows)
            yc = rows - (yc - rows + 1) ;
#else
          if (yc < 0)
            yc = 0;
          else if (yc >= rows)
            yc = rows - 1;
#endif

          inPix = IMAGEFseq_pix(inImage, 0, yc, frame);
          for (x = -whalf; x <= whalf; x++) {
            xc = x0 + x + dx;
#if 0
            if (xc < 0)
              xc = -xc ;
            else if (xc >= cols)
              xc = cols - (xc - cols + 1) ;
#else
            if (xc < 0)
              xc = 0;
            else if (xc >= cols)
              xc = cols - 1;
#endif

#if 0
            *sptr++ = *IMAGEFseq_pix(inImage, xc, yc, frame) ;
#else
            *sptr++ = *(inPix + xc);
#endif
          }
        }
        qsort(sort_array, wsq, sizeof(float), compare_sort_array);
        *outPix++ = sort_array[median_index];
      }
    }
  }
  return (0);
}
static int compare_sort_array(const void *pf1, const void *pf2)
{
  float f1, f2;

  f1 = *(float *)pf1;
  f2 = *(float *)pf2;

  /*  return(f1 > f2 ? 1 : f1 == f2 ? 0 : -1) ;*/
  if (f1 > f2)
    return (1);
  else if (f1 < f2)
    return (-1);

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageConvolveGaussianByte(IMAGE *Isrc, IMAGE *gImage, IMAGE *Iout, int dst_frameno)
{
  static IMAGE *Itmp = NULL;
  int ksize;
  float *kernel;
  ubyte *buf;

  if (!ImageCheckSize(Isrc, Itmp, 0, 0, 0)) {
    if (Itmp) ImageFree(&Itmp);
    Itmp = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, 1);
  }
  if (!Iout) Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, 1);

  ImageSetSize(Itmp, Isrc->rows, Isrc->cols);

  kernel = IMAGEFpix(gImage, 0, 0);
  ksize = gImage->cols;
  ImageConvolve1dByte(Isrc, Itmp, kernel, ksize, IMAGE_VERTICAL);

  buf = IMAGEpix(Iout, 0, 0);
  Iout->image = IMAGEseq_pix(Iout, 0, 0, dst_frameno);
  ImageConvolve1dByte(Itmp, Iout, kernel, ksize, IMAGE_HORIZONTAL);

  Iout->image = buf;
  return (Iout);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageConvolveGaussian(IMAGE *Isrc, IMAGE *gImage, IMAGE *Iout, int dst_frameno)
{
  static IMAGE *Itmp = NULL;
  int ksize;
  float *kernel, *buf;

  if (Isrc->pixel_format == PFBYTE) return (ImageConvolveGaussianByte(Isrc, gImage, Iout, dst_frameno));

  if (!ImageCheckSize(Isrc, Itmp, 0, 0, 0)) {
    if (Itmp) ImageFree(&Itmp);
    Itmp = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);
  }
  if (!Iout) Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);

  ImageSetSize(Itmp, Isrc->rows, Isrc->cols);

  kernel = IMAGEFpix(gImage, 0, 0);
  ksize = gImage->cols;
  ImageConvolve1d(Isrc, Itmp, kernel, ksize, IMAGE_VERTICAL);

  buf = IMAGEFpix(Iout, 0, 0);
  Iout->image = (ubyte *)IMAGEFseq_pix(Iout, 0, 0, dst_frameno);
  ImageConvolve1d(Itmp, Iout, kernel, ksize, IMAGE_HORIZONTAL);

  Iout->image = (ubyte *)buf;
  return (Iout);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageCircularConvolveGaussian(IMAGE *Isrc, IMAGE *gImage, IMAGE *Iout, int dst_frameno)
{
  static IMAGE *Itmp = NULL;
  int ksize;
  float *kernel, *buf;

  if (Isrc->pixel_format != PFFLOAT)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "ImageCircularConvolveGaussian: type %d"
                 "not supported",
                 Isrc->pixel_format));

  if (!ImageCheckSize(Isrc, Itmp, 0, 0, 0)) {
    if (Itmp) ImageFree(&Itmp);
    Itmp = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);
  }
  if (!Iout) Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);

  ImageSetSize(Itmp, Isrc->rows, Isrc->cols);

  kernel = IMAGEFpix(gImage, 0, 0);
  ksize = gImage->cols;
  ImageCircularConvolve1d(Isrc, Itmp, kernel, ksize, IMAGE_VERTICAL);

  buf = IMAGEFpix(Iout, 0, 0);
  Iout->image = (ubyte *)IMAGEFseq_pix(Iout, 0, 0, dst_frameno);
  ImageCircularConvolve1d(Itmp, Iout, kernel, ksize, IMAGE_HORIZONTAL);

  Iout->image = (ubyte *)buf;
  return (Iout);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageConvolveGaussianFrames(IMAGE *Isrc, IMAGE *gImage, IMAGE *Idst)
{
  int frame, src_frames, dst_frames;
  ubyte *src_buf, *dst_buf;

  if (!Idst) Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);

  src_frames = Isrc->num_frame;
  dst_frames = Idst->num_frame;
  Isrc->num_frame = Idst->num_frame = 1;
  src_buf = Isrc->image;
  dst_buf = Idst->image;
  for (frame = 0; frame < src_frames; frame++) {
    ImageConvolveGaussian(Isrc, gImage, Idst, 0);
    Isrc->image += Isrc->sizeimage;
    if (Isrc != Idst) Idst->image += Idst->sizeimage;
  }

  Isrc->image = src_buf;
  Idst->image = dst_buf;
  Isrc->num_frame = src_frames;
  Idst->num_frame = dst_frames;

  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void ImageConvolve1d(IMAGE *I, IMAGE *J, float k[], int len, int axis)
{
  int x, y, width, height, halflen;
  int i;
  float *outPix;
  float *ki, total, *inBase;
  static int *xi_LUT = NULL, LUT_width, LUT_height, *yi_LUT = NULL, LUT_len = 0;

  width = I->cols;
  height = I->rows;

  halflen = len / 2;

  if ((LUT_len != len) || (LUT_width != width)) {
    LUT_width = width;
    if (xi_LUT) free(xi_LUT - LUT_len / 2);
    xi_LUT = (int *)calloc(width + 2 * halflen, sizeof(int));
    if (!xi_LUT) ErrorExit(ERROR_NO_MEMORY, "ImageConvolve1d: could not allocate LUT\n");
    xi_LUT += halflen;
    for (i = -halflen; i < width + halflen; i++) {
      if (i < 0)
        xi_LUT[i] = 0;
      else if (i >= width)
        xi_LUT[i] = width - 1;
      else
        xi_LUT[i] = i;
    }
  }

  if ((LUT_len != len) || (LUT_height != height)) {
    LUT_height = height;
    if (yi_LUT) free(yi_LUT - LUT_len / 2);
    yi_LUT = (int *)calloc(height + 2 * halflen, sizeof(int));
    if (!yi_LUT) ErrorExit(ERROR_NO_MEMORY, "ImageConvolve1d: could not allocate LUT\n");
    yi_LUT += halflen;
    for (i = -halflen; i < height + halflen; i++) {
      if (i < 0)
        yi_LUT[i] = 0;
      else if (i >= height)
        yi_LUT[i] = height - 1;
      else
        yi_LUT[i] = i;
    }
  }
  LUT_len = len;

  outPix = IMAGEFpix(J, 0, 0);
  if (axis == IMAGE_HORIZONTAL) {
    for (y = 0; y < height; y++) {
      inBase = IMAGEFpix(I, 0, y);
      for (x = 0; x < width; x++) {
        total = 0.0f;

        for (ki = k, i = 0; i < len; i++) total += *ki++ * *(inBase + xi_LUT[x + i - halflen]);

        *outPix++ = total;
      }
    }
  }
  else {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        inBase = IMAGEFpix(I, x, 0);
        total = 0.0f;

        for (ki = k, i = 0; i < len; i++) total += *ki++ * *(inBase + yi_LUT[y + i - halflen] * width);

        *outPix++ = total;
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void ImageCircularConvolve1d(IMAGE *I, IMAGE *J, float k[], int len, int axis)
{
  int x, y, width, height, halflen, xi, yi;
  int i;
  float *outPix;
  float *ki, total, *inBase;

  width = I->cols;
  height = I->rows;

  halflen = len / 2;

  outPix = IMAGEFpix(J, 0, 0);
  if (axis == IMAGE_HORIZONTAL) {
    for (y = 0; y < height; y++) {
      inBase = IMAGEFpix(I, 0, y);
      for (x = 0; x < width; x++) {
        total = 0.0f;

        for (ki = k, i = 0; i < len; i++) {
          xi = x + i - halflen;
          if (xi < 0) /* use circular topology */
            xi += width;
          if (xi >= width) xi -= width;
          total += *ki++ * *(inBase + xi);
        }

        *outPix++ = total;
      }
    }
  }
  else {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        inBase = IMAGEFpix(I, x, 0);
        total = 0.0f;

        for (ki = k, i = 0; i < len; i++) {
          yi = y + i - halflen;
          if (yi < 0) /* use circular topology */
            yi += height;
          if (yi >= height) yi -= height;
          total += *ki++ * *(inBase + yi * width);
        }

        *outPix++ = total;
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void ImageConvolve1dByte(IMAGE *I, IMAGE *J, float k[], int len, int axis)
{
  int x, y, width, height, halflen;
  int i;
  float *ki, total;
  ubyte *inBase, *outPix;
  static int *xi_LUT = NULL, LUT_width, LUT_height, *yi_LUT = NULL, LUT_len = 0;

  width = I->cols;
  height = I->rows;

  halflen = len / 2;

  if ((LUT_len != len) || (LUT_width != width)) {
    LUT_width = width;
    if (xi_LUT) free(xi_LUT - LUT_len / 2);
    xi_LUT = (int *)calloc(width + 2 * halflen, sizeof(int));
    if (!xi_LUT) ErrorExit(ERROR_NO_MEMORY, "ImageConvolve1d: could not allocate LUT\n");
    xi_LUT += halflen;
    for (i = -halflen; i < width + halflen; i++) {
      if (i < 0)
        xi_LUT[i] = 0;
      else if (i >= width)
        xi_LUT[i] = width - 1;
      else
        xi_LUT[i] = i;
    }
  }

  if ((LUT_len != len) || (LUT_height != height)) {
    LUT_height = height;
    if (yi_LUT) free(yi_LUT - LUT_len / 2);
    yi_LUT = (int *)calloc(height + 2 * halflen, sizeof(int));
    if (!yi_LUT) ErrorExit(ERROR_NO_MEMORY, "ImageConvolve1d: could not allocate LUT\n");
    yi_LUT += halflen;
    for (i = -halflen; i < height + halflen; i++) {
      if (i < 0)
        yi_LUT[i] = 0;
      else if (i >= height)
        yi_LUT[i] = height - 1;
      else
        yi_LUT[i] = i;
    }
  }
  LUT_len = len;

  outPix = IMAGEpix(J, 0, 0);
  if (axis == IMAGE_HORIZONTAL) {
    for (y = 0; y < height; y++) {
      inBase = IMAGEpix(I, 0, y);
      for (x = 0; x < width; x++) {
        total = 0.0f;

        for (ki = k, i = 0; i < len; i++) total += *ki++ * (float)*(inBase + xi_LUT[x + i - halflen]);

        *outPix++ = (ubyte)total;
      }
    }
  }
  else {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        inBase = IMAGEpix(I, x, 0);
        total = 0.0f;

        for (ki = k, i = 0; i < len; i++) total += *ki++ * (float)*(inBase + yi_LUT[y + i - halflen] * width);

        *outPix++ = (ubyte)total;
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              use the kernel described by Burt and Adelson to reduce
              the pyramid 1 level.
----------------------------------------------------------------------*/
#define KERNEL_SIZE 5
#define K_A 0.4f
static float kernel[KERNEL_SIZE] = {0.25f - K_A / 2.0f, .25f, K_A, 0.25f, 0.25f - K_A / 2.0f};

IMAGE *ImageReduce(IMAGE *Isrc, IMAGE *Idst)
{
  int rows, cols;
  static IMAGE *Itmp = NULL;

  rows = Isrc->rows;
  cols = Isrc->cols;
  if (!ImageCheckSize(Isrc, Itmp, rows, cols, 0) && Itmp) {
    ImageFree(&Itmp);
    Itmp = NULL;
  }

  if (!Itmp) {
    Itmp = ImageAlloc(rows, cols, PFFLOAT, 1);
    if (!Itmp) return (NULL);
  }
  else {
    ImageSetSize(Itmp, rows, cols);
    ImageClearArea(Itmp, 0, 0, -1, -1, 0.0f, -1);
  }

  rows /= 2;
  cols /= 2;

  if (!ImageCheckSize(Isrc, Idst, rows, cols, 0)) {
    if (Idst) ImageFree(&Idst);
    Idst = ImageAlloc(rows, cols, PFFLOAT, 1);
  }
  else
    ImageSetSize(Idst, rows, cols);

  /* blur vertically */
  ImageConvolve1d(Isrc, Itmp, kernel, KERNEL_SIZE, IMAGE_VERTICAL);
  ImageReduce1d(Itmp, Idst, kernel, KERNEL_SIZE, IMAGE_HORIZONTAL);
#if 0
  {
    char str[100] ;
    sprintf(str, "tmp%d.hipl", rows*2) ;
    ImageWrite(Itmp, str) ;
    sprintf(str, "out%d.hipl", rows*2) ;
    ImageWrite(Idst, str) ;
  }
#endif
  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             use k[] to scale the image down by 2.
----------------------------------------------------------------------*/
void ImageReduce1d(IMAGE *I, IMAGE *J, float k[], int len, int axis)
{
  int x, y, i, Jwidth, Jheight, xi, yi, halflen, Iwidth, Iheight;
  float total;

  Jwidth = J->cols;
  Jheight = J->rows;
  Iwidth = I->cols;
  Iheight = I->rows;

  halflen = (len - 1) / 2;
  if (axis == IMAGE_HORIZONTAL) {
    for (y = 0; y < Jheight; y++) {
      yi = 2 * y;
      for (x = 0; x < Jwidth; x++) {
        total = 0.0f;

        if (x >= 3 && x <= 6) i = 0;

        for (i = 0; i < len; i++) {
          /* Neumann boundary conditions */
          xi = 2 * x + i - halflen;
#if 0
          if (xi < 0)
            xi = -xi ;
          else if (xi >= Iwidth)
            xi = Iwidth - (xi - Iwidth+1) ;
#else
          if (xi < 0)
            xi = 0;
          else if (xi >= Iwidth)
            xi = Iwidth - 1;
#endif

          total = total + k[i] * *IMAGEFpix(I, xi, yi);
        }
        *IMAGEFpix(J, x, y) = total;
      }
    }
  }
  else {
    for (y = 0; y < Jheight; y++) {
      for (x = 0; x < Jwidth; x++) {
        total = 0.0f;
        xi = 2 * x;

#if 0
        if (((x == 9) && (y == 127)) ||
            ((y == 9) && (x == 127)))
          i = 0 ;
#endif

        for (i = 0; i < len; i++) {
          /* Neumann boundary conditions */
          yi = 2 * y + i - halflen;
#if 0
          if (yi < 0)
            yi = -yi ;
          else if (yi >= Iheight)
            yi = Iheight - (yi - Iheight+1) ;
#else
          if (yi < 0)
            yi = 0;
          else if (yi >= Iheight)
            yi = Iheight - 1;
#endif

          total = total + k[i] * *IMAGEFpix(I, xi, yi);
        }
        *IMAGEFpix(J, x, y) = total;
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             construct a gaussian bump which tails to 0. Returns
             an image which is ((4*xsigma)+1) x ((4*ysigma)+1)
----------------------------------------------------------------------*/
IMAGE *ImageGaussian(float xsigma, float ysigma)
{
  IMAGE *image;
  float norm, ytwo_sigma, xtwo_sigma, fx, fy, k;
  int x, y, xlen, ylen, xhalf, yhalf;

  /* build the kernel in k */
  xlen = (int)nint(6.0f * xsigma) + 1;
  ylen = (int)nint(6.0f * ysigma) + 1;
  if (EVEN(xlen)) xlen++;
  if (EVEN(ylen)) ylen++;
  xhalf = xlen / 2;
  yhalf = ylen / 2;
  image = ImageAlloc(xlen, ylen, PFFLOAT, 1);

  norm = 0.0f;
  xtwo_sigma = 2.0f * xsigma;
  ytwo_sigma = 2.0f * ysigma;

  for (x = 0; x < xlen; x++) {
    fx = (float)(x - xhalf);
    if (fabs(fx) <= xtwo_sigma)
      k = (float)exp((double)(-fx * fx / (xtwo_sigma * xsigma)));
    else if (xtwo_sigma < (float)fabs(fx) && (float)fabs(fx) <= 4.0f * xsigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0 - fabs(fx) / (double)xsigma, 4.0);
    else
      k = 0;

    for (y = 0; y < ylen; y++) *IMAGEFpix(image, x, y) = k;
  }

  for (y = 0; y < ylen; y++) {
    fy = (float)(y - yhalf);
    if (fabs(fy) <= ytwo_sigma)
      k = (float)exp((double)(-fy * fy / (ytwo_sigma * ysigma)));
    else if (ytwo_sigma < fabs(fy) && fabs(fy) <= 4 * ysigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0 - fabs(fy) / (double)ysigma, 4.0);
    else
      k = 0;

    for (x = 0; x < xlen; x++) {
      *IMAGEFpix(image, x, y) *= k;
      norm += *IMAGEFpix(image, x, y);
    }
  }

  /* normalize kernel to sum to 1 */
  for (x = 0; x < xlen; x++) {
    for (y = 0; y < ylen; y++) *IMAGEFpix(image, x, y) /= norm;
  }

  return (image);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             construct a splined gaussian bump which tails to 0.
             Returns an image which is (8*sigma)+1
             (Nitzberg and Shiota, 1993)
----------------------------------------------------------------------*/
IMAGE *ImageGaussian1d(float sigma, int max_len)
{
  IMAGE *image;
  float norm, two_sigma, fx, k;
  int x, half, len;

  /* build the kernel in k */
  len = (int)nint(8.0f * sigma) + 1;
  if (ISEVEN(len)) /* ensure it's even */
    len++;
  if (max_len && (max_len < len)) len = max_len;
  half = len / 2;
  image = ImageAlloc(1, len, PFFLOAT, 1);

  norm = 0.0f;
  two_sigma = 2.0f * sigma;

  for (norm = 0.0f, x = 0; x < len; x++) {
    fx = (float)(x - half);
    if (fabs(fx) <= two_sigma)
      k = (float)exp((double)(-fx * fx / (two_sigma * sigma)));
    else if (two_sigma < fabs(fx) && fabs(fx) <= 4.0f * sigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0f - fabs(fx) / (double)sigma, 4.0);
    else
      k = 0;

    *IMAGEFpix(image, x, 0) = k;
    norm += k;
  }

  /* normalize kernel to sum to 1 */
  for (x = 0; x < len; x++) *IMAGEFpix(image, x, 0) /= norm;

  return (image);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int ImageSobel(IMAGE *Isrc, IMAGE *gradImage, IMAGE *dxImage, IMAGE *dyImage)
{
  static IMAGE *xImage = NULL, *yImage = NULL;
  IMAGE *Iin;
  int x, y, rows, cols;
  float *xpix, *ypix, *gradpix = NULL, xval, yval, gval;

  if (Isrc->pixel_format != PFFLOAT) {
    Iin = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);
    ImageCopy(Isrc, Iin);
  }
  else
    Iin = Isrc;

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!dxImage) {
    if (!ImageCheckSize(Isrc, xImage, 0, 0, 0)) {
      if (xImage) ImageFree(&xImage);
      xImage = ImageAlloc(rows, cols, PFFLOAT, 1);
    }
    else
      ImageSetSize(xImage, rows, cols);

    dxImage = xImage;
  }

  if (!dyImage) {
    if (!ImageCheckSize(Isrc, yImage, 0, 0, 0)) {
      if (yImage) ImageFree(&yImage);
      yImage = ImageAlloc(rows, cols, PFFLOAT, 1);
    }
    else
      ImageSetSize(yImage, rows, cols);

    dyImage = yImage;
  }

  ImageSetSize(dxImage, rows, cols);
  ImageSetSize(dyImage, rows, cols);
  ImageSobelY(Isrc, dyImage);
  ImageSobelX(Isrc, dxImage);

  if (gradImage) {
    ImageSetSize(gradImage, rows, cols);
    gradpix = IMAGEFpix(gradImage, 0, 0);

    xpix = IMAGEFpix(dxImage, 0, 0);
    ypix = IMAGEFpix(dyImage, 0, 0);
    for (y = 0; y < rows; y++) {
      for (x = 0; x < cols; x++) {
        xval = *xpix++;
        yval = *ypix++;
        gval = (float)sqrt((double)(xval * xval + yval * yval));
        *gradpix++ = gval;
      }
    }
  }

  if (Iin != Isrc) ImageFree(&Iin);

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
/*
   -0.25    0    0.25
   -0.50    0    0.50
   -0.25    0    0.25
*/
#define FAST_SOBEL 1
#if !FAST_SOBEL
int ImageSobelX(IMAGE *Isrc, IMAGE *xImage)
{
  float *tl_pix, *ml_pix, *bl_pix, *tr_pix, *mr_pix, *br_pix, *outPtr;
  int rows, cols, row, col;

  rows = Isrc->rows;
  cols = Isrc->cols;
  outPtr = IMAGEFpix(xImage, 1, 1);
  tl_pix = IMAGEFpix(Isrc, 0, 0);
  ml_pix = IMAGEFpix(Isrc, 0, 1);
  bl_pix = IMAGEFpix(Isrc, 0, 2);
  tr_pix = IMAGEFpix(Isrc, 2, 0);
  mr_pix = IMAGEFpix(Isrc, 2, 1);
  br_pix = IMAGEFpix(Isrc, 2, 2);

  /* don't apply sobel to outer ring to pixels to avoid border effects */
  rows--;
  cols--;
  for (row = 1; row < rows; row++) {
    for (col = 1; col < cols; col++) {
      *outPtr++ = -.25f * *tl_pix++ - .5f * *ml_pix++ - .25f * *bl_pix++ + .25f * *tr_pix++ + .5f * *mr_pix++ +
                  .25f * *br_pix++;
    }
    outPtr += 2;
    tl_pix += 2;
    ml_pix += 2;
    bl_pix += 2;
    tr_pix += 2;
    mr_pix += 2;
    br_pix += 2;
  }

  return (0);
}
#else
/*----------------------------------------------------------------------
Parameters:

Description:
use overlapping windows to speed up sobel calculation

-0.25    0    0.25
-0.50    0    0.50
-0.25    0    0.25
----------------------------------------------------------------------*/
int ImageSobelX(IMAGE *Isrc, IMAGE *xImage)
{
  float *tr_pix, *mr_pix, *br_pix, *outPtr, left, middle, right;
  int rows, cols, row, col;

  rows = Isrc->rows;
  cols = Isrc->cols;
  outPtr = IMAGEFpix(xImage, 1, 1);

  tr_pix = IMAGEFpix(Isrc, 0, 0);
  mr_pix = IMAGEFpix(Isrc, 0, 1);
  br_pix = IMAGEFpix(Isrc, 0, 2);

  /* don't apply sobel to outer ring to pixels to avoid border effects */
  rows--;
  cols--;
  for (row = 1; row < rows; row++) {
    left = *tr_pix++ + 2.0f * *mr_pix++ + *br_pix++;
    middle = *tr_pix++ + 2.0f * *mr_pix++ + *br_pix++;

    for (col = 1; col < cols; col++) {
      right = *tr_pix++ + 2.0f * *mr_pix++ + *br_pix++;
      *outPtr++ = (right - left) * .125f;
      left = middle;
      middle = right;
    }
    outPtr += 2;
  }

  return (0);
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
   -0.25   -.50   -0.25
    0       0      0
    0.25    .50    0.25
----------------------------------------------------------------------*/
#if !FAST_SOBEL
int ImageSobelY(IMAGE *Isrc, IMAGE *yImage)
{
  float *tl_pix, *tm_pix, *tr_pix, *bl_pix, *bm_pix, *br_pix, *outPtr;
  int rows, cols, row, col;

  rows = Isrc->rows;
  cols = Isrc->cols;
  outPtr = IMAGEFpix(yImage, 1, 1);
  tl_pix = IMAGEFpix(Isrc, 0, 0);
  tm_pix = IMAGEFpix(Isrc, 1, 0);
  tr_pix = IMAGEFpix(Isrc, 2, 0);
  bl_pix = IMAGEFpix(Isrc, 0, 2);
  bm_pix = IMAGEFpix(Isrc, 1, 2);
  br_pix = IMAGEFpix(Isrc, 2, 2);

  /* don't apply sobel to outer ring to pixels to avoid border effects */
  rows--;
  cols--;
  for (row = 1; row < rows; row++) {
    for (col = 1; col < cols; col++) {
      *outPtr++ = -.25f * *tl_pix++ - .5f * *tm_pix++ - .25f * *tr_pix++ + .25f * *bl_pix++ + .5f * *bm_pix++ +
                  .25f * *br_pix++;
    }
    outPtr += 2;
    tl_pix += 2;
    tm_pix += 2;
    tr_pix += 2;
    bl_pix += 2;
    bm_pix += 2;
    br_pix += 2;
  }

  return (0);
}
#else
/*----------------------------------------------------------------------
Parameters:

Description:
use overlapping windows to speed up sobel calculation

-0.25   -.50   -0.25
0       0      0
0.25    .50    0.25
----------------------------------------------------------------------*/
int ImageSobelY(IMAGE *Isrc, IMAGE *yImage)
{
  float *tr_pix, *br_pix, *outPtr, left, middle, right;
  int rows, cols, row, col;

  rows = Isrc->rows;
  cols = Isrc->cols;
  outPtr = IMAGEFpix(yImage, 1, 1);

  tr_pix = IMAGEFpix(Isrc, 0, 0);
  br_pix = IMAGEFpix(Isrc, 0, 2);

  /* don't apply sobel to outer ring of pixels to avoid border effects */
  rows--;
  cols--;
  for (row = 1; row < rows; row++) {
    left = *br_pix++ - *tr_pix++;
    middle = *br_pix++ - *tr_pix++;

    for (col = 1; col < cols; col++) {
      right = *br_pix++ - *tr_pix++;
      *outPtr++ = (right + 2.0f * middle + left) * .125f;
      left = middle;
      middle = right;
    }
    outPtr += 2;
  }

  return (0);
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageXDerivative(IMAGE *Isrc, IMAGE *xImage)
{
  if (!xImage) xImage = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);

  ImageConvolve3x3(Isrc, sx, xImage);

  return (xImage);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageYDerivative(IMAGE *Isrc, IMAGE *yImage)
{
  if (!yImage) yImage = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1);

  ImageConvolve3x3(Isrc, sy, yImage);

  return (yImage);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int ImageConvolve3x3(IMAGE *Isrc, float kernel[], IMAGE *outImage)
{
  int rows, cols, x, y, xk, yk, k, xi, yi, frame;
  float *fkpix, sum, *fopix, *fipix;

  rows = Isrc->rows;
  cols = Isrc->cols;

  for (frame = 0; frame < Isrc->num_frame; frame++) {
    switch (Isrc->pixel_format) {
      case PFFLOAT:
        fopix = (float *)IMAGEFseq_pix(outImage, 0, 0, frame);
        for (y = 0; y < rows; y++) {
          for (x = 0; x < cols; x++, fopix++) {
            fkpix = kernel;
            for (sum = 0.0f, k = 0, yk = -1; yk <= 1; yk++) {
              yi = y + yk; /* image coordinate */
              if (yi < 0)
                yi = 0;
              else if (yi >= rows)
                yi = rows - 1;

              for (xk = -1; xk <= 1; xk++, k++, fkpix++) {
                xi = x + xk; /* image coordinate */
                if (xi < 0)
                  xi = 0;
                else if (xi >= cols)
                  xi = cols - 1;
                fipix = IMAGEFseq_pix(Isrc, xi, yi, frame);
                sum = sum + *fipix * *fkpix;
              }
            }
            *fopix = sum;
          }
        }
        break;
      default:
        fprintf(stderr, "ImageConvolve3x3: unsupported pixel format %d\n", Isrc->pixel_format);
        exit(-1);
        break;
    }
  }

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int xoffsets[] = {0, -1, 0, 1, 0};
static int yoffsets[] = {-1, 0, 0, 0, 1};
#define ONE_EIGHTH (1.0f / 8.0f)
static float weights[] = {ONE_EIGHTH, ONE_EIGHTH, -0.5f, ONE_EIGHTH, ONE_EIGHTH};
#define LAPLACIAN_POINTS (sizeof(xoffsets) / sizeof(xoffsets[0]))

IMAGE *ImageLaplacian(IMAGE *Isrc, IMAGE *outImage)
{
  int rows, cols, x, y, xi, yi;
  float *fkpix, sum, *fopix, fival;

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!outImage) outImage = ImageAlloc(rows, cols, PFFLOAT, 1);

  switch (Isrc->pixel_format) {
    case PFFLOAT:
      fopix = (float *)IMAGEFpix(outImage, 0, 0);
      for (y = 0; y < rows; y++) {
        for (x = 0; x < cols; x++, fopix++) {
          fkpix = weights;
          sum = 0.0f;
          for (unsigned int i = 0; i < LAPLACIAN_POINTS; i++) {
            yi = y + yoffsets[i]; /* image coordinate */
            if (yi < 0)
              yi = 0;
            else if (yi >= rows)
              yi = rows - 1;

            xi = x + xoffsets[i]; /* image coordinate */
            if (xi < 0)
              xi = 0;
            else if (xi >= cols)
              xi = cols - 1;
            fival = *IMAGEFpix(Isrc, xi, yi);
            sum += fival * *fkpix++;
          }
          *fopix = sum;
        }
      }
      break;
    default:
      fprintf(stderr, "ImageLaplacian: unsupported pixel format %d\n", Isrc->pixel_format);
      exit(-1);
      break;
  }

  return (outImage);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageSigmaFilter(IMAGE *Isrc, int wsize, float nsigma, IMAGE *Ioffset, IMAGE *Idst)
{
  float *sort_array;
  int x0, y0, rows, cols, x, y, whalf, yc, dx, dy, frame, wsq, w, npix;
  float *sptr, *outPix, min_val, max_val, *inPix, val, mean, sigma, sigma_thresh, filter_val;
  IMAGE *Iout, *Iin;
  int xc;

  if (nsigma <= 0.0f) nsigma = 2.0f;
  rows = Isrc->rows;
  cols = Isrc->cols;
  if (Isrc->pixel_format != PFFLOAT) {
    Iin = ImageAlloc(rows, cols, PFFLOAT, 1);
    ImageCopy(Isrc, Iin);
    ImageValRange(Isrc, &min_val, &max_val);
  }
  else
    Iin = Isrc;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else
    Iout = Idst;

  wsq = wsize * wsize;
  whalf = (wsize - 1) / 2;

  /* create a static array for sorting pixels in */
  sort_array = (float *)calloc(wsq, sizeof(float));

  for (frame = 0; frame < Iin->num_frame; frame++) {
    outPix = IMAGEFseq_pix(Iout, 0, 0, frame);
    for (y0 = 0; y0 < rows; y0++) {
      for (x0 = 0; x0 < cols; x0++) {
        /*
                 x and y are in window coordinates, while xc and yc are in image
                 coordinates.
         */
        if (Ioffset) {
          dx = nint(*IMAGEFpix(Ioffset, x0, y0));
          dy = nint(*IMAGEFseq_pix(Ioffset, x0, y0, 1));
        }
        else
          dx = dy = 0;

        for (sptr = sort_array, y = -whalf; y <= whalf; y++) {
          /* reflect across the boundary */
          yc = y + y0 + dy;
          if (yc < 0)
            yc = 0;
          else if (yc >= rows)
            yc = rows - 1;

          inPix = IMAGEFseq_pix(Iin, 0, yc, frame);
          for (x = -whalf; x <= whalf; x++) {
            xc = x0 + x + dx;
            if (xc < 0)
              xc = 0;
            else if (xc >= cols)
              xc = cols - 1;

            *sptr++ = *(inPix + xc);
          }
        }

        /* calculate mean and standard deviation in window */
        mean = sigma = 0.0f;
        for (w = 0, sptr = sort_array; w < wsq; w++) mean += *sptr++;
        mean /= wsq;
        for (w = 0, sptr = sort_array; w < wsq; w++) {
          val = *sptr++;
          val = (val - mean);
          val *= val;
          sigma += val;
        }
        sigma = sqrt(sigma) / wsq;

        /* calculate average of all pixels within nsigma std of mean */
        sigma_thresh = nsigma * sigma;
        filter_val = 0.0f;
        npix = 0;
        for (w = 0, sptr = sort_array; w < wsq; w++) {
          val = *sptr++;
          if (fabs(val - mean) <= sigma_thresh) /* include in filter */
          {
            filter_val += val;
            npix++;
          }
        }
        if (npix)
          *outPix++ = filter_val / npix;
        else
          *outPix++ = sort_array[(int)(wsq / 2.0f)];
      }
    }
  }

  free(sort_array);
  if (Idst != Iout) {
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  if (Isrc != Iin) ImageFree(&Iin);

  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int imageMeanFilter3x3(IMAGE *Isrc, IMAGE *Idst);
static int imageMeanFilter2x2(IMAGE *Isrc, IMAGE *Idst);
IMAGE *ImageMeanFilter(IMAGE *Isrc, int wsize, IMAGE *Idst)
{
  IMAGE *Iout;
  int rows, cols;

  if (Isrc->pixel_format != PFFLOAT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ImageMeanFilter: unsupported input format %d", Isrc->pixel_format));

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else
    Iout = Idst;

  switch (wsize) {
    case 3:
      imageMeanFilter3x3(Isrc, Idst);
      break;
    case 2:
      imageMeanFilter2x2(Isrc, Idst);
      break;
    default:
      ErrorReturn(Idst, (ERROR_UNSUPPORTED, "ImageMeanFilter: unsupported filter size %d", wsize));
      break; /* never used */
  }

  if (Iout != Idst) {
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int imageMeanFilter3x3(IMAGE *Isrc, IMAGE *Idst)
{
  float left, middle, right, *top, *center, *bottom, *out;
  int rows, cols, row, col, cols_minus_1, rows_minus_1;

  rows = Isrc->rows;
  cols = Isrc->cols;

  /* do 4 corner points separately and 1st and last rows and cols separately */

  /* top right corner */
  left = *IMAGEFpix(Isrc, cols - 2, 0) + *IMAGEFpix(Isrc, cols - 2, 1);
  middle = *IMAGEFpix(Isrc, cols - 1, 0) + *IMAGEFpix(Isrc, cols - 1, 1);
  *IMAGEFpix(Idst, cols - 1, 0) = (left + middle) / 4.0f;

  /* bottom right corner */
  left = *IMAGEFpix(Isrc, cols - 2, rows - 2) + *IMAGEFpix(Isrc, cols - 2, rows - 1);
  middle = *IMAGEFpix(Isrc, cols - 1, rows - 2) + *IMAGEFpix(Isrc, cols - 1, rows - 1);
  *IMAGEFpix(Idst, cols - 1, rows - 1) = (left + middle) / 4.0f;

  /* top left corner, also initializes for first row */
  left = *IMAGEFpix(Isrc, 0, 0) + *IMAGEFpix(Isrc, 0, 1);
  middle = *IMAGEFpix(Isrc, 1, 0) + *IMAGEFpix(Isrc, 1, 1);
  *IMAGEFpix(Idst, 0, 0) = (left + middle) / 4.0f;

  /* skip 1st and last elements */
  cols_minus_1 = cols - 1;
  rows_minus_1 = rows - 1;

  /* first do top row */
  out = IMAGEFpix(Idst, 1, 0);
  center = IMAGEFpix(Isrc, 2, 0);
  bottom = IMAGEFpix(Isrc, 2, 1);

  for (col = 1; col < cols_minus_1; col++) {
    right = *center++ + *bottom++;
    *out++ = (left + middle + right) / 6.0f;
    left = middle;
    middle = right;
  }

  /* bottom left corner and initialize for bottom row */
  left = *IMAGEFpix(Isrc, 0, rows - 2) + *IMAGEFpix(Isrc, 0, rows - 1);
  middle = *IMAGEFpix(Isrc, 1, rows - 2) + *IMAGEFpix(Isrc, 1, rows - 1);
  *IMAGEFpix(Idst, 0, rows - 1) = (left + middle) / 4.0f;

  /* do bottom row */
  out = IMAGEFpix(Idst, 1, rows - 1);
  center = IMAGEFpix(Isrc, 2, rows - 1);
  top = IMAGEFpix(Isrc, 2, rows - 2);

  for (col = 1; col < cols_minus_1; col++) {
    right = *center++ + *top++;
    *out++ = (left + middle + right) / 6.0f;
    left = middle;
    middle = right;
  }

  /* do 1st column */
  left = *IMAGEFpix(Isrc, 0, 0) + *IMAGEFpix(Isrc, 1, 0);
  middle = *IMAGEFpix(Isrc, 0, 1) + *IMAGEFpix(Isrc, 1, 1);
  out = IMAGEFpix(Idst, 0, 1);
  center = IMAGEFpix(Isrc, 0, 2);
  bottom = IMAGEFpix(Isrc, 1, 2);

  for (row = 1; row < rows_minus_1; row++) {
    right = *center + *bottom;
    *out = (left + middle + right) / 6.0f;
    out += cols;
    center += cols;
    bottom += cols;
    left = middle;
    middle = right;
  }

  /*
     the naming conventions are messed up for the first and last
     column. Everything is basically transposed with left referring to top,
     etc... Sorry.
   */
  /* do last column */
  left = *IMAGEFpix(Isrc, cols - 2, 0) + *IMAGEFpix(Isrc, cols - 1, 0);
  middle = *IMAGEFpix(Isrc, cols - 2, 1) + *IMAGEFpix(Isrc, cols - 1, 1);

  center = IMAGEFpix(Isrc, cols - 1, 2);
  bottom = IMAGEFpix(Isrc, cols - 2, 2);

  out = IMAGEFpix(Idst, cols - 1, 1);

  for (row = 1; row < rows_minus_1; row++) {
    right = *center + *bottom;
    *out = (left + middle + right) / 6.0f;
    out += cols;
    center += cols;
    bottom += cols;
    left = middle;
    middle = right;
  }

  /* now do the image */
  out = IMAGEFpix(Idst, 1, 1);
  top = IMAGEFpix(Isrc, 2, 0);
  center = IMAGEFpix(Isrc, 2, 1);
  bottom = IMAGEFpix(Isrc, 2, 2);

  for (row = 1; row < rows_minus_1; row++) {
    left = *IMAGEFpix(Isrc, 0, row - 1) + *IMAGEFpix(Isrc, 0, row) + *IMAGEFpix(Isrc, 0, row + 1);
    middle = *IMAGEFpix(Isrc, 1, row - 1) + *IMAGEFpix(Isrc, 1, row) + *IMAGEFpix(Isrc, 1, row + 1);

    for (col = 1; col < cols_minus_1; col++) {
      right = *top++ + *center++ + *bottom++;
      *out++ = (left + middle + right) / 9.0f;
      left = middle;
      middle = right;
    }
    top += 2;
    bottom += 2;
    center += 2;
    out += 2;
  }

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              convolve with a 2x2 averaging filter using overlap for
              efficient implementation
----------------------------------------------------------------------*/
static int imageMeanFilter2x2(IMAGE *Isrc, IMAGE *Idst)
{
  float left, right, *top, *bottom, *out;
  int rows, cols, row, col, cols_minus_1, rows_minus_1;

  rows = Isrc->rows;
  cols = Isrc->cols;

  /*
     do last column and bottom right point separately using backwards
     averaging instead of forward
     */

  /* bottom right corner, set it to average of 4 neighbors */
  left = *IMAGEFpix(Isrc, cols - 2, rows - 2) + *IMAGEFpix(Isrc, cols - 2, rows - 1);
  right = *IMAGEFpix(Isrc, cols - 1, rows - 2) + *IMAGEFpix(Isrc, cols - 1, rows - 1);
  *IMAGEFpix(Idst, cols - 1, rows - 1) = (left + right) * 0.25f;

  /* do last column */

  /* skip last column and last row */
  cols_minus_1 = cols - 1;
  rows_minus_1 = rows - 1;

  /* do last column */
  left = *IMAGEFpix(Isrc, cols - 2, 0) + *IMAGEFpix(Isrc, cols - 1, 0);

  bottom = IMAGEFpix(Isrc, cols - 1, 1); /* actually right */
  top = IMAGEFpix(Isrc, cols - 2, 1);    /* actually left */
  out = IMAGEFpix(Idst, cols - 1, 0);

  for (row = 0; row < rows_minus_1; row++) {
    right = *top + *bottom;
    *out = (left + right) * 0.25f;
    out += cols;
    top += cols;
    bottom += cols;
    left = right;
  }

  /* do last row */
  out = IMAGEFpix(Idst, 0, rows - 1);
  top = IMAGEFpix(Isrc, 0, rows - 2);
  bottom = IMAGEFpix(Isrc, 0, rows - 1);
  left = *top++ + *bottom++;

  for (col = 0; col < cols_minus_1; col++) {
    right = *top++ + *bottom++;
    *out++ = (left + right) * 0.25f;
    left = right;
  }

  /* now do the image */
  out = IMAGEFpix(Idst, 0, 0);
  top = IMAGEFpix(Isrc, 0, 0);
  bottom = IMAGEFpix(Isrc, 0, 1);

  for (row = 0; row < rows_minus_1; row++) {
    left = *top++ + *bottom++;

    for (col = 0; col < cols_minus_1; col++) {
      right = *top++ + *bottom++;
      *out++ = (left + right) * 0.25f;
      left = right;
    }
    out++;
  }

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             gray-scale 3x3 dilation - local max
----------------------------------------------------------------------*/
IMAGE *ImageGreyDilate(IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE *Iout;
  int rows, cols, x, y, xk, yk, xi, yi, frame;
  float image_fmin, image_fmax, fmax, *fopix, *fipix;

  if (Isrc->pixel_format != PFFLOAT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ImageDilate: unsupported input format %d", Isrc->pixel_format));

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else
    Iout = Idst;

  rows = Isrc->rows;
  cols = Isrc->cols;

  ImageValRange(Isrc, &image_fmin, &image_fmax);
  for (frame = 0; frame < Isrc->num_frame; frame++) {
    switch (Isrc->pixel_format) {
      case PFFLOAT:
        fopix = (float *)IMAGEFseq_pix(Iout, 0, 0, frame);
        for (y = 0; y < rows; y++) {
          for (x = 0; x < cols; x++, fopix++) {
            fmax = image_fmin;
            for (yk = -1; yk <= 1; yk++) {
              yi = y + yk; /* image coordinate */
              if (yi < 0)
                yi = 0;
              else if (yi >= rows)
                yi = rows - 1;

              for (xk = -1; xk <= 1; xk++) {
                xi = x + xk; /* image coordinate */
                if (xi < 0)
                  xi = 0;
                else if (xi >= cols)
                  xi = cols - 1;
                fipix = IMAGEFseq_pix(Isrc, xi, yi, frame);
                if (*fipix > fmax) fmax = *fipix;
              }
            }
            *fopix = fmax;
          }
        }
        break;
      default:
        ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ImageDilate: unsupported pixel format %d\n", Isrc->pixel_format));
        exit(-1);
        break;
    }
  }

  if (Iout != Idst) {
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageGreyErode(IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE *Iout;
  int rows, cols, x, y, xk, yk, xi, yi, frame;
  float image_fmin, image_fmax, fmin, *fopix, *fipix;

  if (Isrc->pixel_format != PFFLOAT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ImageDilate: unsupported input format %d", Isrc->pixel_format));

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else
    Iout = Idst;

  rows = Isrc->rows;
  cols = Isrc->cols;

  ImageValRange(Isrc, &image_fmin, &image_fmax);
  for (frame = 0; frame < Isrc->num_frame; frame++) {
    switch (Isrc->pixel_format) {
      case PFFLOAT:
        fopix = (float *)IMAGEFseq_pix(Iout, 0, 0, frame);
        for (y = 0; y < rows; y++) {
          for (x = 0; x < cols; x++, fopix++) {
            fmin = image_fmax;
            for (yk = -1; yk <= 1; yk++) {
              yi = y + yk; /* image coordinate */
              if (yi < 0)
                yi = 0;
              else if (yi >= rows)
                yi = rows - 1;

              for (xk = -1; xk <= 1; xk++) {
                xi = x + xk; /* image coordinate */
                if (xi < 0)
                  xi = 0;
                else if (xi >= cols)
                  xi = cols - 1;
                fipix = IMAGEFseq_pix(Isrc, xi, yi, frame);
                if (*fipix < fmin) fmin = *fipix;
              }
            }
            *fopix = fmin;
          }
        }
        break;
      default:
        ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ImageErode: unsupported pixel format %d\n", Isrc->pixel_format));
        break;
    }
  }

  if (Iout != Idst) {
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  return (Idst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *ImageCorrelate(IMAGE *Itemplate, IMAGE *Isrc, int zeropad, IMAGE *Icorr)
{
  IMAGE *Iconj, *Ifcorr, *Ifsrc, *Ireal, *Iimag, *Iftmp;
  int ecode;

#if 0
  if (zeropad)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "ImageCorrelate: zero padding unsupported")) ;
#else
  if (zeropad) Isrc = ImageZeroPad(Isrc, NULL);
#endif

  /* check to see if the template as already been FTed */
  if (Itemplate->pixel_format != PFCOMPLEX && Itemplate->pixel_format != PFDBLCOM)
    Iftmp = ImageDFT(Itemplate, NULL);
  else
    Iftmp = Itemplate;

  Iconj = ImageConjugate(Iftmp, NULL);
  Ifsrc = ImageDFT(Isrc, NULL);
  Ifcorr = ImageMul(Iconj, Ifsrc, NULL);
  Icorr = ImageInverseDFT(Ifcorr, Icorr);

  if (Icorr->pixel_format == PFCOMPLEX) {
    /* flipquad can't handle complex images, do it separately */
    Ireal = ImageAlloc(Icorr->rows, Icorr->cols, PFFLOAT, 1);
    Iimag = ImageAlloc(Icorr->rows, Icorr->cols, PFFLOAT, 1);
    ImageSplit(Icorr, Ireal, Iimag);
    ecode = h_flipquad(Ireal, Ireal);
    if (ecode != HIPS_OK) ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode);
    ecode = h_flipquad(Iimag, Iimag);
    if (ecode != HIPS_OK) ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode);
    ImageCombine(Ireal, Iimag, Icorr);
    ImageFree(&Ireal);
    ImageFree(&Iimag);
  }
  else {
    ecode = h_flipquad(Icorr, Icorr);
    if (ecode != HIPS_OK) ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode);
  }

  if (Iftmp != Itemplate) /* allocated an image for fft */
    ImageFree(&Iftmp);
  ImageFree(&Iconj);
  ImageFree(&Ifcorr);
  ImageFree(&Ifsrc);

  return (Icorr);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a spatial correlation of Ikernel with Isrc around the
           region of the point x0,y0. The size of the correlation region is
           given by wsize
------------------------------------------------------*/
IMAGE *ImageCorrelateRegion(IMAGE *Isrc, IMAGE *Ikernel, IMAGE *Idst, int row0, int col0, int wsize)
{
  int col_offset, row_offset, row, col, rows, cols, whalf, krow0, kcol0, krow, kcol, drow, dcol;
  // int rowstart, rowend, colstart, colend;
  CPIX *src, *kernel;
  float sreal, simag, kreal, kimag, val, *dst, total;

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  kcol0 = col0 - Isrc->cols / 2;
  krow0 = row0 - Isrc->rows / 2;

  whalf = (wsize - 1) / 2;
  // rowstart = MAX(row0 - whalf, 0);
  // colstart = MAX(col0 - whalf, 0);
  // rowend = MIN(rows - 1, row0 + wsize);
  // colend = MIN(cols - 1, col0 + wsize);

  for (row_offset = -whalf; row_offset <= whalf; row_offset++) {
    drow = row0 + row_offset;
    if (drow < 0 || drow >= rows) continue;

    dst = IMAGEFpix(Idst, col0 - whalf, drow);
    for (col_offset = -whalf; col_offset <= whalf; col_offset++, dst++) {
      dcol = col0 + col_offset;
      if (dcol < 0 || dcol >= cols) continue;
      total = 0.0f;
      for (row = 0; row < rows; row++) {
        krow = row + row_offset + krow0;
        if (krow < 0 || krow >= rows) continue;

        src = IMAGECpix(Isrc, 0, row);
        for (col = 0; col < cols; col++, src++) {
          kcol = col + col_offset + kcol0;
          if (kcol < 0 || kcol >= cols) continue;
          kernel = IMAGECpix(Ikernel, krow, kcol);
          kreal = kernel->real;
          kimag = kernel->imag;
          sreal = src->real;
          simag = src->imag;
          val = kreal * sreal + kimag * simag; /* real part */
          total += val;
        }
      }
      *dst = total;
    }
  }

  return (Idst);
}
IMAGE *ImageLOGFilter(IMAGE *Isrc, float sigma, IMAGE *Idst)
{
  IMAGE *Ig, *Itmp;

  Ig = ImageGaussian1d(sigma, 0);

  Itmp = ImageConvolveGaussian(Isrc, Ig, NULL, 0);
  Idst = ImageLaplacian(Itmp, Idst);
#if 0
  ImageWrite(Itmp, "g.hipl") ;
  ImageWrite(Idst, "l.hipl") ;
#endif
  ImageFree(&Itmp);
  ImageFree(&Ig);
  return (Idst);
}
IMAGE *ImageDOGFilter(IMAGE *Isrc, float psigma, float nsigma, IMAGE *Idst)
{
  IMAGE *Igp, *Ign, *Itmpp, *Itmpn;

  Igp = ImageGaussian1d(psigma, 0);
  Ign = ImageGaussian1d(nsigma, 0);

  Itmpp = ImageConvolveGaussian(Isrc, Igp, NULL, 0);
  Itmpn = ImageConvolveGaussian(Isrc, Ign, NULL, 0);
  Idst = ImageSubtract(Itmpp, Itmpn, Idst);
#if 1
  ImageWrite(Itmpp, "gp.hipl");
  ImageWrite(Itmpn, "gn.hipl");
  ImageWrite(Idst, "dog.hipl");
#endif
  ImageFree(&Itmpp);
  ImageFree(&Itmpn);
  ImageFree(&Igp);
  ImageFree(&Ign);
  return (Idst);
}
