/*
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

/*
 *       FILE NAME:   offset.c
 *
 *       DESCRIPTION: routines for generating displacement vector fields.
 *                    Split from image.c
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        7/1/96
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
#include <unistd.h> /* for unlink */

#include "hips.h"

#include "diag.h"
#include "error.h"
#include "image.h"
#include "machine.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "utils.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#define MAX_STEPS 30
IMAGE *ImageCalculateOffset(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Ioffset)
{
  static IMAGE *Iorient = NULL, *Idir = NULL;
  Timer then;
  int msec;

  if (Gdiag & DIAG_TIMER) then.reset();
  Iorient = ImageOffsetOrientation(Ix, Iy, wsize, Iorient);

  if ((Gdiag & DIAG_WRITE) && (Gdiag & DIAG_VERBOSE)) ImageWrite(Iorient, "orient.hipl");

  if (Gdiag & DIAG_TIMER) {
    msec = then.milliseconds();
    fprintf(stderr, "orientation took  %2.3f sec\n", (float)msec / 1000.0f);
    then.reset();
  }

  Idir = ImageOffsetDirection(Ix, Iy, wsize, Iorient, Idir);

  if ((Gdiag & DIAG_WRITE) && (Gdiag & DIAG_VERBOSE)) ImageWrite(Idir, "dir.hipl");

  if (Gdiag & DIAG_TIMER) {
    msec = then.milliseconds();
    fprintf(stderr, "direction took  %2.3f sec\n", (float)msec / 1000.0f);
  }

  Ioffset = ImageOffsetMagnitude(Idir, Ioffset, MAX_STEPS);
  if ((Gdiag & DIAG_WRITE) && (Gdiag & DIAG_VERBOSE)) ImageWrite(Ioffset, "offset.hipl");

  return (Ioffset);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#define NS_FSCALE 100000.0f /* to scale to approx. same range as new stuff */

IMAGE *ImageCalculateNitShiOffset(IMAGE *Ix, IMAGE *Iy, int wsize, float mu, float c, IMAGE *Ioffset)
{
  int x0, y0, rows, cols, x, y, whalf, x_plus_off, y_plus_off;
  float vx, vy, vsq, c1, *g, *xpix, *ypix;
  static float *gaussian = NULL;
  static int w = 0;
  float gauss, fxpix, fypix, dot_product;
  int xc, yc;

  if ((Gdiag & DIAG_SHOW) && (Gdiag & DIAG_VERBOSE))
    fprintf(stderr, "ImageCalculateNitshiOffset: mu = %2.4f, c = %2.4f, wsize = %d\n", mu, c, wsize);
  rows = Ix->rows;
  cols = Ix->cols;

  if (!Ioffset) Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);

  mu *= mu;
  vsq = 0.0f; /* prevent compiler warning */

  whalf = (wsize - 1) / 2;
  c1 = c;

  /* create a local gaussian window */
  if (wsize != w) {
    free(gaussian);
    gaussian = NULL;
    w = wsize;
  }

  if (!gaussian) /* allocate a gaussian bump */
  {
    float den, norm;

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

  xpix = IMAGEFpix(Ioffset, 0, 0);
  ypix = IMAGEFseq_pix(Ioffset, 0, 0, 1);
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, xpix++, ypix++) {
      /*
        x and y are in window coordinates, while xc and yc are in image
        coordinates.
      */
      /* first build offset direction */
      vx = vy = 0.0f;
      for (g = gaussian, y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0;
        if ((yc < 0) || (yc >= rows)) {
          g += wsize;
          continue;
        }

        for (x = -whalf; x <= whalf; x++, g++) {
          xc = x0 + x;
          if ((xc < 0) || (xc >= cols)) continue;

          fxpix = *IMAGEFpix(Ix, xc, yc);
          fypix = *IMAGEFpix(Iy, xc, yc);
          dot_product = x * fxpix + y * fypix;
          gauss = *g;
          dot_product *= gauss;
          vx += (dot_product * fxpix);
          vy += (dot_product * fypix);
        }
      }

      /* calculated phi(V), only needed for original NitShi algorithm */
      vsq = vx * vx + vy * vy;

      vx = vx * c1 / (float)sqrt((double)(mu * mu + vsq));
      vy = vy * c1 / (float)sqrt((double)(mu * mu + vsq));
      x_plus_off = x0 - vx;
      y_plus_off = y0 - vy;
      if (x_plus_off < 0)
        vx = x0;
      else if (x_plus_off >= cols)
        vx = x0 - cols + 1;
      if (y_plus_off < 0)
        vy = y0;
      else if (y_plus_off >= rows)
        vy = y0 - rows + 1;
      *xpix = -vx;
      *ypix = -vy;
    }
  }

  if ((Gdiag & DIAG_WRITE) && (Gdiag & DIAG_VERBOSE)) ImageWrite(Ioffset, "nitshi_offset.hipl");
  return (Ioffset);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             normalize offset distances by searching for isopotential
             surfaces in offset direction, then modifying offset distance
             to be the distance to the isopotential surface.
----------------------------------------------------------------------*/
IMAGE *ImageNormalizeOffsetDistances(IMAGE *Isrc, IMAGE *Idst, int maxsteps)
{
  float *src_xpix, *src_ypix, *dst_xpix, *dst_ypix, slope, xf, yf, dot;
  int x0, y0, rows, cols, delta, i;
  int x, y;
  float dx, dy, odx, ody;

  if (!Idst) Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, Isrc->num_frame);

  rows = Isrc->rows;
  cols = Isrc->cols;

  src_xpix = IMAGEFpix(Isrc, 0, 0);
  src_ypix = IMAGEFseq_pix(Isrc, 0, 0, 1);
  dst_xpix = IMAGEFpix(Idst, 0, 0);
  dst_ypix = IMAGEFseq_pix(Idst, 0, 0, 1);

  /* for each point in the image */
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, src_xpix++, src_ypix++, dst_xpix++, dst_ypix++) {
      /*
            search for the first point in the offset direction who's dot product
            with the current offset vector is below some threshold.
      */
      dx = *src_xpix;
      dy = *src_ypix;
      if (FZERO(dx) && (FZERO(dy))) continue;

      if (fabs(dx) > fabs(dy)) /* use unit steps in x direction */
      {
        delta = nint(dx / (float)fabs(dx));
        slope = delta * dy / dx;
        for (i = 0, yf = (float)y0 + slope, x = x0 + delta; i <= maxsteps; x += delta, yf += slope, i++) {
          y = nint(yf);
          if (y <= 0 || y >= (rows - 1) || x <= 0 || x >= (cols - 1)) break;

          odx = *IMAGEFpix(Isrc, x, y);
          ody = *IMAGEFseq_pix(Isrc, x, y, 1);
          dot = odx * dx + ody * dy;
          if (dot <= 0) break;
        }
        x -= delta;
        y = nint(yf - slope);
      }
      else /* use unit steps in y direction */
      {
        delta = nint(dy / (float)fabs(dy));
        slope = delta * dx / dy;
        for (i = 0, xf = (float)x0 + slope, y = y0 + delta; i < maxsteps; y += delta, xf += slope, i++) {
          x = nint(xf);
          if (y <= 0 || y >= (rows - 1) || x <= 0 || x >= (cols - 1)) break;

          odx = *IMAGEFpix(Isrc, x, y);
          ody = *IMAGEFseq_pix(Isrc, x, y, 1);
          dot = odx * dx + ody * dy;
          if (dot <= 0) break;
        }
        y -= delta;
        x = nint(xf - slope);
      }

      *dst_xpix = x - x0;
      *dst_ypix = y - y0;
    }
  }

  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              smooth offsets by using winner-take-all algorithm in
              a neighborhood in the direction orthogonal to the
              offset.
----------------------------------------------------------------------*/


IMAGE *ImageSmoothOffsets(IMAGE *Isrc, IMAGE *Idst, int wsize)
{
  float *src_xpix, *src_ypix, *dst_xpix, *dst_ypix, slope, dx, dy, f, xf, yf, *wdx, *wdy, *wphase, *wmag, dist, mag;
  int x0, y0, rows, cols, x = 0, y, delta, i, whalf, *wx, *wy;

  wsize = 5;
  if (!Idst) Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, Isrc->num_frame);

  ImageCopy(Isrc, Idst);
  whalf = (wsize - 1) / 2;
  /*
  allocate five windows of the size specified by the user. Two to hold
  the offset vectors, one for the magnitude, one for the phase, and 1 for
  the voting weight.
  */
  wdx = (float *)calloc(wsize, sizeof(float));
  wdy = (float *)calloc(wsize, sizeof(float));
  wphase = (float *)calloc(wsize, sizeof(float));
  wmag = (float *)calloc(wsize, sizeof(float));
  wx = (int *)calloc(wsize, sizeof(float));
  wy = (int *)calloc(wsize, sizeof(float));

  rows = Isrc->rows;
  cols = Isrc->cols;

  src_xpix = IMAGEFpix(Isrc, 0, 0);
  src_ypix = IMAGEFseq_pix(Isrc, 0, 0, 1);
  dst_xpix = IMAGEFpix(Idst, 0, 0);
  dst_ypix = IMAGEFseq_pix(Idst, 0, 0, 1);

  /* for each point in the image */
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, src_xpix++, src_ypix++, dst_ypix++, dst_xpix++) {
      /* fill the offset vector array */
      dx = *src_xpix;
      dy = *src_ypix;

      /* calculate orthogonal slope = -dx/dy */
      f = dx;
      dx = dy;
      dy = -f;
      mag = (float)hypot(dx, dy);

      if (ISTINY(mag)) /* don't know what direction to search in */
        continue;
      else if (fabs(dx) > fabs(dy)) /* use unit steps in x direction */
      {
        delta = nint(dx / (float)fabs(dx));
        slope = delta * dy / dx; /* orthogonal slope */

        yf = (float)y0 - (float)whalf * slope;
        x = x0 - whalf * delta;
        for (i = 0; i < wsize; x += delta, yf += slope, i++) {
          y = nint(yf);
          wx[i] = x;
          wy[i] = y;
          if (y <= 0 || y >= (rows - 1) || x <= 0 || x >= (cols - 1))
            wdx[i] = wdy[i] = wmag[i] = wphase[i] = 0.0f;
          else {
            dx = *IMAGEFpix(Isrc, x, y);
            dy = *IMAGEFseq_pix(Isrc, x, y, 1);
            wdx[i] = dx;
            wdy[i] = dy;
            wmag[i] = (float)hypot((double)dx, (double)dy);
            wphase[i] = (float)latan2((double)dy, (double)dx);
          }
        }
      }
      else /* use unit steps in y direction */
      {
        delta = nint(dy / (float)fabs(dy));
        slope = delta * dx / dy; /* orthogonal slope */
        xf = (float)x0 - (float)whalf * slope;
        y = y0 - whalf * delta;
        for (i = 0; i < wsize; y += delta, xf += slope, i++) {
          wx[i] = x;
          wy[i] = y;
          x = nint(xf);
          if (y <= 0 || y >= (rows - 1) || x <= 0 || x >= (cols - 1))
            wdx[i] = wdy[i] = wmag[i] = wphase[i] = 0.0f;
          else {
            dx = *IMAGEFpix(Isrc, x, y);
            dy = *IMAGEFseq_pix(Isrc, x, y, 1);
            wdx[i] = dx;
            wdy[i] = dy;
            wmag[i] = (float)hypot((double)dx, (double)dy);
            wphase[i] = (float)latan2((double)dy, (double)dx);
          }
        }
      }

/*
check both left and right neighbors in direction orthogonal to our
offset direction. If our neighbors  left and right neighbors are coherent
in terms of direction (within some threshold of each other), and the
neighbor's diirection is not bracketed by them, then modify it.
*/
#define MIN_ANGLE_DIST RADIANS(30.0f)

      if (!ISTINY(wmag[0])) /* check 'left' neighbor */
      {
        dist = angleDistance(wphase[0], wphase[2]);
        if (dist < MIN_ANGLE_DIST) {

          dx = (wdx[2] + wdx[0]) / 2.0f;
          dy = (wdy[2] + wdy[0]) / 2.0f;
          /*
          check to see that the computed angle is not too close to being
          exactly 180 degrees out of alignment with the current one which
          would indicate we are close to the border of a coherent edge.
          */
          dist = angleDistance(wphase[2], wphase[1]);
          if ((dist - PI) > MIN_ANGLE_DIST) {
            *IMAGEFpix(Idst, wx[1], wy[1]) = dx;
            *IMAGEFseq_pix(Idst, wx[1], wy[1], 1) = dy;
          }
        }
      }

      if (!ISTINY(wmag[4])) /* check 'right' neighbor */
      {
        dist = angleDistance(wphase[4], wphase[2]);
        if (dist < MIN_ANGLE_DIST) {

          dx = (wdx[2] + wdx[4]) / 2.0f;
          dy = (wdy[2] + wdy[4]) / 2.0f;
          /*
          check to see that the computed angle is not too close to being
          exactly 180 degrees out of alignment with the current one which
          would indicate we are close to the border of a coherent edge.
          */
          dist = angleDistance(wphase[2], wphase[3]);
          if ((dist - PI) > MIN_ANGLE_DIST) {
            *IMAGEFpix(Idst, wx[3], wy[3]) = dx;
            *IMAGEFseq_pix(Idst, wx[3], wy[3], 1) = dy;
          }
        }
      }
    }
  }

  free(wx);
  free(wy);
  free(wdx);
  free(wdy);
  free(wphase);
  free(wmag);
  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              apply an offset vector to a filtered image.
----------------------------------------------------------------------*/
IMAGE *ImageApplyOffset(IMAGE *Isrc, IMAGE *Ioffset, IMAGE *Idst)
{
  int x, y, rows, cols, dx, dy;
  float *dst, *src, *dx_pix, *dy_pix;
  IMAGE *Iout, *Iin;

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else
    Iout = Idst;
  if (Isrc->pixel_format != PFFLOAT) {
    Iin = ImageAlloc(rows, cols, PFFLOAT, 1);
    ImageCopy(Isrc, Iin);
  }
  else
    Iin = Isrc;

  if (!ImageCheckSize(Isrc, Idst, 0, 0, 0)) ErrorReturn(NULL, (ERROR_SIZE, "ImageApplyOffset: dst not big enough"));

  dst = IMAGEFpix(Iout, 0, 0);
  dx_pix = IMAGEFpix(Ioffset, 0, 0);
  dy_pix = IMAGEFseq_pix(Ioffset, 0, 0, 1);

  for (y = 0; y < rows; y++) {
    for (x = 0; x < cols; x++) {
      dx = (int)*dx_pix++;
      dy = (int)*dy_pix++;
      src = IMAGEFpix(Iin, x + dx, y + dy);
      *dst++ = *src;
    }
  }

  if (Iin != Isrc) ImageFree(&Iin);
  if (Iout != Idst) {
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  return (Idst);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
              use an offset vector field to specify edge locations.
----------------------------------------------------------------------*/
IMAGE *ImageOffsetMedialAxis(IMAGE *Ioffset, IMAGE *Iedge)
{
  int x, y, rows, cols, dx, dy;
  float *dx_pix, *dy_pix;
  UCHAR *edge;
  IMAGE *Iout;

  rows = Ioffset->rows;
  cols = Ioffset->cols;

  if (!Iedge) Iedge = ImageAlloc(rows, cols, PFBYTE, 1);

  if (Iedge->pixel_format != PFBYTE)
    Iout = ImageAlloc(rows, cols, PFBYTE, 1);
  else
    Iout = Iedge;


  /* assume everything is an edge */
  ImageClearArea(Iout, -1, -1, -1, -1, 0.0f, -1);
  dx_pix = IMAGEFpix(Ioffset, 0, 0);
  dy_pix = IMAGEFseq_pix(Ioffset, 0, 0, 1);

  for (y = 0; y < rows; y++) {
    for (x = 0; x < cols; x++) {
      dx = (int)*dx_pix++;
      dy = (int)*dy_pix++;
      edge = IMAGEpix(Iout, x + dx, y + dy);
      (*edge)++; /* count # of times used */
    }
  }

  if (Iout != Iedge) {
    ImageCopy(Iout, Iedge);
    ImageFree(&Iout);
  }
  return (Iedge);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageCalculateOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Ioffset)
{
  int x0, y0, rows, cols, x, y, whalf;
  float vx, vy, *g, *xpix, *ypix;
  static float *gaussian = NULL;
  static int w = 0;
  float gauss, fxpix, fypix, dot_product;
  int xc, yc;

  rows = Ix->rows;
  cols = Ix->cols;

  if (!Ioffset) Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);

  whalf = (wsize - 1) / 2;

  /* create a local gaussian window */
  if (wsize != w) {
    free(gaussian);
    gaussian = NULL;
    w = wsize;
  }

  if (!gaussian) /* allocate a gaussian bump */
  {
    float den, norm;

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

  xpix = IMAGEFpix(Ioffset, 0, 0);
  ypix = IMAGEFseq_pix(Ioffset, 0, 0, 1);
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, xpix++, ypix++) {
      /*
        x and y are in window coordinates, while xc and yc are in image
        coordinates.
      */
      /* first build offset direction */
      vx = vy = 0.0f;
      for (g = gaussian, y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0;
        if ((yc < 0) || (yc >= rows)) {
          g += wsize;
          continue;
        }

        for (x = -whalf; x <= whalf; x++, g++) {
          xc = x0 + x;
          if ((xc < 0) || (xc >= cols)) continue;

          fxpix = *IMAGEFpix(Ix, xc, yc);
          fypix = *IMAGEFpix(Iy, xc, yc);
          dot_product = x * fxpix + y * fypix;
          gauss = *g;
          dot_product *= gauss;
          vx += (dot_product * fxpix);
          vy += (dot_product * fypix);
        }
      }

      *xpix = -vx;
      *ypix = -vy;
    }
  }

  return (Ioffset);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              use an offset field to scale an image

              first generate an offset field using the full resolution, but
              only performing the offset calculation on pixel centers of the
              scaled image. Then do length computation using gradient reversal
              in the full image as the criterion.
----------------------------------------------------------------------*/
#define SMOOTH_SIGMA 4.0f
#define WSIZE 3
#define WHALF (WSIZE - 1) / 2

/* global for debugging, make local if this routine is ever really used. */

#define MEDIAN_INDEX ((3 * 3 - 1) / 2)
static int compare_sort_array(const void *pf1, const void *pf2);

IMAGE *ImageOffsetScale(IMAGE *Isrc, IMAGE *Idst)
{
  static IMAGE *Ioffset = NULL, *Ioffset2 = NULL, *Ix = NULL, *Iy = NULL, *Ismooth = NULL;
  static IMAGE *Igauss = NULL;
  int srows, scols, drows, dcols, xs, ys;
  int x0, y0, x, y, ystep, xstep, maxsteps, i, idx, idy;
  float vx, vy, *g, *xpix, *ypix, dx, dy, ox, oy, delta, slope, xf, yf, odx, ody, *dpix, *spix, sort_array[3 * 3],
      *sptr;
  static float *gaussian = NULL;
  float gauss, fxpix, fypix, dot_product;
  int xc, yc;

  srows = Isrc->rows;
  scols = Isrc->cols;
  drows = Idst->rows;
  dcols = Idst->cols;

  xstep = nint((double)scols / (double)dcols);
  ystep = nint((double)srows / (double)drows);

  if (Ioffset && ((Ioffset->rows != drows) || (Ioffset->cols != dcols))) {
    ImageFree(&Ioffset);
    ImageFree(&Ioffset2);
  }
  if (Ix && ((Ix->rows != srows) || (Ix->cols != scols))) {
    ImageFree(&Ix);
    ImageFree(&Iy);
    ImageFree(&Ismooth);
  }

  if (!Ix) {
    Ix = ImageAlloc(srows, scols, PFFLOAT, 1);
    Iy = ImageAlloc(srows, scols, PFFLOAT, 1);
    Ismooth = ImageAlloc(srows, scols, PFFLOAT, 1);
  }
  if (!Ioffset) {
    Ioffset = ImageAlloc(drows, dcols, PFFLOAT, 2);
    Ioffset2 = ImageAlloc(drows, dcols, PFFLOAT, 2);
  }
  if (!Igauss) Igauss = ImageGaussian1d(SMOOTH_SIGMA, 5);

  ImageConvolveGaussian(Isrc, Igauss, Ismooth, 0);
  /* ImageWrite(Ismooth, "smooth.hipl") ;*/
  ImageSobel(Ismooth, NULL, Ix, Iy);

  /* now calculate offset image */
  if (!gaussian) /* allocate a gaussian bump */
  {
    float den, norm;

    gaussian = (float *)calloc(WSIZE * WSIZE, sizeof(float));
    den = WSIZE * WSIZE + WSIZE + 1;
    norm = 0.0f;
    for (g = gaussian, y = 0; y < WSIZE; y++) {
      yc = y - WHALF;
      for (x = 0; x < WSIZE; x++, g++) {
        xc = x - WHALF;
        *g = (float)exp(-36.0 * sqrt((double)(xc * xc + yc * yc)) / (double)den);
        norm += *g;
      }
    }

    /* normalize gaussian */
    for (g = gaussian, y = 0; y < WSIZE; y++) {
      for (x = 0; x < WSIZE; x++, g++) *g /= norm;
    }
  }

  xpix = IMAGEFpix(Ioffset, 0, 0);
  ypix = IMAGEFseq_pix(Ioffset, 0, 0, 1);
  /*
    x0 and y0 are the coordinates of the center of the window in the
    unscaled image.
    xc and yc are the coordinates of the window used for the offset
    computation in the unscaled coordinates.
    x and y are local window coordinates (i.e. -1 .. 1)
  */
  for (y0 = nint(0.5 * (double)ystep); y0 < srows; y0 += ystep) {
    for (x0 = nint(0.5 * (double)xstep); x0 < scols; x0 += xstep, xpix++, ypix++) {
      /*
        x and y are in window coordinates, while xc and yc are in image
        coordinates.
      */
      /* build offset direction */
      vx = vy = 0.0f;
      for (g = gaussian, y = -WHALF; y <= WHALF; y++) {
        /* reflect across the boundary */
        yc = y + y0;
        if ((yc < 0) || (yc >= srows)) {
          g += WSIZE;
          continue;
        }

        for (x = -WHALF; x <= WHALF; x++, g++) {
          xc = x0 + x;
          if ((xc < 0) || (xc >= scols)) continue;

          fxpix = *IMAGEFpix(Ix, xc, yc);
          fypix = *IMAGEFpix(Iy, xc, yc);
          dot_product = x * fxpix + y * fypix;
          gauss = *g;
          dot_product *= gauss;
          vx += (dot_product * fxpix);
          vy += (dot_product * fypix);
        }
      }

      *xpix = -vx;
      *ypix = -vy;
    }
  }


  ImageCopy(Ioffset, Ioffset2);
  /*
    now normalize offset lengths by searching for reversal in gradient field
    in offset direction
   */
  maxsteps = MAX(xstep, ystep) * 2;
  xpix = IMAGEFpix(Ioffset2, 0, 0);
  ypix = IMAGEFseq_pix(Ioffset2, 0, 0, 1);
  for (ys = 0, y0 = nint(0.5 * (double)ystep); y0 < srows; y0 += ystep, ys++) {
    for (xs = 0, x0 = nint(0.5 * (double)xstep); x0 < scols; x0 += xstep, xpix++, ypix++, xs++) {
      ox = *xpix; /* initial offset direction */
      oy = *ypix;

      dx = *IMAGEFpix(Ix, x0, y0); /* starting gradient values */
      dy = *IMAGEFpix(Iy, x0, y0);

#define SMALL 0.00001f

      if ((fabs(ox) < SMALL) && (fabs(oy) < SMALL))
      continue;

      if (fabs(ox) > fabs(oy)) /* use unit steps in x direction */
      {
        delta = nint(ox / (float)fabs(ox));
        slope = delta * oy / ox;
        for (i = 0, yf = (float)y0 + slope, x = x0 + delta; i <= maxsteps; x += delta, yf += slope, i++) {
          y = nint(yf);
          if (y <= 0 || y >= (srows - 1) || x <= 0 || x >= (scols - 1)) break;

          odx = *IMAGEFpix(Ix, x, y);
          ody = *IMAGEFpix(Iy, x, y);
          dot_product = odx * dx + ody * dy;
          if (dot_product <= 0) break;
        }
        x -= delta;
        y = nint(yf - slope);
      }
      else /* use unit steps in y direction */
      {
        delta = nint(oy / (float)fabs(oy));
        slope = delta * ox / oy;
        for (i = 0, xf = (float)x0 + slope, y = y0 + delta; i < maxsteps; y += delta, xf += slope, i++) {
          x = nint(xf);
          if (y <= 0 || y >= (srows - 1) || x <= 0 || x >= (scols - 1)) break;

          odx = *IMAGEFpix(Ix, x, y);
          ody = *IMAGEFpix(Iy, x, y);
          dot_product = odx * dx + ody * dy;
          if (dot_product <= 0) break;
        }
        y -= delta;
        x = nint(xf - slope);
      }

      *xpix = x - x0;
      *ypix = y - y0;
    }
  }


  /* now use the offset field to scale the image with a median filter */
  xpix = IMAGEFpix(Ioffset2, 0, 0);
  ypix = IMAGEFseq_pix(Ioffset2, 0, 0, 1);
  dpix = IMAGEFpix(Idst, 0, 0);

  /* apply median filter */
  for (y0 = nint(0.5 * (double)ystep); y0 < srows; y0 += ystep) {
    for (x0 = nint(0.5 * (double)xstep); x0 < scols; x0 += xstep, xpix++, ypix++) {
      idx = nint(*xpix); /*  offset direction */
      idy = nint(*ypix);

      for (sptr = sort_array, y = -WHALF; y <= WHALF; y++) {
        yc = y + y0 + idy;
        if (yc < 0)
          yc = 0;
        else if (yc >= srows)
          yc = srows - 1;

        spix = IMAGEFpix(Isrc, 0, yc);
        for (x = -WHALF; x <= WHALF; x++) {
          xc = x0 + x + idx;
          if (xc < 0)
            xc = 0;
          else if (xc >= scols)
            xc = scols - 1;
          *sptr++ = *(spix + xc);
        }
      }
      qsort(sort_array, 3 * 3, sizeof(float), compare_sort_array);
      *dpix++ = sort_array[MEDIAN_INDEX];
    }
  }



  return (Idst);
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
IMAGE *ImageOffsetOrientation(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient)
{
  if (!Iorient) Iorient = ImageAlloc(Ix->rows, Ix->cols, PFFLOAT, 2);
  Iorient->num_frame = 1;
  ImageMeanFilter(Ix, 3, Iorient);
  Iorient->image += Iorient->sizeimage;
  ImageMeanFilter(Iy, 3, Iorient);
  Iorient->image -= Iorient->sizeimage;
  Iorient->num_frame = 2;

  return (Iorient);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient, IMAGE *Ioffset)
{
  int x0, y0, rows, cols, x, y, whalf, xc, yc, yoff, off;
  float *xpix, *ypix, dx, dy, *or_xpix, *or_ypix, *oxpix, *oypix, dir, ox, oy, dot;

  rows = Ix->rows;
  cols = Ix->cols;

  if (!Ioffset) Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);

  if (!ImageCheckSize(Ix, Ioffset, 0, 0, 2)) {
    ImageFree(&Ioffset);
    Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);
  }

  whalf = (wsize - 1) / 2;
  xpix = IMAGEFpix(Ix, 0, 0);
  ypix = IMAGEFpix(Iy, 0, 0);
  or_xpix = IMAGEFpix(Iorient, 0, 0);
  or_ypix = IMAGEFseq_pix(Iorient, 0, 0, 1);
  oxpix = IMAGEFpix(Ioffset, 0, 0);
  oypix = IMAGEFseq_pix(Ioffset, 0, 0, 1);
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, xpix++, ypix++) {
      /*
        Now calculate the orientation for this point by averaging local gradient
        orientation within the specified window.

        x and y are in window coordinates, while xc and yc are in image
        coordinates.
      */
      /* calculate orientation vector */
      ox = *or_xpix++;
      oy = *or_ypix++;
      dir = 0.0f;
      for (y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0;
        if ((yc < 0) || (yc >= rows)) continue;

        yoff = y * cols;
        for (x = -whalf; x <= whalf; x++) {
          xc = x0 + x;
          if ((xc < 0) || (xc >= cols)) continue;

          off = yoff + x;
          dx = *(xpix + off);
          dy = *(ypix + off);
          dot = dx * ox + dy * oy;
          if (dot < 0.0f) dot = 0.0f;
          dir += (x * ox + y * oy) * dot;
        }
      }

      if (FZERO(dir))
        ox = oy = 0.0f;
      else if (dir > 0.0f) /* flip by 180 */
      {
        ox = -ox;
        oy = -oy;
      }
      *oxpix++ = ox;
      *oypix++ = oy;
    }
  }
  return (Ioffset);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *ImageOffsetDirectionMap(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient, IMAGE *Idir, IMAGE *Ioffset)
{
  int x0, y0, rows, cols, x, y, whalf, xc, yc, yoff, off;
  float *xpix, *ypix, dx, dy, *or_xpix, *or_ypix, *oxpix, *oypix, dir, ox, oy, dot;

  rows = Ix->rows;
  cols = Ix->cols;

  if (!Ioffset) Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);

  if (!ImageCheckSize(Ix, Ioffset, 0, 0, 2)) {
    ImageFree(&Ioffset);
    Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);
  }

  whalf = (wsize - 1) / 2;
  xpix = IMAGEFpix(Ix, 0, 0);
  ypix = IMAGEFpix(Iy, 0, 0);
  or_xpix = IMAGEFpix(Iorient, 0, 0);
  or_ypix = IMAGEFseq_pix(Iorient, 0, 0, 1);
  oxpix = IMAGEFpix(Ioffset, 0, 0);
  oypix = IMAGEFseq_pix(Ioffset, 0, 0, 1);
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, xpix++, ypix++) {
      /*
        Now calculate the direction for this point by averaging local gradient
        orientation within the specified window.

        x and y are in window coordinates, while xc and yc are in image
        coordinates.
      */
      /* calculate orientation vector */
      ox = *or_xpix++;
      oy = *or_ypix++;
      dir = 0.0f;
      for (y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0;
        if ((yc < 0) || (yc >= rows)) continue;

        yoff = y * cols;
        for (x = -whalf; x <= whalf; x++) {
          xc = x0 + x;
          if ((xc < 0) || (xc >= cols)) continue;

          off = yoff + x;
          dx = *(xpix + off);
          dy = *(ypix + off);
          dot = dx * ox + dy * oy;
          if (dot < 0.0f) dot = 0.0f;
          dir += (x * ox + y * oy) * dot;
        }
      }

      if (FZERO(dir))
        *IMAGEFpix(Idir, x0, y0) = 0.0f;
      else if (dir > 0.0f)
        *IMAGEFpix(Idir, x0, y0) = -1.0f;
      else
        *IMAGEFpix(Idir, x0, y0) = 1.0f;

      if (FZERO(dir))
        ox = oy = 0.0f;
      else if (dir > 0.0f) /* flip by 180 */
      {
        ox = -ox;
        oy = -oy;
      }
      *oxpix++ = ox;
      *oypix++ = oy;
    }
  }
  return (Ioffset);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             use a Bresenham line drawing algorithm to do search
----------------------------------------------------------------------*/
#define FSCALE 100.0f

IMAGE *ImageOffsetMagnitude(IMAGE *Isrc, IMAGE *Idst, int maxsteps)
{
  int rows, cols, x, y, ax, ay, sx, sy, x1, y1, dx, dy, odx, ody, d, xn, yn, steps, dot = 0, xold, yold;
  float *src_xpix, *src_ypix, *dst_xpix, *dst_ypix, *oxpix, *oypix;

  if (!Idst) Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, Isrc->num_frame);

  rows = Isrc->rows;
  cols = Isrc->cols;

  src_xpix = IMAGEFpix(Isrc, 0, 0);
  src_ypix = IMAGEFseq_pix(Isrc, 0, 0, 1);
  dst_xpix = IMAGEFpix(Idst, 0, 0);
  dst_ypix = IMAGEFseq_pix(Idst, 0, 0, 1);
  for (y1 = 0; y1 < rows; y1++) {
    for (x1 = 0; x1 < cols; x1++) {
      /* do a Bresenham algorithm do find the offset line at this point */
      dx = nint(*src_xpix * FSCALE);
      dy = nint(*src_ypix * FSCALE);
      xold = x = x1;
      yold = y = y1;
      ax = ABS(dx) << 1;
      sx = SGN(dx);
      ay = ABS(dy) << 1;
      sy = SGN(dy);

      oxpix = src_xpix++;
      oypix = src_ypix++;

      if (ax > ay) /* x dominant */
      {
        d = ay - (ax >> 1);
        for (steps = 0; steps < maxsteps; steps++) {
          odx = nint(*oxpix * FSCALE);
          ody = nint(*oypix * FSCALE);
          dot = odx * dx + ody * dy;

          if (dot <= 0) break;
          if (d >= 0) /* move in y direction */
          {
            yn = y + sy;
            if (yn < 0 || yn >= rows) break;
            xn = x;
            oxpix += (sy * cols);
            oypix += (sy * cols);
            d -= ax;
          }
          else /* move in x direction */
          {
            xn = x + sx;
            if (xn < 0 || xn >= cols) break;
            yn = y;
            oxpix += sx;
            oypix += sx;
            d += ay;
          }

          xold = x;
          yold = y;
          x = xn;
          y = yn;
        }
      }
      else /* y dominant */
      {
        d = ax - (ay >> 1);
        for (steps = 0; steps < maxsteps; steps++) {
          odx = nint(*oxpix * FSCALE);
          ody = nint(*oypix * FSCALE); /* offset vector at this point */
          dot = odx * dx + ody * dy;
          if (dot <= 0) /* vector field has reversed or changed directions */
            break;
          if (d >= 0) /* move only in x direction */
          {
            xn = x + sx;
            if (xn < 0 || xn >= cols) break;
            yn = y;
            oxpix += sx;
            oypix += sx;
            d -= ay;
          }
          else /* only move in y direction */
          {
            yn = y + sy;
            if (yn < 0 || yn >= rows) break;
            xn = x;
            oypix += (sy * cols);
            oxpix += (sy * cols);
            d += ax;
          }

          xold = x;
          yold = y;
          x = xn;
          y = yn;
        }
      }

      *dst_xpix++ = (float)(xold - x1);
      *dst_ypix++ = (float)(yold - y1);
    }
  }

  return (Idst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              calculate direction using Nitzberg-Shiota formula
----------------------------------------------------------------------*/
IMAGE *ImageNitshiOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient, IMAGE *Ioffset)
{
  int x0, y0, rows, cols, x, y, whalf, xc, yc, yoff, off;
  float *xpix, *ypix, dx, dy, *or_xpix, *or_ypix, *oxpix, *oypix, dir, dot, ox, oy, dirx, diry;

  rows = Ix->rows;
  cols = Ix->cols;

  if (!Ioffset) Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);

  if (!ImageCheckSize(Ix, Ioffset, 0, 0, 2)) {
    ImageFree(&Ioffset);
    Ioffset = ImageAlloc(rows, cols, PFFLOAT, 2);
  }

  whalf = (wsize - 1) / 2;
  xpix = IMAGEFpix(Ix, 0, 0);
  ypix = IMAGEFpix(Iy, 0, 0);
  or_xpix = IMAGEFpix(Iorient, 0, 0);
  or_ypix = IMAGEFseq_pix(Iorient, 0, 0, 1);
  oxpix = IMAGEFpix(Ioffset, 0, 0);
  oypix = IMAGEFseq_pix(Ioffset, 0, 0, 1);
  for (y0 = 0; y0 < rows; y0++) {
    for (x0 = 0; x0 < cols; x0++, xpix++, ypix++) {
      /*
        Now calculate the orientation for this point by averaging local gradient
        orientation within the specified window.

        x and y are in window coordinates, while xc and yc are in image
        coordinates.
      */
      /* calculate orientation vector */
      ox = *or_xpix++;
      oy = *or_ypix++;

      dirx = diry = dir = 0.0f;

      for (y = -whalf; y <= whalf; y++) {
        /* reflect across the boundary */
        yc = y + y0;
        if ((yc < 0) || (yc >= rows)) continue;

        yoff = y * cols;
        for (x = -whalf; x <= whalf; x++) {
          xc = x0 + x;
          if ((xc < 0) || (xc >= cols)) continue;

          off = yoff + x;
          dx = *(xpix + off);
          dy = *(ypix + off);
          dot = (x * dx + y * dy);
          dirx += dot * dx;
          diry += dot * dy;
        }
      }
      dir = dirx * ox + diry * oy;

      if (ISTINY(dir))
        ox = oy = 0.0f;
      else if (dir > 0.0f) /* flip by 180 */
      {
        ox = -ox;
        oy = -oy;
      }
      *oxpix++ = ox;
      *oypix++ = oy;
    }
  }
  return (Ioffset);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
              only calculate directions at points along search vector.
----------------------------------------------------------------------*/

static int imageOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient, int x0, int y0);
static int imageOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Iorient, int x0, int y0)
{
  int rows, cols, x, y, whalf, xc, yc, yoff, off, d;
  float *xpix, *ypix, dx, dy, *or_xpix, *or_ypix, dir, ox, oy;

  rows = Ix->rows;
  cols = Ix->cols;

  whalf = (wsize - 1) / 2;
  xpix = IMAGEFpix(Ix, x0, y0);
  ypix = IMAGEFpix(Iy, x0, y0);
  or_xpix = IMAGEFpix(Iorient, x0, y0);
  or_ypix = IMAGEFseq_pix(Iorient, x0, y0, 1);

  /*
    Now calculate the orientation for this point by averaging local gradient
    orientation within the specified window.

    x and y are in window coordinates, while xc and yc are in image
    coordinates.
  */
  /* calculate orientation vector */
  ox = *or_xpix++;
  oy = *or_ypix++;

  dir = 0.0f;
  for (y = -whalf; y <= whalf; y++) {
    /* reflect across the boundary */
    yc = y + y0;
    if ((yc < 0) || (yc >= rows)) continue;

    yoff = y * cols;
    for (x = -whalf; x <= whalf; x++) {
      xc = x0 + x;
      if ((xc < 0) || (xc >= cols)) continue;

      off = yoff + x;
      dx = *(xpix + off);
      dy = *(ypix + off);
      dir += (x * ox + y * oy) * fabs(dx * ox + dy * oy);
    }
  }

  if (FZERO(dir))
    d = 0;
  else if (dir > 0.0f) /* flip by 180 */
    d = -1;
  else
    d = 1;

  return (d);
}

IMAGE *ImageOffsetDirectionMagnitude(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy, int wsize, IMAGE *Idst, int maxsteps)
{
  int rows, cols, x, y, ax, ay, sx, sy, x1, y1, dx, dy, odx, ody, d, xn, yn, steps, dir, dot;
  float *src_xpix, *src_ypix, *dst_xpix, *dst_ypix, *oxpix, *oypix, fdir;
  byte *calculated;
  static IMAGE *Icalculated = NULL;

  rows = Isrc->rows;
  cols = Isrc->cols;

  if (!Icalculated || (rows != Icalculated->rows || cols != Icalculated->cols)) {
    if (Icalculated) ImageFree(&Icalculated);
    Icalculated = ImageAlloc(rows, cols, PFBYTE, 1);
  }
  else
    ImageClear(Icalculated);

  if (!Idst) Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, Isrc->num_frame);

  src_xpix = IMAGEFpix(Isrc, 0, 0);
  src_ypix = IMAGEFseq_pix(Isrc, 0, 0, 1);
  dst_xpix = IMAGEFpix(Idst, 0, 0);
  dst_ypix = IMAGEFseq_pix(Idst, 0, 0, 1);
  calculated = IMAGEpix(Icalculated, 0, 0);
  for (y1 = 0; y1 < rows; y1++) {
    for (x1 = 0; x1 < cols; x1++, calculated++) {
      /* do a Bresenham algorithm do find the offset line at this point */
      if (*calculated == 0) {
        dir = imageOffsetDirection(Ix, Iy, wsize, Isrc, x1, y1);
        fdir = (float)dir;
        *calculated = 1;
        *IMAGEFpix(Idst, x1, y1) = *src_xpix * fdir;
        *IMAGEFseq_pix(Idst, x1, y1, 1) = *src_ypix * fdir;
      }
      dx = nint(*IMAGEFpix(Idst, x1, y1) * FSCALE);
      dy = nint(*IMAGEFseq_pix(Idst, x1, y1, 1) * FSCALE);
      x = x1;
      y = y1;
      ax = ABS(dx) << 1;
      sx = SGN(dx);
      ay = ABS(dy) << 1;
      sy = SGN(dy);

      oxpix = src_xpix++;
      oypix = src_ypix++;

      if (ax > ay) /* x dominant */
      {
        d = ay - (ax >> 1);
        for (steps = 0; steps < maxsteps; steps++) {
          if (!*IMAGEpix(Icalculated, x, y)) {
            dir = imageOffsetDirection(Ix, Iy, wsize, Isrc, x, y);
            fdir = (float)dir;
            *IMAGEpix(Icalculated, x, y) = 1;
            *IMAGEFpix(Idst, x, y) = *oxpix * fdir;
            *IMAGEFseq_pix(Idst, x, y, 1) = *oypix * fdir;
          }
          odx = nint(*IMAGEFpix(Idst, x, y) * FSCALE);
          ody = nint(*IMAGEFseq_pix(Idst, x, y, 1) * FSCALE);
          dot = odx * dx + ody * dy;
          if (dot <= 0) break;
          if (d >= 0) {
            yn = y + sy;
            if (yn < 0 || yn >= rows) break;
            oxpix += (sy * cols);
            oypix += (sy * cols);
            d -= ax;
          }
          else
            yn = y;
          oxpix += sx;
          oypix += sx;
          xn = x + sx;
          if (xn < 0 || xn >= cols) break;

          x = xn;
          y = yn;
          d += ay;
        }
      }
      else /* y dominant */
      {
        d = ax - (ay >> 1);
        for (steps = 0; steps < maxsteps; steps++) {
          if (!*IMAGEpix(Icalculated, x, y)) {
            dir = imageOffsetDirection(Ix, Iy, wsize, Isrc, x, y);
            fdir = (float)dir;
            *IMAGEpix(Icalculated, x, y) = 1;
            *IMAGEFpix(Idst, x, y) = *oxpix * fdir;
            *IMAGEFseq_pix(Idst, x, y, 1) = *oypix * fdir;
          }
          odx = nint(*IMAGEFpix(Idst, x, y) * FSCALE);
          ody = nint(*IMAGEFseq_pix(Idst, x, y, 1) * FSCALE);
          dot = odx * dx + ody * dy;
          if (dot <= 0) break;
          if (d >= 0) {
            xn = x + sx;
            if (xn < 0 || xn >= cols) break;
            oxpix += sx;
            oypix += sx;
            d -= ay;
          }
          else
            xn = x;
          yn = y + sy;
          if (yn < 0 || yn >= rows) break;

          x = xn;
          y = yn;
          oypix += (sy * cols);
          oxpix += (sy * cols);
          d += ax;
        }
      }
      *dst_xpix++ = (float)(x - x1);
      *dst_ypix++ = (float)(y - y1);
    }
  }

  return (Idst);
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
              apply an offset vector to a filtered image.
----------------------------------------------------------------------*/
IMAGE *ImageFilterMinMax(IMAGE *Imin, IMAGE *Imax, IMAGE *Idir, IMAGE *Ioffset, IMAGE *Idst)
{
  int x, y, rows, cols, dx, dy, dir;
  float *dst, src, *dx_pix, *dy_pix, *dir_pix;
  IMAGE *Iout;

  rows = Imin->rows;
  cols = Imin->cols;

  if (!Idst) Idst = ImageAlloc(rows, cols, PFFLOAT, 1);

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1);
  else

    Iout = Idst;
  if (!ImageCheckSize(Imin, Idst, 0, 0, 0)) ErrorReturn(NULL, (ERROR_SIZE, "ImageApplyOffset: dst not big enough"));

  dst = IMAGEFpix(Iout, 0, 0);
  dir_pix = IMAGEFpix(Idir, 0, 0);
  dx_pix = IMAGEFpix(Ioffset, 0, 0);
  dy_pix = IMAGEFseq_pix(Ioffset, 0, 0, 1);

  for (y = 0; y < rows; y++) {
    for (x = 0; x < cols; x++) {
      if (x == 7 && y == 50) DiagBreak();

      dir = (int)*dir_pix++;
      dx = (int)*dx_pix++;
      dy = (int)*dy_pix++;
      if (dir > 0) /* moving in gradient direction */
        src = *IMAGEFpix(Imax, x + dx, y + dy);
      else if (dir < 0)
        src = *IMAGEFpix(Imin, x + dx, y + dy);
      else
        src = (*IMAGEFpix(Imin, x + dx, y + dy) + *IMAGEFpix(Imax, x + dx, y + dy)) / 2;
      *dst++ = src;
    }
  }

  if (Iout != Idst) {
    ImageCopy(Iout, Idst);
    ImageFree(&Iout);
  }
  return (Idst);
}
