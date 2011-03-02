/**
 * @file  kinput.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
 *    $Revision: 1.8 $
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hips.h"
#include "kinput.h"
#include "macros.h"

/* debugging stuff */
#define TEST    0

#if USE_PYRAMID
static float kinputGetPoint(HIPSPyramid *pyr, int x, int y, int level) ;
#endif


#define SIGMA_SCALE  1.5f

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
KINPUT *
KinputAlloc(int rows, int cols, int nscales, int input_size, float sigma,
            float sigma_scale_factor, int abs_gradient)
{
  KINPUT  *kinput ;
  int     scale ;

  if (nscales > MAX_SCALES)
    nscales = MAX_SCALES ;

  kinput = (KINPUT *)calloc(1, sizeof(KINPUT)) ;
  kinput->nscales = nscales ;
  kinput->ninputs = nscales * input_size * input_size * 2 ;
  kinput->parms.abs_gradient = abs_gradient ;
  if (FZERO(sigma_scale_factor))
    sigma_scale_factor = kinput->parms.sigma_scale_factor = SIGMA_SCALE_FACTOR;
  else
    kinput->parms.sigma_scale_factor = sigma_scale_factor ;

#if !USE_PYRAMID
  kinput->xInputs = ImageAlloc(rows, cols, PFFLOAT, nscales) ;
  kinput->yInputs = ImageAlloc(rows, cols, PFFLOAT, nscales) ;
  for (scale = 1 ; scale < nscales ; scale++)
  {
    kinput->gImages[scale] =
      ImageGaussian1d(sigma*sigma_scale_factor*(float)(scale), 0) ;
    kinput->parms.sigmas[scale] = (float)scale * sigma_scale_factor * sigma ;
  }
#endif

  kinput->inputs = (float *)calloc(kinput->ninputs, sizeof(float)) ;

  /* fill out parameter block */
  kinput->parms.nscales = nscales ;
  kinput->parms.input_size = input_size ;

  return(kinput) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
KinputFree(KINPUT **pkinput)
{
  KINPUT  *kinput ;
  int     scale ;

  kinput = *pkinput ;
  *pkinput = NULL ;

#if USE_PYRAMID
  ImageFreePyramid(kinput->xpyr) ;
  ImageFreePyramid(kinput->ypyr) ;
  if (kinput->xImage)
    ImageFree(&kinput->xImage) ;
  if (kinput->yImage)
    ImageFree(&kinput->yImage) ;
#else
  ImageFree(&kinput->xInputs) ;
  ImageFree(&kinput->yInputs) ;
  for (scale = 1 ; scale < kinput->nscales ; scale++)
    ImageFree(&kinput->gImages[scale]) ;
#endif

  free(kinput->inputs) ;

  free(kinput) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
KinputInit(KINPUT *kinput, IMAGE *image)
{
  int scale ;

#if USE_PYRAMID
  /* first compute derivatives */
  kinput->xImage = ImageXDerivative(image, kinput->xImage) ;
  kinput->yImage = ImageYDerivative(image, kinput->yImage) ;
  kinput->xpyr = ImagePyramid(kinput->xImage, kinput->xpyr, kinput->nscales) ;
  kinput->ypyr = ImagePyramid(kinput->yImage, kinput->ypyr, kinput->nscales) ;
#else

  /* 1st image in inputs is derivative at inner scale */
  ImageXDerivative(image, kinput->xInputs) ;
  ImageYDerivative(image, kinput->yInputs) ;

  /* now compute blurred derivates at different scales */
  for (scale = 1 ; scale < kinput->nscales ; scale++)
  {
    ImageConvolveGaussian(kinput->xInputs, kinput->gImages[scale],
                          kinput->xInputs, scale) ;
    ImageConvolveGaussian(kinput->yInputs, kinput->gImages[scale],
                          kinput->yInputs, scale) ;
  }
#if 0
  ImageWrite(image, "i.hipl") ;
  ImageWriteFrames(kinput->xInputs, "xin.hipl", 0, kinput->xInputs->num_frame) ;
  ImageWriteFrames(kinput->yInputs, "yin.hipl", 0, kinput->yInputs->num_frame) ;
#endif
#endif
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
KinputVector(KINPUT *kinput, int x0, int y0)
{
  int    x, y, x1, y1, rows, cols, scale, neg_half, half_in ;
  float  *in, xval, yval ;

  in = kinput->inputs ;
  neg_half = (kinput->parms.input_size - 1) / 2 ;
  half_in = kinput->parms.input_size / 2 ;

#if USE_PYRAMID
  /* int index ; */
  rows = kinput->xImage->rows ;
  cols = kinput->xImage->cols ;

  /*
    the input vector consists of pairs of x and y derivative values at each
    level of the pyramid at this point.
  */
  for (index = scale = 0 ; scale < kinput->nscales ; scale++)
  {
    for (y = y0-(neg_half ; y <= y0+half_in ; y++)
           {
             for (x = x0-(neg_half ; x <= x0+half_in ; x++)
                      {
                        xval = kinputGetPoint(kinput->xpyr, x, y, scale) ;
                          if (kinput->parms.abs_gradient)
                            xval = fabs(xval) ;
                          kinput->inputs[index++] = xval ;
                          yval = kinputGetPoint(kinput->ypyr, x, y, scale) ;
                          if (kinput->parms.abs_gradient)
                            yval = fabs(yval) ;
                          kinput->inputs[index++] = yval ;
                        }
                      }
                    }
#else
  /*
  the input vector consists of pairs of x and y derivative values at each
  level of blurring at this point.
  */
  rows = kinput->xInputs->rows ;
  cols = kinput->xInputs->cols ;

  for (scale = 0 ; scale < kinput->nscales ; scale++)
  {
    for (y1 = -neg_half; y1 <= half_in ; y1++)
    {
      for (x1 = -neg_half ; x1 <= half_in ; x1++)
      {
        x = x0 + x1 ;
        if (x < 0)
          x = 0 ;
        if (x >= cols)
          x = cols - 1 ;

        y = y0 + y1 ;
        if (y < 0)
          y = 0 ;
        if (y >= rows)
          y = rows - 1 ;

        xval = *IMAGEFseq_pix(kinput->xInputs, x, y, scale) ;
        yval = *IMAGEFseq_pix(kinput->yInputs, x, y, scale) ;

#define DONT_USE_INNER_SCALE 0
#if DONT_USE_INNER_SCALE
        if (!scale)
          xval = yval = 0.0f ;
#endif

        *in++ = xval ;
        *in++ = yval ;
      }
    }
  }
#endif

                    return(0) ;
}
#if USE_PYRAMID

static float
kinputGetPoint(IMAGEPyramid *pyr, int x, int y, int level)
{
  float val;

  int   x1, x2, y1, y2, cols, rows, cols0, rows0 ;
  float dx, dy, xc, yc, *image, *f1p, *f2p, *f3p, *f4p, xpct, ypct;

  cols0 = pyr->images[0]->cols ;
  rows0 = pyr->images[0]->rows ;
  if (x < 0)
    x = 0 ;
  if (y < 0)
    y = 0 ;
  if (x >= cols0)
    x = cols0-1 ;
  if (y >= rows0)
    y = rows0-1 ;

  image = IMAGEFpix(pyr->images[level], 0, 0) ;
  cols = pyr->images[level]->cols ;
  rows = pyr->images[level]->rows ;

  xpct  = (float)x / (float)(cols0-1) ;
  ypct  = (float)y / (float)(rows0-1) ;
  xc = xpct * (float)(cols-1) ;
  yc = ypct * (float)(rows-1) ;

#define BILINEAR 1
#if BILINEAR

  /* do bilinear interpolation on 4 surrounding points */
  x1 = (int)floor(xc) ;
  x2 = (int)ceil(xc) ;
  y1 = (int)floor(yc) ;
  y2 = (int)ceil(yc) ;

#if 0
  if (x2 >= cols)
    x2 = cols-1 ;
  if (y2 >= rows)
    y2 = rows - 1 ;
#endif

  if ((x1 < 0) || (y1 < 0) || (x2 >= cols) ||
      (y2 >= pyr->images[level]->rows))
  {
    fprintf(stderr,
            "kinputGetPoint(%d, %d, %d): point (%d,%d,%d,%d) out of bounds\n",
            x, y, level, x1, y1, x2, y2) ;
    return(0.0f) ;
  }

  /* all other distances are 1-a or 1-b */
  dy = yc - (float)y1 ;
  dx = xc - (float)x1 ;

  f1p = image + y1 * cols + x1 ;    /* dy, dx */
  f2p = image + y2 * cols + x1 ;    /* 1-dy, dx */
  f3p = image + y1 * cols + x2 ;    /* dy, 1-dx */
  f4p = image + y2 * cols + x2 ;    /* 1-dy, 1-dx */
  val = (1.0 - dy) * ((1.0 - dx) * *f1p + dx * *f3p) +
        dy * ((1.0 - dx) * *f2p + dx * *f4p) ;

#else
  x1 = nint(xc) ;
  y1 = nint(yc) ;
  val = *IMAGEFpix(pyr->images[level], x1, y1) ;
#endif
#if TEST
  if (x == TCOL && y == TROW)
    fprintf(stderr, "inputs(%2.1f, %2.1f) = %2.4f\n", xc, yc, val) ;
#endif
  return(val) ;
}
#endif
