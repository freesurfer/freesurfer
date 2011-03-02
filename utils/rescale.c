/**
 * @file  rescale.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.3 $
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


/*******************************************************************
 * Name:    rescale.c
 * Author:  Douglas N. Greve, 5/14/96
 * Purpose: rescales the pixels of an image to within the given
 *          maximum and minimum.  As a bonus, it also returns
 *          the row and column of the min and max.
 *
 *          Input  Formats: float
 *          Output Formats: byte
 *
 * Notes:   1. The input and output do not have to be of the same
 *          format.
 *          2. MinPoint[0] = row ...
 *          3. This is different than h_linscale() in two ways:
 *             a.  The desired output range is specified instead
 *                 of the scale and offset.
 *             b.  The output image does not have to be float
 *                 (as in linscale).
 ******************************************************************/
#include <hips.h>
#include "minmaxrc.h"
#include "rescale.h"
#include "macros.h"

/* y = Offset + Slope * x, y=New, x=Old  ***/
static int h_GetScaleParams(float OldMin,  float OldMax,
                            float NewMin,  float NewMax,
                            float *Offset, float *Slope)
{
  float Del;
  int Warning=1;

  Del = OldMax - OldMin;
  if (Del < .00000001F)
  {
    *Offset = 0.0F;
    *Slope  = 1.0F;
    Warning = 1;
  }
  else
  {
    *Slope = (NewMax - NewMin)/Del;
    *Offset = (NewMin - *Slope * OldMin);
    Warning = 0;
  }

  return(Warning);
}
/*************************************************************/
int h_rescale(struct header *phdSrc,
              float NewMin, float NewMax,
              int MinPoint[2], int MaxPoint[2],
              struct header *phdDst)
{
  Pixelval Min,Max;
  float Offset, Slope;
  int err;

  /** find the present min and max and where they are **/
  err = h_minmaxrc(phdSrc, &Min, MinPoint, &Max, MaxPoint);
  if (err) return(err);

  /****** switch over the input format *************/
  switch (phdSrc->pixel_format)
  {

  case PFFLOAT:
    /* for float input, switch over the output format */
    switch (phdDst->pixel_format)
    {

    case PFBYTE:
      h_GetScaleParams(Min.v_float, Max.v_float,
                       NewMin, NewMax, &Offset, &Slope);
      return(h_rescale_fb(phdSrc, Slope, Offset, phdDst));
      break;

    case PFFLOAT:
      h_GetScaleParams(Min.v_float, Max.v_float,
                       NewMin, NewMax, &Offset, &Slope);
      return(h_rescale_ff(phdSrc, Slope, Offset, phdDst));
      break;

    default:
      return(perr(HE_FMTSUBR,"h_rescale: Destination",
                  hformatname(phdDst->pixel_format)));
    }/* end switch over destination pixel format for float input*/
    break;

  default:
    return(perr(HE_FMTSUBR,"h_rescale: Source",
                hformatname(phdSrc->pixel_format)));
  }/********** end switch over source pixel format *************/
  return(HIPS_OK) ;
}

/********************************************************/
int h_rescale_fb(struct header *phdSrc,
                 float Slope, float Offset,
                 struct header *phdDst)
{
  register INT32B n;
  register float *pSrc;
  register byte  *pDst;

  pSrc = (float *) phdSrc->image;
  pDst = (byte  *) phdDst->image;

  for (n=0;n<phdSrc->numpix;n++)
    *pDst++ = (byte)( *pSrc++ * Slope + Offset);

  return(HIPS_OK);
}

/********************************************************/
int h_rescale_ff(struct header *phdSrc,
                 float Slope, float Offset,
                 struct header *phdDst)
{
  register INT32B n;
  register float *pSrc;
  register float  *pDst;

  pSrc = (float *) phdSrc->image;
  pDst = (float  *) phdDst->image;

  for (n=0;n<phdSrc->numpix;n++)
    *pDst++ = (float)( *pSrc++ * Slope + Offset);

  return(HIPS_OK);
}

