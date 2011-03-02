/**
 * @file  minmaxrc.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:45 $
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


/**************************************************************
 * Name:    minmaxrc.c
 * Author:  Douglas N. Greve, 5/14/96
 * Purpose: hips-based routines to find the minimum and maximum
 *          pixel values in an image as well as their row and
 *          column.  The hips function h_minmax() does not
 *          return the row and col of the min and max.
 * Notes:   1. MinPoint[0] = row of minimum...
 ***************************************************************/

#include <hips.h>
#include "minmaxrc.h"

/************************************************************/
int h_minmaxrc(struct header *phdSrc,
               Pixelval *Min, int MinPoint[2],
               Pixelval *Max, int MaxPoint[2])
{
  switch (phdSrc->pixel_format)
  {
  case PFBYTE:
    return(h_minmaxrc_b(phdSrc,(byte *)   Min, MinPoint,
                        (byte *)   Max, MaxPoint));
  case PFSHORT:
    return(h_minmaxrc_s(phdSrc, (short *) Min, MinPoint,
                        (short *) Max, MaxPoint));
  case PFINT:
    return(h_minmaxrc_i(phdSrc, (int *)    Min, MinPoint,
                        (int * )  Max, MaxPoint));
  case PFFLOAT:
    return(h_minmaxrc_f(phdSrc, (float *) Min, MinPoint,
                        (float *) Max, MaxPoint));
  default:
    return(perr(HE_FMTSUBR,"h_minmaxrc:",
                hformatname(phdSrc->pixel_format)));
  }
}

/************************************************************/
int h_minmaxrc_b(struct header *phdSrc,
                 byte *Min, int MinPoint[2],
                 byte *Max, int MaxPoint[2])
{
  register int r,c;
  register byte x,*p;

  p = phdSrc->image;

  *Min = *p;
  *Max = *p;
  MinPoint[0] = 0;
  MinPoint[1] = 0;
  MaxPoint[0] = 0;
  MaxPoint[1] = 0;

  for (r=0;r<phdSrc->rows;r++)
  {
    for (c=0;c<phdSrc->cols;c++)
    {
      x = *p++;
      if (x > *Max)
      {
        *Max = x;
        MaxPoint[0] = r;
        MaxPoint[1] = c;
      }
      else if (x < *Min)
      {
        *Min = x;
        MinPoint[0] = r;
        MinPoint[1] = c;
      }
    }
  }
  return(HIPS_OK);
}
/************************************************************/
int h_minmaxrc_s(struct header *phdSrc,
                 short *Min, int MinPoint[2],
                 short *Max, int MaxPoint[2])
{
  register int r,c;
  register short x,*p;

  p = (short *)phdSrc->image;

  *Min = *p;
  *Max = *p;
  MinPoint[0] = 0;
  MinPoint[1] = 0;
  MaxPoint[0] = 0;
  MaxPoint[1] = 0;

  for (r=0;r<phdSrc->rows;r++)
  {
    for (c=0;c<phdSrc->cols;c++)
    {
      x = *p++;
      if (x > *Max)
      {
        *Max = x;
        MaxPoint[0] = r;
        MaxPoint[1] = c;
      }
      else if (x < *Min)
      {
        *Min = x;
        MinPoint[0] = r;
        MinPoint[1] = c;
      }
    }
  }
  return(HIPS_OK);
}
/************************************************************/
int h_minmaxrc_i(struct header *phdSrc,
                 int *Min, int MinPoint[2],
                 int *Max, int MaxPoint[2])
{
  register int r,c;
  register int x,*p;

  p = (int *)phdSrc->image;

  *Min = *p;
  *Max = *p;
  MinPoint[0] = 0;
  MinPoint[1] = 0;
  MaxPoint[0] = 0;
  MaxPoint[1] = 0;

  for (r=0;r<phdSrc->rows;r++)
  {
    for (c=0;c<phdSrc->cols;c++)
    {
      x = *p++;
      if (x > *Max)
      {
        *Max = x;
        MaxPoint[0] = r;
        MaxPoint[1] = c;
      }
      else if (x < *Min)
      {
        *Min = x;
        MinPoint[0] = r;
        MinPoint[1] = c;
      }
    }
  }
  return(HIPS_OK);
}
/************************************************************/
int h_minmaxrc_f(struct header *phdSrc,
                 float *Min, int MinPoint[2],
                 float *Max, int MaxPoint[2])
{
  register int r,c;
  register float x,*p;

  p = (float *)phdSrc->image;

  *Min = *p;
  *Max = *p;
  MinPoint[0] = 0;
  MinPoint[1] = 0;
  MaxPoint[0] = 0;
  MaxPoint[1] = 0;

  for (r=0;r<phdSrc->rows;r++)
  {
    for (c=0;c<phdSrc->cols;c++)
    {
      x = *p++;
      if (x > *Max)
      {
        *Max = x;
        MaxPoint[0] = r;
        MaxPoint[1] = c;
      }
      else if (x < *Min)
      {
        *Min = x;
        MinPoint[0] = r;
        MinPoint[1] = c;
      }
    }
  }
  return(HIPS_OK);
}
