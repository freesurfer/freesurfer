/*
 *       FILE NAME:   mri.c
 *
 *       DESCRIPTION: utilities for MRI  data structure
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        1/8/97
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#define USE_ELECTRIC_FENCE  1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include <errno.h>

#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "filter.h"
#include "box.h"
#include "region.h"
#include "nr.h"
#include "mritransform.h"
#include "utils.h"
#include "matrix.h"

extern int errno;

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define DEBUG_POINT(x,y,z)  (((x==8&&y==9) || (x==9&&y==8)) &&((z)==15))

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static long mris_alloced = 0 ;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRImatch(MRI *mri1, MRI *mri2)
{
  return(
         (mri1->width == mri2->width) &&
         (mri1->height == mri2->height) &&
         (mri1->depth == mri2->depth) &&
         (mri1->type == mri2->type)
         ) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIscalarMul(MRI *mri_src, MRI *mri_dst, float scalar)
{
  int     width, height, depth, x, y, z, frame ;
  BUFTYPE *psrc, *pdst ;
  float   *pfsrc, *pfdst, dval ;
  short   *pssrc, *psdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_UCHAR:
          psrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            dval = *psrc++ * scalar ;
            if (dval < 0)
              dval = 0 ;
            if (dval > 255)
              dval = 255 ;
            *pdst++ = dval ;
          }
          break ;
        case MRI_FLOAT:
          pfsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
          pfdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pfdst++ = *pfsrc++ * scalar ;
          break ;
        case MRI_SHORT:
          pssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          psdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *psdst++ = (short)nint((float)*pssrc++ * scalar) ;
          break ;
        default:
          ErrorReturn(NULL, 
                      (ERROR_UNSUPPORTED, 
                       "MRIscalarMul: unsupported type %d", mri_src->type)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIvalRange(MRI *mri, float *pmin, float *pmax)
{
  int      width, height, depth, x, y, z ;
  float    fmin, fmax, *pf, val ;
  BUFTYPE  *pb ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  switch (mri->type)
  {
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pf = &MRIFvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = *pf++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_INT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRIIvox(mri, x, y, z) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRISvox(mri, x, y, z) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pb = &MRIvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)*pb++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "MRIvalRange: unsupported type %d",
                 mri->type)) ;
  }

  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIvalRangeRegion(MRI *mri, float *pmin, float *pmax, MRI_REGION *region)
{
  int      width, height, depth, x, y, z, x0, y0, z0 ;
  float    fmin, fmax, *pf, val ;
  BUFTYPE  *pb ;

  width = region->x + region->dx ;
  if (width > mri->width)
    width = mri->width ;
  height = region->y + region->dy ;
  if (height > mri->height)
    height = mri->height ;
  depth = region->z + region->dz ;
  if (depth > mri->depth)
    depth = mri->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  switch (mri->type)
  {
  case MRI_FLOAT:
    for (z = z0 ; z < depth ; z++)
    {
      for (y = y0 ; y < height ; y++)
      {
        pf = &MRIFvox(mri, x0, y, z) ;
        for (x = x0 ; x < width ; x++)
        {
          val = *pf++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (z = z0 ; z < depth ; z++)
    {
      for (y = y0 ; y < height ; y++)
      {
        pb = &MRIvox(mri, x0, y, z) ;
        for (x = x0 ; x < width ; x++)
        {
          val = (float)*pb++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "MRIvalRangeRegion: unsupported type %d",
                 mri->type)) ;
  }

  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
MRIclipRegion(MRI *mri, MRI_REGION *reg_src, MRI_REGION *reg_clip)
{
  int  x2, y2, z2 ;

  x2 = MIN(mri->width-1, reg_src->x + reg_src->dx - 1) ;
  y2 = MIN(mri->height-1, reg_src->y + reg_src->dy - 1) ;
  z2 = MIN(mri->depth-1, reg_src->z + reg_src->dz - 1) ;
  reg_clip->x = MAX(0, reg_src->x) ;
  reg_clip->y = MAX(0, reg_src->y) ;
  reg_clip->z = MAX(0, reg_src->z) ;
  reg_clip->dx = x2 - reg_clip->x + 1 ;
  reg_clip->dy = y2 - reg_clip->y + 1 ;
  reg_clip->dz = z2 - reg_clip->z + 1 ;
  return(reg_clip) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIvalScale(MRI *mri_src, MRI *mri_dst, float flo, float fhi) 
{
  int      width, height, depth, x, y, z ;
  float    fmin, fmax, *pf_src, *pf_dst, val, scale ;
  short    *ps_src, *ps_dst ;
  BUFTYPE  *pb_src, *pb_dst ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIvalRange(mri_src, &fmin, &fmax) ;
  scale = (fhi - flo) / (fmax - fmin) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  switch (mri_src->type)
  {
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pf_src = &MRIFvox(mri_src, 0, y, z) ;
        pf_dst = &MRIFvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = *pf_src++ ;
          *pf_dst++ = (val - fmin) * scale + flo ;
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        ps_src = &MRISvox(mri_src, 0, y, z) ;
        ps_dst = &MRISvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)(*ps_src++) ;
          *ps_dst++ = (short)nint((val - fmin) * scale + flo) ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pb_src = &MRIvox(mri_src, 0, y, z) ;
        pb_dst = &MRIvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)*pb_src++ ;
          *pb_dst++ = (BUFTYPE)nint((val - fmin) * scale + flo) ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, "MRIvalScale: unsupported type %d",
                 mri_src->type)) ;
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIconfThresh(MRI *mri_src, MRI *mri_probs, MRI *mri_classes, MRI *mri_dst, 
                     float thresh, int min_target, int max_target)
{
  int      x, y, z, width, height, depth, class ;
  float    *pprobs, prob ;
  BUFTYPE  *pclasses, *pdst, *psrc, src ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_classes, NULL) ;

  width = mri_classes->width ;
  height = mri_classes->height ;
  depth = mri_classes->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pprobs = &MRIFvox(mri_probs, 0, y, z) ;
      pclasses = &MRIvox(mri_classes, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        src = *psrc++ ;
        prob = *pprobs++ ;
        class = (int)*pclasses++ ;
        if (prob >= thresh && ((class >= min_target) && (class <= max_target)))
          *pdst++ = src ;
        else if ((class >= min_target) && (class <= max_target))
          *pdst++ = 25 ;
        else
          *pdst++ = 0 ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIboundingBoxNbhd(MRI *mri, int thresh, int wsize,MRI_REGION *box)
{
  int      width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi, xk, yk, zk,
           whalf, in_brain ;
  BUFTYPE  *psrc ;
  float    *pfsrc ;
  short    *pssrc ;

  whalf = (wsize-1)/2 ;
  box->dx = width = mri->width ;
  box->dy = height = mri->height ;
  box->dz = depth = mri->depth ;

  x1 = y1 = z1 = 0 ;
  box->x = width-1 ;
  box->y = height-1 ;
  box->z = depth-1 ;
  switch (mri->type)
  {
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        psrc = &MRIvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*psrc++ > thresh)
          {
            in_brain = 1 ;
            for (zk = -whalf ; in_brain && zk <= whalf ; zk++)
            {
              zi = mri->zi[z+zk] ;
              for (yk = -whalf ; in_brain && yk <= whalf ; yk++)
              {
                yi = mri->yi[y+yk] ;
                for (xk = -whalf ; in_brain && xk <= whalf ; xk++)
                {
                  xi = mri->xi[x+xk] ;
                  if (MRIvox(mri, xi, yi, zi) < thresh)
                    in_brain = 0 ;
                }
              }
            }
            if (in_brain)
            {
              if (x < box->x)
                box->x = x ;
              if (y < box->y)
                box->y = y ;
              if (z < box->z)
                box->z = z ;
              if (x > x1)
                x1 = x ;
              if (y > y1)
                y1 = y ;
              if (z > z1)
                z1 = z ;
            }
          }
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pssrc = &MRISvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pssrc++ > thresh)
          {
            in_brain = 1 ;
            for (zk = -whalf ; in_brain && zk <= whalf ; zk++)
            {
              zi = mri->zi[z+zk] ;
              for (yk = -whalf ; in_brain && yk <= whalf ; yk++)
              {
                yi = mri->yi[y+yk] ;
                for (xk = -whalf ; in_brain && xk <= whalf ; xk++)
                {
                  xi = mri->xi[x+xk] ;
                  if (MRISvox(mri, xi, yi, zi) < thresh)
                    in_brain = 0 ;
                }
              }
            }
            if (in_brain)
            {
              if (x < box->x)
                box->x = x ;
              if (y < box->y)
                box->y = y ;
              if (z < box->z)
                box->z = z ;
              if (x > x1)
                x1 = x ;
              if (y > y1)
                y1 = y ;
              if (z > z1)
                z1 = z ;
            }
          }
        }
      }
    }
    break ;
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pfsrc = &MRIFvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pfsrc++ > thresh)
          {
            in_brain = 1 ;
            for (zk = -whalf ; in_brain && zk <= whalf ; zk++)
            {
              zi = mri->zi[z+zk] ;
              for (yk = -whalf ; in_brain && yk <= whalf ; yk++)
              {
                yi = mri->yi[y+yk] ;
                for (xk = -whalf ; in_brain && xk <= whalf ; xk++)
                {
                  xi = mri->xi[x+xk] ;
                  if (MRIFvox(mri, xi, yi, zi) < thresh)
                    in_brain = 0 ;
                }
              }
            }
            if (in_brain)
            {
            if (x < box->x)
              box->x = x ;
            if (y < box->y)
              box->y = y ;
            if (z < box->z)
              box->z = z ;
            if (x > x1)
              x1 = x ;
            if (y > y1)
              y1 = y ;
            if (z > z1)
              z1 = z ;
            }
          }
        }
      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "MRIboundingBoxNbd: unsupported type %d",
                 mri->type)) ;
    break ;
  }
  box->x -= whalf+1 ; box->y -= whalf+1 ; box->z -= whalf+1 ;
  x1 += whalf+1 ; y1 += whalf+1 ; z1 += whalf+1 ; 
  box->x = MAX(0,box->x) ; box->y = MAX(0,box->y) ; box->z = MAX(0,box->z) ; 
  x1 = MIN(width-1,x1) ; y1 = MIN(height-1, y1) ; z1 = MIN(depth-1, z1) ;
  box->dx = x1 - box->x + 1 ;
  box->dy = y1 - box->y + 1 ;
  box->dz = z1 - box->z + 1 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MIN_DARK 10

int
MRIfindApproximateSkullBoundingBox(MRI *mri, int thresh,MRI_REGION *box)
{
  int      width, height, depth, x, y, z, x1, y1, z1, ndark, max_dark, start ;
  double   means[3] ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  MRIcenterOfMass(mri, means, thresh) ;


  /* search for left edge */
  ndark = max_dark = 0 ; 
  y = nint(means[1]) ; z = nint(means[2]) ;
  for (start = x1 = x = nint(means[0]) ; x >= 0 ; x--)
  {
    if (MRIvox(mri, x, y, z) < thresh)
    {
      if (!ndark)
        start = x ;
      ndark++ ;
    }
    else
    {
      if (ndark > max_dark)
      {
        max_dark = ndark ; x1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    x1 = start ;
  }
  if (max_dark < MIN_DARK)
    x1 = 0 ;
  box->x = x1 ;

  /* search for right edge */
  ndark = max_dark = 0 ; 
  y = nint(means[1]) ; z = nint(means[2]) ;
  for (start = x1 = x = nint(means[0]) ; x < width ; x++)
  {
    if (MRIvox(mri, x, y, z) < thresh)
    {
      if (!ndark)
        start = x ;
      ndark++ ;
    }
    else
    {
      if (ndark >= max_dark)
      {
        max_dark = ndark ; x1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    x1 = start ;
  }
  if (max_dark < MIN_DARK)
    x1 = mri->width-1 ;
  box->dx = x1 - box->x + 1 ;

  /* search for inferior edge */
  ndark = max_dark = 0 ; 
  x = nint(means[0]) ; z = nint(means[2]) ;
  for (start = y1 = y = nint(means[1]) ; y >= 0 ; y--)
  {
    if (MRIvox(mri, x, y, z) < thresh)
    {
      if (!ndark)
        start = y ;
      ndark++ ;
    }
    else
    {
      if (ndark >= max_dark)
      {
        max_dark = ndark ; y1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    y1 = start ;
  }
  if (max_dark < MIN_DARK)
    y1 = 0 ;
  box->y = y1 ;

  /* search for superior edge */
  ndark = max_dark = 0 ; 
  x = nint(means[0]) ; z = nint(means[2]) ;
  for (start = y = y1 = nint(means[1]) ; y < height ; y++)
  {
    if (MRIvox(mri, x, y, z) < thresh)
    {
      if (!ndark)
        start = y ;
      ndark++ ;
    }
    else
    {
      if (ndark >= max_dark)
      {
        max_dark = ndark ; y1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    y1 = start ;
  }
  if (max_dark < MIN_DARK)
    y1 = mri->height-1 ;
  box->dy = y1 - box->y + 1 ;

  /* search for posterior edge */
  ndark = max_dark = 0 ; 
  x = nint(means[0]) ; y = nint(means[1]) ;
  for (z1 = start = z = nint(means[2]) ; z >= 0 ; z--)
  {
    if (MRIvox(mri, x, y, z) < thresh)
    {
      if (!ndark)
        start = z ;
      ndark++ ;
    }
    else
    {
      if (ndark >= max_dark)
      {
        max_dark = ndark ; z1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    z1 = start ;
  }
  if (max_dark < MIN_DARK)
    z1 = 0 ;
  box->z = z1 ;

  /* search for anterior edge */
  ndark = max_dark = 0 ; 
  x = nint(means[0]) ; y = nint(means[1]) ;
  for (start = z = nint(means[2]) ; z < depth ; z++)
  {
    if (MRIvox(mri, x, y, z) < thresh)
    {
      if (!ndark)
        start = z ;
      ndark++ ;
    }
    else
    {
      if (ndark >= max_dark)
      {
        max_dark = ndark ; z1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    z1 = start ;
  }
  if (max_dark < MIN_DARK)
    z1 = mri->depth-1 ;
  box->dz = z1 - box->z + 1 ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIboundingBox(MRI *mri, int thresh, MRI_REGION *box)
{
  int      width, height, depth, x, y, z, x1, y1, z1 ;
  BUFTYPE  *psrc ;
  float    *pfsrc ;

  box->dx = width = mri->width ;
  box->dy = height = mri->height ;
  box->dz = depth = mri->depth ;

  x1 = y1 = z1 = 0 ;
  box->x = width-1 ;
  box->y = height-1 ;
  box->z = depth-1 ;
  switch (mri->type)
  {
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        psrc = &MRIvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*psrc++ > thresh)
          {
            if (x < box->x)
              box->x = x ;
            if (y < box->y)
              box->y = y ;
            if (z < box->z)
              box->z = z ;
            if (x > x1)
              x1 = x ;
            if (y > y1)
              y1 = y ;
            if (z > z1)
              z1 = z ;
          }
        }
      }
    }
    break ;
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pfsrc = &MRIFvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pfsrc++ > thresh)
          {
            if (x < box->x)
              box->x = x ;
            if (y < box->y)
              box->y = y ;
            if (z < box->z)
              box->z = z ;
            if (x > x1)
              x1 = x ;
            if (y > y1)
              y1 = y ;
            if (z > z1)
              z1 = z ;
          }
        }
      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, "MRIboundingBox: unsupported type %d",
                 mri->type)) ;
    break ;
  }
  box->dx = x1 - box->x + 1 ;
  box->dy = y1 - box->y + 1 ;
  box->dz = z1 - box->z + 1 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIcheckSize(MRI *mri_src, MRI *mri_check, int width, int height, int depth)
{
  if (!mri_check)
    return(0) ;

  if (!width)
    width = mri_src->width ;
  if (!height)
    height = mri_src->height ;
  if (!depth)
    depth = mri_src->depth ;
  
  if (width != mri_check->width ||
      height != mri_check->height ||
      depth != mri_check->depth)
    return(0) ;

  return(1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRItransformRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *src_region,
                           MRI_REGION *dst_region)
{
  Real  xw, yw, zw, xt, yt, zt, xv, yv, zv ;

  if (mri_src->slice_direction != mri_dst->slice_direction)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRItransformRegion(%s): slice directions must match",
                 mri_src->fname)) ;
    
  if (!mri_src->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRItransformRegion(%s): no transform loaded",
                 mri_src->fname)) ;
  if (!mri_dst->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRItransformRegion(%s): no transform loaded",
                 mri_dst->fname)) ;
/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri_src->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIregionToTalairachRegion: unsupported slice direction %d",
                 mri_src->slice_direction)) ;
  }

  xv = (Real)src_region->x ;
  yv = (Real)src_region->y ;
  zv = (Real)src_region->z ;
  MRIvoxelToWorld(mri_src, xv, yv, zv, &xw, &yw, &zw) ;
  transform_point(mri_src->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  transform_point(mri_dst->inverse_linear_transform, xt, yt, zt, &xw,&yw,&zw);
  MRIworldToVoxel(mri_dst, xw, yw, zw, &xv, &yv, &zv) ;
  dst_region->x = nint(xv) ;
  dst_region->y = nint(yv) ;
  dst_region->z = nint(zv) ;

  xv = (Real)(src_region->x + src_region->dx - 1) ;
  yv = (Real)(src_region->y + src_region->dy - 1) ;
  zv = (Real)(src_region->z + src_region->dz - 1) ;
  MRIvoxelToWorld(mri_src, xv, yv, zv, &xw, &yw, &zw) ;
  transform_point(mri_src->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  transform_point(mri_dst->inverse_linear_transform, xt, yt, zt, &xw,&yw,&zw);
  MRIworldToVoxel(mri_dst, xw, yw, zw, &xv, &yv, &zv) ;
  dst_region->dx = nint(xv - (Real)dst_region->x) + 1 ;
  dst_region->dy = nint(yv - (Real)dst_region->y) + 1 ;
  dst_region->dz = nint(zv - (Real)dst_region->z) + 1 ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIvoxelToVoxel(MRI *mri_src, MRI *mri_dst, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt)
{
  Real  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri_src->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToVoxel(%s): no transform loaded", mri_src->fname));
#endif

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri_src->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToVoxel: unsupported slice direction %d",
                 mri_src->slice_direction)) ;
  }

  if (!mri_src->linear_transform || !mri_dst->inverse_linear_transform)
  {
    /* 
       if either doesn't have a transform defined, assume they are in
       the same coordinate system.
       */
    *pxt = xv ; *pyt = yv ; *pzt = zv ;
  }
  else
  {
    MRIvoxelToWorld(mri_src, xv, yv, zv, &xw, &yw, &zw) ;
    if (mri_src->linear_transform)
      transform_point(mri_src->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
    else
    {
      xt = xw ; yt = yw ; zt = zw ;
    }
    if (mri_dst->inverse_linear_transform)
      transform_point(mri_dst->inverse_linear_transform, xt,yt,zt,&xw,&yw,&zw);
    else
    {
      xw = xt ; yw = yt ; zw = zt ;
    }
    MRIworldToVoxel(mri_dst, xw, yw, zw, pxt, pyt, pzt) ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIvoxelToTalairachVoxel(MRI *mri, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt)
{
  Real  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif
/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 mri->slice_direction)) ;
  }

  MRIvoxelToWorld(mri, xv, yv, zv, &xw, &yw, &zw) ;
  if (mri->linear_transform)
    transform_point(mri->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  { xt = xw ; yt = yw ; zt = zw ; }
  MRIworldToVoxel(mri, xt, yt, zt, pxt, pyt, pzt) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIvoxelToTalairach(MRI *mri, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt)
{
  Real  xw, yw, zw ;

#if 0
  if (!mri->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 mri->slice_direction)) ;
  }

  MRIvoxelToWorld(mri, xv, yv, zv, &xw, &yw, &zw) ;
  if (mri->linear_transform)
    transform_point(mri->linear_transform, xw, yw, zw, pxt, pyt, pzt) ;
  else
  { *pxt = xw ; *pyt = yw ; *pzt = zw ; }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRItalairachToVoxel(MRI *mri, Real xt, Real yt, Real zt,
                               Real *pxv, Real *pyv, Real *pzv)
{
  Real  xw, yw, zw ;

#if 0
  if (!mri->inverse_linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRItalairachToVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 mri->slice_direction)) ;
  }

  if (mri->inverse_linear_transform)
    transform_point(mri->inverse_linear_transform, xt, yt, zt, &xw, &yw,&zw);
  else
  { xw = xt ; yw = yt ; zw = zt ; }
  MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRItalairachVoxelToVoxel(MRI *mri, Real xtv, Real ytv, Real ztv,
                               Real *pxv, Real *pyv, Real *pzv)
{
  Real  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri->inverse_linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRItalairachVoxelToVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 mri->slice_direction)) ;
  }

  MRIvoxelToWorld(mri, xtv, ytv, ztv, &xt, &yt, &zt) ;
  if (mri->inverse_linear_transform)
    transform_point(mri->inverse_linear_transform, xt, yt, zt, &xw, &yw,&zw);
  else
  { xw = xt ; yw = yt ; zw = zt ; }
  MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRItalairachVoxelToWorld(MRI *mri, Real xtv, Real ytv, Real ztv,
                               Real *pxw, Real *pyw, Real *pzw)
{
  Real  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri->inverse_linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRItalairachVoxelToVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri->slice_direction)
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, 
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 mri->slice_direction)) ;
  }

  MRIvoxelToWorld(mri, xtv, ytv, ztv, &xt, &yt, &zt) ;
  if (mri->inverse_linear_transform)
    transform_point(mri->inverse_linear_transform, xt, yt, zt, &xw, &yw,&zw);
  else
  { xw = xt ; yw = yt ; zw = zt ; }
  *pxw = xw ; *pyw = yw ; *pzw = zw ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIvoxelToWorld(MRI *mri, Real xv, Real yv, Real zv, 
                Real *pxw, Real *pyw, Real *pzw)
{
  int   ip, jp, kp ;

  switch (mri->slice_direction)
  {
  default:
  case MRI_UNDEFINED:
    /*      ras for MRIvox(mri, i, j, k)    */
    ip = xv - (mri->width-1) / 2;
    jp = yv - (mri->height-1) / 2;
    kp = zv - (mri->depth-1) / 2;

    *pxw = 
      mri->x_r * mri->xsize * ip + 
      mri->y_r * mri->ysize * jp + 
      mri->z_r * mri->zsize * kp + mri->c_r;
    *pyw = 
      mri->x_a * mri->xsize * ip + 
      mri->y_a * mri->ysize * jp + 
      mri->z_a * mri->zsize * kp + mri->c_a;
    *pzw = 
      mri->x_s * mri->xsize * ip + 
      mri->y_s * mri->ysize * jp + 
      mri->z_s * mri->zsize * kp + mri->c_s;
    break ;

  case MRI_CORONAL:
    /* z coordinate system is inverted relative to Talairach space */
    
    /* transform to new origin */
#if 0
    *pxw = (Real)mri->xend - xv * mri->xsize ;
    *pyw = zv * mri->zsize + (Real)mri->zstart ;
    *pzw = -(yv * mri->ysize + (Real)mri->ystart) ;
#else
    trans_SetBounds ( mri->xstart, mri->xend, mri->ystart, mri->yend, 
                      mri->zstart, mri->zend );
    trans_SetResolution ( mri->xsize, mri->ysize, mri->zsize );
    trans_VoxelToRAS(xv, yv, zv, pxw, pyw, pzw) ;
#endif
    /*  fprintf(stderr, "zw = (%d - %2.0f + 1) + %2.1f = %2.1f\n",
        mri->depth, zv, mri->zstart, *pzw) ;*/
    break ;
#if 0
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToWorld: unsupported slice direction %d", 
                 mri->slice_direction)) ;
    break ;
#endif
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIworldToTalairachVoxel(MRI *mri, Real xw, Real yw, Real zw,
                         Real *pxv, Real *pyv, Real *pzv)
{
  Real  xt, yt, zt ;

  transform_point(mri->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  MRIworldToVoxel(mri, xt, yt, zt, pxv, pyv, pzv) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int   MRIworldToVoxelIndex(MRI *mri, Real xw, Real yw, Real zw,
                int *pxv, int *pyv, int *pzv)
{
  switch (mri->slice_direction)
  {
  case MRI_CORONAL:
#if 0
    *pxv = ((Real)mri->xend - xw) / mri->xsize ;
    *pzv = (yw - (Real)mri->zstart) / mri->zsize ;
    *pyv = (-zw - (Real)mri->ystart) / mri->ysize ;
#else
    trans_SetBounds ( mri->xstart, mri->xend, mri->ystart, mri->yend, 
                      mri->zstart, mri->zend );
    trans_SetResolution ( mri->xsize, mri->ysize, mri->zsize );
    trans_RASToVoxelIndex(xw, yw, zw, pxv, pyv, pzv) ;
#endif
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED,
                 "MRIworldToVoxel: unsupported slice direction %d", 
                 mri->slice_direction)) ;
    break ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIworldToVoxel(MRI *mri, Real xw, Real yw, Real zw, 
                Real *pxv, Real *pyv, Real *pzv)
{
  switch (mri->slice_direction)
  {
  case MRI_UNDEFINED:
#if 1
    {
      MATRIX *m_R, *m_R_inv ;
      VECTOR *v_w, *v_v ;
      m_R = MatrixAlloc(3,3,MATRIX_REAL) ;
      *MATRIX_RELT(m_R,1,1) = mri->x_r*mri->xsize ;
      *MATRIX_RELT(m_R,1,2) = mri->y_r*mri->ysize ;
      *MATRIX_RELT(m_R,1,3) = mri->z_r*mri->zsize ;

      *MATRIX_RELT(m_R,2,1) = mri->x_a*mri->xsize ;
      *MATRIX_RELT(m_R,2,2) = mri->y_a*mri->ysize ;
      *MATRIX_RELT(m_R,2,3) = mri->z_a*mri->zsize ;

      *MATRIX_RELT(m_R,3,1) = mri->x_s*mri->xsize ;
      *MATRIX_RELT(m_R,3,2) = mri->y_s*mri->ysize ;
      *MATRIX_RELT(m_R,3,3) = mri->z_s*mri->zsize ;
      
      m_R_inv = MatrixInverse(m_R, NULL) ;
      if (!m_R_inv)
      {
        MatrixPrint(stderr, m_R) ;
        ErrorExit(ERROR_BADPARM, "MRIRASToVoxel: noninvertible xform") ;
      }
      MatrixFree(&m_R) ;
      v_w = VectorAlloc(3, MATRIX_REAL) ;
      V3_LOAD(v_w, xw, yw, zw) ;
      v_v = MatrixMultiply(m_R_inv, v_w, NULL) ;

      *pxv = V3_X(v_v) + (mri->width-1)/2 ;
      *pyv = V3_Y(v_v) + (mri->height-1)/2 ;
      *pzv = V3_Z(v_v) + (mri->depth-1)/2 ;
      VectorFree(&v_v) ; VectorFree(&v_w) ; MatrixFree(&m_R_inv) ;
    }
#else
    *pxv = 
      (mri->width-1)/2 +
      ((mri->x_r * (xw-mri->c_r)) +
       (mri->y_r * (yw-mri->c_a)) +
       (mri->z_r * (zw-mri->c_s)))
      / mri->xsize ;

    *pyv = 
      (mri->height-1)/2 +
      ((mri->x_a * (xw-mri->c_r)) +
       (mri->y_a * (yw-mri->c_a)) +
       (mri->z_a * (zw-mri->c_s)))
      / mri->ysize ;
    *pzv = 
      (mri->depth-1)/2 +
      ((mri->x_s * (xw-mri->c_r)) +
       (mri->y_s * (yw-mri->c_a)) +
       (mri->z_s * (zw-mri->c_s)))
      / mri->zsize ;
#endif
    break;
  case MRI_CORONAL:
#if 0
    *pxv = ((Real)mri->xend - xw) / mri->xsize ;
    *pzv = (yw - (Real)mri->zstart) / mri->zsize ;
    *pyv = (-zw - (Real)mri->ystart) / mri->ysize ;
#else
    trans_SetBounds ( mri->xstart, mri->xend, mri->ystart, mri->yend, 
                      mri->zstart, mri->zend );
    trans_SetResolution ( mri->xsize, mri->ysize, mri->zsize );
    trans_RASToVoxel(xw, yw, zw, pxv, pyv, pzv) ;
#endif
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED,
                 "MRIworldToVoxel: unsupported slice direction %d", 
                 mri->slice_direction)) ;
    break ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIinitHeader(MRI *mri)
{
  mri->ptype = 2 ;

  /* most of these are in mm */
  mri->imnr0 = 1 ;
  mri->imnr1 = mri->depth ;
  mri->fov = mri->width ;
  mri->thick = 1 ;
  mri->xsize = 1 ;
  mri->ysize = 1 ;
  mri->zsize = 1 ;
  mri->ps = 1 ;
  mri->xstart = -mri->width/2 ;
  mri->xend = mri->width/2 ;
  mri->ystart = -mri->height/2 ;
  mri->yend = mri->height/2 ;
  mri->zstart = -mri->depth/2 ;
  mri->zend = mri->depth/2 ;
  mri->x_r = mri->x_a = mri->x_s = 0.0;
  mri->y_r = mri->y_a = mri->y_s = 0.0;
  mri->z_r = mri->z_a = mri->z_s = 0.0;
  mri->c_r = mri->c_a = mri->c_s = 0.0;
  mri->ras_good_flag = 0;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          change the direction of slices
------------------------------------------------------*/
MRI *
MRIextract(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                  int dx, int dy, int dz)
{
  return(MRIextractInto(mri_src, mri_dst, x0, y0, z0, dx, dy, dz, 0, 0, 0)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Extract a cubic region of an MR image and return it to the caller
------------------------------------------------------*/
MRI *
MRIextractRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *region)
{
  return(MRIextractInto(mri_src, mri_dst, region->x, region->y, region->z, 
                        region->dx, region->dy, region->dz, 0, 0, 0)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Extract a cubic region of an MR image and return it to the caller
------------------------------------------------------*/
MRI *
MRIextractIntoRegion(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                            MRI_REGION *region)
{
  return(MRIextractInto(mri_src, mri_dst, x0, y0, z0, region->dx, region->dy, 
                        region->dz, region->x, region->y, region->z)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Extract a cubic region of an MR image and return it to the caller
------------------------------------------------------*/
MRI *
MRIextractInto(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                  int dx, int dy, int dz, int x1, int y1, int z1)
{
  int  width, height, depth, ys, zs, yd, zd, bytes, frame, xsize,ysize,zsize,
       dst_alloced = 0 ;

  width = mri_src->width ;
  depth = mri_src->depth ;
  height = mri_src->height ;

  if (z0 >= depth || y0 >= height || x0 >= width)
    ErrorReturn(NULL,
                (ERROR_BADPARM, 
                 "MRIextractInto: bad src location (%d, %d, %d)", x0,y0,z0));
  if (x0 < 0)
    x0 = 0 ;
  if (y0 < 0)
    y0 = 0 ;
  if (z0 < 0)
    z0 = 0 ;
  if (x0+dx > width)
    dx = (width - x0) ;
  if (y0+dy > height)
    dy = (height - y0) ;
  if (z0+dz > depth)
    dz = (depth - z0) ;
  if (x1 < 0)
    x1 = 0 ;
  if (y1 < 0)
    y1 = 0 ;
  if (z1 < 0)
    z1 = 0 ;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(dx, dy, dz, mri_src->type, mri_src->nframes) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->imnr0 = z0 + mri_src->imnr0 - z1 ;
    mri_dst->imnr1 = mri_dst->imnr0 + dz - 1 ;
    dst_alloced = 1 ;
  }

  if (mri_src->type != mri_dst->type)
  {
    MRIfree(&mri_dst) ;
    ErrorReturn(NULL, 
                (ERROR_BADPARM, 
                 "MRIextractInto: src and dst types must match"));
  }

  if (x1+dx > mri_dst->width)
    dx = (mri_dst->width - x1) ;
  if (y1+dy > mri_dst->height)
    dy = (mri_dst->height - y1) ;
  if (z1+dz > mri_dst->depth)
    dz = (mri_dst->depth - z1) ;

  xsize = mri_src->xsize ;
  ysize = mri_src->ysize ;
  zsize = mri_src->zsize ;

  if (dst_alloced)
  {
    mri_dst->xstart += x0*xsize ;
    mri_dst->xend = mri_dst->xstart + dx*xsize ;
    mri_dst->ystart += y0*ysize ;
    mri_dst->yend = mri_dst->ystart + dy*ysize ;
    mri_dst->zstart += z0*zsize ;
    mri_dst->zend = mri_dst->zstart + dz*zsize  ;
  }

  bytes = dx ;
  switch (mri_src->type)
  {
  case MRI_FLOAT:
    bytes *= sizeof(float) ;
    break ;
  case MRI_LONG:
    bytes *= sizeof(long) ;
    break ;
  case MRI_INT:
    bytes *= sizeof(int) ;
    break ;
  default:
    break ;
  }

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (zd = z1, zs = z0 ; zs < z0+dz ; zs++, zd++)
    {
      for (yd = y1, ys = y0 ; ys < y0+dy ; ys++, yd++)
      {
        switch (mri_src->type)
        {
        case MRI_UCHAR:
          memcpy(&MRIseq_vox(mri_dst, x1, yd, zd,frame), 
                 &MRIseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        case MRI_FLOAT:
          memcpy(&MRIFseq_vox(mri_dst, x1, yd, zd,frame), 
                 &MRIFseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        }
      }
    }
  }
#if 0
  mri_dst->imnr0 = z0 + mri_src->imnr0 - z1;
  mri_dst->imnr1 = mri_dst->imnr0+mri_dst->depth-1 ;
#endif
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          change the direction of slices
------------------------------------------------------*/
MRI *
MRIreslice(MRI *mri_src, MRI *mri_dst, int slice_direction)
{
  int     width, height, depth, x1, x2, x3 ;
  BUFTYPE *psrc, val, *pdst ;

  if (slice_direction == mri_src->slice_direction)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
    return(mri_dst) ;
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if ((mri_src->slice_direction == MRI_SAGITTAL && 
       slice_direction == MRI_CORONAL) || 
      (mri_src->slice_direction == MRI_CORONAL && 
       slice_direction == MRI_SAGITTAL))
  {
/*
   coronal images are back to front of the head, thus the depth axis
   points from the nose to the back of the head, with x from neck to 
   crown, and y from ear to ear.
*/
    /* x1 --> x3
       x2 --> x2
       x3 --> x1
       */
    if (!mri_dst)
    {
      mri_dst = MRIalloc(depth, height, width, mri_src->type) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else if (depth != mri_dst->width || height != mri_dst->height ||
             width != mri_dst->depth)
    {
      ErrorReturn(NULL, 
                  (ERROR_BADPARM,
                   "MRIreslice: invalid destination size (%d, %d, %d)",
                   mri_dst->width, mri_dst->height, mri_dst->depth)) ;
    }
    
    for (x3 = 0 ; x3 < depth ; x3++)
    {
      for (x2 = 0 ; x2 < height ; x2++)
      {
        psrc = &MRIvox(mri_src, 0, x2, x3) ;
        for (x1 = 0 ; x1 < width ; x1++)
        {
          /* swap so in place transformations are possible */
          mri_dst->slices[x1][x2][x3] = *psrc++ ;
        }
      }
    }
  }
  else
    if ((mri_src->slice_direction == MRI_HORIZONTAL && 
         slice_direction == MRI_CORONAL) || 
        (mri_src->slice_direction == MRI_CORONAL && 
         slice_direction == MRI_HORIZONTAL))
    {
/*
   horizontal images are top to bottom of the head, thus the depth axis
   points from the top of the head to the neck, with x from ear to ear
   and y from nose to back of head.
*/
    /* x3 --> x2
       x2 --> x3
       x1 --> x1
       */
    if (!mri_dst)
    {
      mri_dst = MRIalloc(width, depth, height, mri_src->type) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else if (depth != mri_dst->height || height != mri_dst->depth ||
             width != mri_dst->width)
      ErrorReturn(NULL, 
                  (ERROR_BADPARM,
                   "MRIreslice: invalid destination size (%d, %d, %d)",
                   mri_dst->width, mri_dst->height, mri_dst->depth)) ;
    
    for (x3 = 0 ; x3 < depth ; x3++)
    {
      for (x2 = 0 ; x2 < height ; x2++)
      {
        psrc = &MRIvox(mri_src, 0, x2, x3) ;
        pdst = &MRIvox(mri_dst, 0, x3, x2) ;
        for (x1 = 0 ; x1 < width ; x1++)
        {
          /* swap so in place transformations are possible */
          *pdst++ = *psrc++ ;
        }
      }
    }
  }
  else
    if ((mri_src->slice_direction == MRI_SAGITTAL && 
         slice_direction == MRI_HORIZONTAL))
    {
/*
   horizontal images are top to bottom of the head, thus the depth axis
   points from the top of the head to the neck, with x from ear to ear
   and y from nose to back of head.
*/
    /* x3 --> x2
       x1 --> x3
       x2 --> x1
       */
    if (!mri_dst)
    {
      mri_dst = MRIalloc(width, depth, height, mri_src->type) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else if (depth != mri_dst->height || height != mri_dst->depth ||
             width != mri_dst->width)
      ErrorReturn(NULL, 
                  (ERROR_BADPARM,
                   "MRIreslice: invalid destination size (%d, %d, %d)",
                   mri_dst->width, mri_dst->height, mri_dst->depth)) ;
    
    for (x3 = 0 ; x3 < depth ; x3++)
    {
      for (x2 = 0 ; x2 < height ; x2++)
      {
        psrc = &MRIvox(mri_src, 0, x2, x3) ;
        for (x1 = 0 ; x1 < width ; x1++)
        {
          /* swap so in place transformations are possible */
          mri_dst->slices[x2][x1][x3] = *psrc++ ;
        }
      }
    }
  }
  else
    if (mri_src->slice_direction == MRI_HORIZONTAL && 
         slice_direction == MRI_SAGITTAL)
    {
/*
   horizontal images are top to bottom of the head, thus the depth axis
   points from the top of the head to the neck, with x from ear to ear
   and y from nose to back of head.
*/
    /* x2 --> x3
       x3 --> x1
       x1 --> x2
       */
    if (!mri_dst)
    {
      mri_dst = MRIalloc(width, depth, height, mri_src->type) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else if (depth != mri_dst->height || height != mri_dst->depth ||
             width != mri_dst->width)
      ErrorReturn(NULL, 
                  (ERROR_BADPARM,
                   "MRIreslice: invalid destination size (%d, %d, %d)",
                   mri_dst->width, mri_dst->height, mri_dst->depth)) ;
    
    for (x3 = 0 ; x3 < depth ; x3++)
    {
      for (x2 = 0 ; x2 < height ; x2++)
      {
        psrc = &MRIvox(mri_src, 0, x2, x3) ;
        for (x1 = 0 ; x1 < width ; x1++)
        {
          /* swap so in place transformations are possible */
          mri_dst->slices[x1][x3][x2] = *psrc++ ;
        }
      }
    }
  }
  else
    switch (mri_src->slice_direction)
  {
  default:
    MRIfree(&mri_dst) ;
    ErrorReturn(NULL,
                (ERROR_BADPARM, 
                 "MRIreslice: mri_src unknown slice direction %d", 
                 mri_src->slice_direction)) ;
    break ;
  case MRI_CORONAL:
/*
   coronal images are back to front of the head, thus the depth axis
   points from the nose to the back of the head, with x from neck to 
   crown, and y from ear to ear.
*/
    switch (slice_direction)
    {
    case MRI_SAGITTAL:
      /* x1 --> x3
         x2 --> x2
         x3 --> x1
         */
      if (!mri_dst)
      {
        mri_dst = MRIalloc(depth, height, width, mri_src->type) ;
        MRIcopyHeader(mri_src, mri_dst) ;
      }
      else if (depth != mri_dst->width || height != mri_dst->height ||
               width != mri_dst->depth)
        ErrorReturn(NULL, 
                    (ERROR_BADPARM,
                     "MRIreslice: invalid destination size (%d, %d, %d)",
                     mri_dst->width, mri_dst->height, mri_dst->depth)) ;

      for (x3 = 0 ; x3 < depth ; x3++)
      {
        for (x2 = 0 ; x2 < height ; x2++)
        {
          psrc = &MRIvox(mri_src, 0, x2, x3) ;
          for (x1 = 0 ; x1 < width ; x1++)
          {
            /* swap so in place transformations are possible */
            val = *psrc++ ;
#if 0
            mri_dst->slices[x3][x2][x1] = mri_src->slices[x1][x2][x3] ;
#endif
            mri_dst->slices[x1][x2][x3] = val ;
          }
        }
      }
      break ;
    case MRI_HORIZONTAL:
      break ;
    }
    break ;
  case MRI_SAGITTAL:
/*
   sagittal images are slices in the plane of the nose, with depth going
   from ear to ear.
*/
    break ;
  }

  mri_dst->slice_direction = slice_direction ;
  mri_dst->ras_good_flag = 0;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Set an MRI intensity values to 0
------------------------------------------------------*/
int
MRIclear(MRI *mri)
{
  int   width, depth, height, bytes, y, z, frame, nframes ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  nframes = mri->nframes ;
  bytes = width ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    bytes *= sizeof(unsigned char) ;
    break ;
  case MRI_BITMAP:
    bytes /= 8 ;
    break ;
  case MRI_FLOAT:
    bytes *= sizeof(float) ;
    break ;
  case MRI_LONG:
    bytes *= sizeof(long) ;
    break ;
  case MRI_INT:
    bytes *= sizeof(int) ;
    break ;
  case MRI_SHORT:
    bytes *= sizeof(short) ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIclear: unsupported input type %d", mri->type)) ;
    break ;
  }

  for (frame = 0 ; frame < nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        memset(mri->slices[z+frame*depth][y], 0, bytes) ;
    }
  }
  
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          find the principle components of a (binary) MRI. The
          eigenvectors are the columns of the matrix mEvectors, the
          eigenvalues are returned in the array evalues and the means
          in means (these last two must be three elements long.)
------------------------------------------------------*/
int
MRIcenterOfMass(MRI *mri,double *means, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, val ;
  long    npoints ;
  double  mx, my, mz, weight ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIcenterOfMass: unsupported input type %d", mri->type)) ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  mx = my = mz = weight = 0.0f ; npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          weight += val ;
          mx += (float)x*val ;
          my += (float)y*val ;
          mz += (float)z*val ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;


  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          find the principle components of a (binary) MRI. The
          eigenvectors are the columns of the matrix mEvectors, the
          eigenvalues are returned in the array evalues and the means
          in means (these last two must be three elements long.)
------------------------------------------------------*/
int
MRIprincipleComponents(MRI *mri, MATRIX *mEvectors, float *evalues, 
                       double *means, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, val ;
  long    npoints ;
  MATRIX  *mCov, *mX, *mXT, *mTmp ;
  double  mx, my, mz, weight ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIprincipleComponents: unsupported input type %d", 
                 mri->type)) ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  mx = my = mz = weight = 0.0f ; npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          weight += val ;
          mx += (float)x*val ;
          my += (float)y*val ;
          mz += (float)z*val ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;     /* zero-mean coordinate vector */
  mXT = NULL ;                              /* transpose of above */
  mTmp = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* tmp matrix for covariance */
  mCov = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* covariance matrix */

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          mX->rptr[1][1] = (x - (int)mx)*val ;
          mX->rptr[2][1] = (y - (int)my)*val ;
          mX->rptr[3][1] = (z - (int)mz)*val ;
          mXT = MatrixTranspose(mX, mXT) ;
          mTmp = MatrixMultiply(mX, mXT, mTmp) ;
          mCov = MatrixAdd(mTmp, mCov, mCov) ;
        }
      }
    }
  }

  if (weight > 0)
    MatrixScalarMul(mCov, 1.0f/weight, mCov) ;

  MatrixEigenSystem(mCov, evalues, mEvectors) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          find the principle components of a (binary) MRI. The
          eigenvectors are the columns of the matrix mEvectors, the
          eigenvalues are returned in the array evalues and the means
          in means (these last two must be three elements long.)
------------------------------------------------------*/
int
MRIbinaryPrincipleComponents(MRI *mri, MATRIX *mEvectors, float *evalues, 
                       double *means, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, val ;
  long    npoints ;
  MATRIX  *mCov, *mX, *mXT, *mTmp ;
  double  mx, my, mz, weight ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIprincipleComponents: unsupported input type %d", 
                 mri->type)) ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  mx = my = mz = weight = 0.0f ; npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          weight++ ;
          mx += (float)x ;
          my += (float)y ;
          mz += (float)z ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;     /* zero-mean coordinate vector */
  mXT = NULL ;                              /* transpose of above */
  mTmp = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* tmp matrix for covariance */
  mCov = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* covariance matrix */

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          mX->rptr[1][1] = (x - (int)mx) ;
          mX->rptr[2][1] = (y - (int)my) ;
          mX->rptr[3][1] = (z - (int)mz) ;
          mXT = MatrixTranspose(mX, mXT) ;
          mTmp = MatrixMultiply(mX, mXT, mTmp) ;
          mCov = MatrixAdd(mTmp, mCov, mCov) ;
        }
      }
    }
  }

  if (weight > 0)
    MatrixScalarMul(mCov, 1.0f/weight, mCov) ;

  MatrixEigenSystem(mCov, evalues, mEvectors) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           threshold an MRI.
------------------------------------------------------*/
MRI *
MRIthresholdRangeInto(MRI *mri_src,MRI *mri_dst,BUFTYPE low_val,BUFTYPE hi_val)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, val ;
  float   *pfsrc, *pfdst, fval ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  switch (mri_src->type)
  {
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        psrc = &MRIvox(mri_src, 0, y, z) ;
        pdst = &MRIvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++, pdst++)
        {
          val = *psrc++ ;
          if (val >=  low_val && val <= hi_val)
            *pdst = val ;
        }
      }
    }
    break ;
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pfsrc = &MRIFvox(mri_src, 0, y, z) ;
        pfdst = &MRIFvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++, pdst++)
        {
          fval = *pfsrc++ ;
          if (fval >=  low_val && fval <= hi_val)
            *pfdst = fval ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, 
                 "MRIthresholdRangeInto: unsupported type %d",
                 mri_src->type)) ;
    break ;
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           threshold an MRI.
------------------------------------------------------*/
MRI *
MRIthreshold(MRI *mri_src, MRI *mri_dst, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val < threshold)
          val = 0 ;
        *pdst++ = val ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           threshold an MRI.
------------------------------------------------------*/
MRI  *
MRIinvertContrast(MRI *mri_src, MRI *mri_dst, float threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
          val = 255 - val ;
        *pdst++ = val ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           threshold an MRI.
------------------------------------------------------*/
MRI *
MRIbinarize(MRI *mri_src, MRI *mri_dst, BUFTYPE threshold, BUFTYPE low_val,
            BUFTYPE hi_val)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val < threshold)
          val = low_val ;
        else
          val = hi_val ;
        *pdst++ = val ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIsubtract(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
        *pdst++ = *p1++ - *p2++ ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIabsdiff(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst, v1, v2 ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        v1 = *p1++ ;
        v2 = *p2++ ;
        if (v1 > v2)
          *pdst++ = v1 - v2 ;
        else
          *pdst++ = v2 - v1 ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIabs(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, b ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = mri_src->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        b = *psrc++ ;
        *pdst++ = abs(b) ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIadd(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
        *pdst++ = *p1++ + *p2++ ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIaverage(MRI *mri_src, int dof, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, src, dst ;
  BUFTYPE *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (!MRIcheckSize(mri_src, mri_dst,0,0,0))
    ErrorReturn(NULL, 
                (ERROR_BADPARM,"MRISaverage: incompatible volume dimensions"));
  if ((mri_src->type != MRI_UCHAR) || (mri_dst->type != MRI_UCHAR))
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED,
                 "MRISaverage: unsupported voxel format %d",mri_src->type));
    
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        src = (int)*pdst * dof ;
        dst = (int)*psrc++ ;
        *pdst++ = (BUFTYPE)((src + dst) / (dof + 1)) ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRImultiply(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
        *pdst++ = *p1++ * *p2++ ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIscaleAndMultiply(MRI *mri1, float scale, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;
  float   out_val ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        out_val = *p1++ * (*p2++/scale) ;
        if (out_val > 255)
          out_val = 255 ;
        else if (out_val < 0)
          out_val = 0 ;
        *pdst++ = (BUFTYPE)nint(out_val) ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIdivide(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        if  (!*p2)
        {
          p2++ ;
          *pdst = 255 ;
        }
        else
          *pdst++ = *p1++ / *p2++ ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Copy one MRI into another (including header info)
------------------------------------------------------*/
MRI *
MRIclone(MRI *mri_src, MRI *mri_dst)
{
  if (!mri_dst)
    mri_dst = 
      MRIallocSequence(mri_src->width, mri_src->height,mri_src->depth, 
                       mri_src->type, mri_src->nframes);

  MRIcopyHeader(mri_src, mri_dst) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Copy one MRI into another (including header info)
------------------------------------------------------*/
MRI *
MRIcloneRoi(MRI *mri_src, MRI *mri_dst)
{
  int  w, h, d ;

  w = mri_src->width - mri_src->roi.x ;
  h = mri_src->height - mri_src->roi.y ;
  d = mri_src->depth - mri_src->roi.z ;
  mri_dst = MRIallocSequence(w, h, d, MRI_FLOAT, mri_src->nframes) ;
  MRIcopyHeader(mri_src, mri_dst) ;
  mri_dst->xstart = mri_src->xstart + mri_src->roi.x * mri_src->xsize ;
  mri_dst->ystart = mri_src->ystart + mri_src->roi.y * mri_src->ysize ;
  mri_dst->zstart = mri_src->zstart + mri_src->roi.z * mri_src->zsize ;
  mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
  mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
  mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Copy one MRI into another (including header info)
------------------------------------------------------*/
MRI *
MRIcopy(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, bytes, x, y, z, frame, val ;
  float   *fdst, *fsrc ;
  BUFTYPE *csrc, *cdst ;
  int     dest_ptype, *isrc;
  short   *ssrc, *sdst ;

  if (mri_src == mri_dst)
    return(mri_dst) ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    if(mri_src->slices)
      mri_dst = MRIallocSequence(width, height, depth, mri_src->type,
                               mri_src->nframes) ;
    else
      mri_dst = MRIallocHeader(width, height, depth, mri_src->type);
  }
  dest_ptype = mri_dst->ptype;
  MRIcopyHeader(mri_src, mri_dst) ;
  mri_dst->ptype = dest_ptype;

  if(!mri_src->slices)
    return(mri_dst);

  if (mri_src->type == mri_dst->type)
  {
    bytes = width ;
    switch (mri_src->type)
    {
    case MRI_UCHAR:
      bytes *= sizeof(BUFTYPE) ;
      break ;
    case MRI_SHORT:
      bytes *= sizeof(short);
      break;
    case MRI_FLOAT:
      bytes *= sizeof(float) ;
      break ;
    case MRI_INT:
      bytes *= sizeof(int) ;
      break ;
    case MRI_LONG:
      bytes *= sizeof(long) ;
      break ;
    }

    for (frame = 0 ; frame < mri_src->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          memcpy(mri_dst->slices[z+frame*depth][y], 
                 mri_src->slices[z+frame*depth][y], bytes) ;
        }
      }
    }
  }
  else
  {
    switch (mri_src->type)
    {
    case MRI_FLOAT:
      switch (mri_dst->type)
      {
      case MRI_SHORT: /* float --> short */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              sdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
              fsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                val = nint(*fsrc++) ;
                *sdst++ = (short)val ;
              }
            }
          }
        }
        break ;
      case MRI_UCHAR:  /* float --> unsigned char */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              cdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
              fsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                val = nint(*fsrc++) ;
                if (val > 255)
                  val = 255 ;
                *cdst++ = (BUFTYPE)val ;
              }
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM, 
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
        }
      break ;
    case MRI_UCHAR:
      switch (mri_dst->type)
      {
      case MRI_FLOAT:   /* unsigned char --> float */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              fdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
              csrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
                *fdst++ = (float)*csrc++ ;
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM, 
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    case MRI_SHORT:
      switch (mri_dst->type)
      {
      case MRI_FLOAT:   /* short --> float */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              fdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
              ssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                if (z == 113 && y == 143 && x == 161)
                  DiagBreak() ;
                *fdst++ = (float)*ssrc++ ;
              }
            }
          }
        }
        break ;
      case MRI_UCHAR:
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              cdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
              ssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                *cdst++ = (float)*ssrc++ ;
              }
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM, 
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    case MRI_INT:
      switch (mri_dst->type)
      {
      case MRI_FLOAT:   /* unsigned char --> float */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              fdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
              isrc = &MRIIseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
                *fdst++ = (float)*isrc++ ;
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM, 
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    default:
      ErrorReturn(NULL,
                  (ERROR_BADPARM, 
                   "MRIcopy: src type %d & dst type %d unsupported",
                   mri_src->type, mri_dst->type)) ;
      break ;  /* in case someone removes the errorreturn */
    }
  }
  return(mri_dst) ;
}
/* 
   make MAX_INDEX way larger than it has to be. This will give
   some headroom for bad (e.g. poorly registered) images without
   sacrificing too much space.
   */
#define MAX_INDEX    500
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          allocate a lookup table that allows indices which
          are outside the image region.
------------------------------------------------------*/
int
MRIallocIndices(MRI *mri)
{
  int width, height, depth, i ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  mri->xi = (int *)calloc(width+2*MAX_INDEX, sizeof(int)) ;
  if (!mri->xi)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRIallocIndices: could not allocate %d elt index array",
              width+2*MAX_INDEX) ;
  mri->yi = (int *)calloc(height+2*MAX_INDEX, sizeof(int)) ;
  if (!mri->yi)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRIallocIndices: could not allocate %d elt index array",
              height+2*MAX_INDEX) ;
  mri->zi = (int *)calloc(depth+2*MAX_INDEX, sizeof(int)) ;
  if (!mri->zi)
    ErrorExit(ERROR_NO_MEMORY, 
              "MRIallocIndices: could not allocate %d elt index array",
              depth+2*MAX_INDEX) ;

/* 
   indexing into these arrays returns valid pixel indices from  
   -MAX_INDEX to width+MAX_INDEX
   */
  mri->xi += MAX_INDEX ;
  mri->yi += MAX_INDEX ;
  mri->zi += MAX_INDEX ;
  for (i = -MAX_INDEX ; i < width+MAX_INDEX ; i++)
  {
    if (i <= 0)
      mri->xi[i] = 0 ;
    else if (i >= width)
      mri->xi[i] = width-1 ;
    else
      mri->xi[i] = i ;
  }
  for (i = -MAX_INDEX ; i < height+MAX_INDEX ; i++)
  {
    if (i <= 0)
      mri->yi[i] = 0 ;
    else if (i >= height)
      mri->yi[i] = height-1 ;
    else
      mri->yi[i] = i ;
  }
  for (i = -MAX_INDEX ; i < depth+MAX_INDEX ; i++)
  {
    if (i <= 0)
      mri->zi[i] = 0 ;
    else if (i >= depth)
      mri->zi[i] = depth-1 ;
    else
      mri->zi[i] = i ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           allocate an MRI data structure as well as space for
           the image data
------------------------------------------------------*/
MRI *
MRIallocSequence(int width, int height, int depth, int type, int nframes)
{
  MRI     *mri ;
  int     slice, row, bpp ;
  BUFTYPE *buf ;

  mris_alloced++ ;

  if ((width <= 0) || (height <= 0) || (depth <= 0))
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIalloc(%d, %d, %d): bad parm",
                 width, height, depth)) ;
#if 1
  mri = MRIallocHeader(width, height, depth, type) ;
  MRIinitHeader(mri) ;
#else
  mri = (MRI *)calloc(1, sizeof(MRI)) ;
  if (!mri)
    ErrorExit(ERROR_NO_MEMORY, "MRIalloc: could not allocate MRI\n") ;

  mri->scale = 1 ;
  mri->height = height ;
  mri->width = width ;
  mri->yinvert = 1 ;
  mri->depth = depth ;
  mri->type = type ;
#endif
  mri->nframes = nframes ;
  MRIallocIndices(mri) ;
  mri->slices = (BUFTYPE ***)calloc(depth*nframes, sizeof(BUFTYPE **)) ;
  if (!mri->slices)
      ErrorExit(ERROR_NO_MEMORY, 
                "MRIalloc: could not allocate %d slices\n", mri->depth) ;

  for (slice = 0 ; slice < depth*nframes ; slice++)
  {
    /* allocate pointer to array of rows */
    mri->slices[slice] = (BUFTYPE **)calloc(mri->height, sizeof(BUFTYPE *)) ;
    if (!mri->slices[slice])
      ErrorExit(ERROR_NO_MEMORY, 
        "MRIalloc(%d, %d, %d): could not allocate %d bytes for %dth slice\n",
        height, width, depth, mri->height*sizeof(BUFTYPE *), slice) ;

#if USE_ELECTRIC_FENCE
    switch (mri->type)
    {
    case MRI_BITMAP:
      bpp = 1 ;
      break ;
    case MRI_UCHAR:
      bpp = 8 ;
      break ;
    case MRI_TENSOR:
    case MRI_FLOAT:
      bpp = sizeof(float) * 8 ;
      break ;
    case MRI_INT:
      bpp = sizeof(int) * 8 ;
      break ;
    case MRI_SHORT:
      bpp = sizeof(short) * 8 ;
      break ;
    case MRI_LONG:
      bpp = sizeof(long) * 8 ;
      break ;
    default:
      ErrorReturn(NULL, 
                  (ERROR_BADPARM,
                   "MRIalloc(%d, %d, %d, %d): unknown type",
                   width, height, depth, mri->type)) ;
      break ;
    }
    bpp /= 8 ;
    buf = (BUFTYPE *)calloc((mri->width*mri->height*bpp), 1) ;
    for (row = 0 ; row < mri->height ; row++)
    {
      mri->slices[slice][row] = buf+(row*mri->width*bpp) ;
    }
#else
    /* allocate each row */
    for (row = 0 ; row < mri->height ; row++)
    {
      switch (mri->type)
      {
      case MRI_BITMAP:
        mri->slices[slice][row] = 
          (BUFTYPE*)calloc(mri->width/8,sizeof(BUFTYPE));
        break ;
      case MRI_UCHAR:
        mri->slices[slice][row] = (BUFTYPE*)calloc(mri->width,sizeof(BUFTYPE));
        break ;
      case MRI_TENSOR:
      case MRI_FLOAT:
        mri->slices[slice][row] = (BUFTYPE *)calloc(mri->width, sizeof(float));
        break ;
      case MRI_INT:
        mri->slices[slice][row] = (BUFTYPE *)calloc(mri->width, sizeof(int)) ;
        break ;
      case MRI_SHORT:
        mri->slices[slice][row] = (BUFTYPE *)calloc(mri->width, sizeof(short));
        break ;
      case MRI_LONG:
        mri->slices[slice][row] = (BUFTYPE *)calloc(mri->width, sizeof(long)) ;
        break ;
      default:
        ErrorReturn(NULL, 
                    (ERROR_BADPARM,
                     "MRIalloc(%d, %d, %d, %d): unknown type",
                     width, height, depth, mri->type)) ;
        break ;
      }

      if (!mri->slices[slice][row])
        ErrorExit(ERROR_NO_MEMORY, 
          "MRIalloc(%d,%d,%d): could not allocate %dth row in %dth slice\n",
         width,height,depth, slice, row) ;
    }
#endif
  }

  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           allocate an MRI data structure but not space for
           the image data
------------------------------------------------------*/
MRI *
MRIallocHeader(int width, int height, int depth, int type)
{
  MRI  *mri ;

  mri = (MRI *)calloc(1, sizeof(MRI)) ;
  if (!mri)
    ErrorExit(ERROR_NO_MEMORY, "MRIalloc: could not allocate MRI\n") ;

  mri->xdir = XDIM;
  mri->ydir = YDIM;
  mri->zdir = ZDIM;
  mri->scale = 1 ;
  mri->roi.dx = mri->width = width ;
  mri->roi.dy = mri->height = height ;
  mri->roi.dz = mri->depth = depth ;
  mri->yinvert = 1 ;
  mri->xsize = mri->ysize = mri->zsize = 1 ;
  mri->type = type ;
  mri->nframes = 1 ;
  mri->xi = mri->yi = mri->zi = NULL ;
  mri->slices = NULL ;
  mri->x_r = mri->x_a = mri->x_s = 0.0;
  mri->y_r = mri->y_a = mri->y_s = 0.0;
  mri->z_r = mri->z_a = mri->z_s = 0.0;
  mri->c_r = mri->c_a = mri->c_s = 0.0;
  mri->ras_good_flag = 0;
  mri->brightness = 1;
  mri->register_mat = NULL;
  mri->subject_name[0] = '\0';
  mri->path_to_t1[0] = '\0';
  mri->fname_format[0] = '\0';
  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           allocate an MRI data structure as well as space for
           the image data
------------------------------------------------------*/
MRI *
MRIalloc(int width, int height, int depth, int type)
{
  return(MRIallocSequence(width, height, depth, type, 1)) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
            Free and MRI data structure and all its attached memory
------------------------------------------------------*/
int
MRIfree(MRI **pmri)
{
  MRI *mri ;
  int slice ;
#if !USE_ELECTRIC_FENCE
  int  row ;
#endif

  mris_alloced-- ;  
  mri = *pmri ;
  if (!mri)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIfree: null pointer\n")) ;

  if (mri->xi)
    free(mri->xi-MAX_INDEX) ;
  if (mri->yi)
    free(mri->yi-MAX_INDEX) ;
  if (mri->zi)
    free(mri->zi-MAX_INDEX) ;

  if (mri->slices)
  {
    for (slice = 0 ; slice < mri->depth*mri->nframes ; slice++)
    {
      if (mri->slices[slice])
      {
#if USE_ELECTRIC_FENCE
        free(mri->slices[slice][0]) ;
#else
        for (row = 0 ; row < mri->height ; row++)
          free(mri->slices[slice][row]) ;
#endif
        free(mri->slices[slice]) ;
      }
    }
    free(mri->slices) ;
  }

  if (mri->free_transform)
    delete_general_transform(&mri->transform) ;

  if(mri->register_mat != NULL)
    MatrixFree(&(mri->register_mat));

  free(mri) ;
  *pmri = NULL ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
            Free and MRI data structure and all its attached memory
------------------------------------------------------*/
int
MRIfreeFrames(MRI *mri, int start_frame)
{
  int slice, row, end_frame ;

  end_frame = mri->nframes-1 ;
  if (!mri)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIfree: null pointer\n")) ;

  if (mri->slices)
  {
    for (slice = start_frame*mri->depth ; 
         slice < mri->depth*end_frame ; slice++)
    {
      if (mri->slices[slice])
      {
        for (row = 0 ; row < mri->height ; row++)
          free(mri->slices[slice][row]) ;
        free(mri->slices[slice]) ;
      }
    }
  }

  mri->nframes -= (end_frame-start_frame+1) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Dump the MRI header to a file
------------------------------------------------------*/
int
MRIdump(MRI *mri, FILE *fp)
{
  fprintf(fp, "%6.6s = %s\n", "fname", mri->fname);
  fprintf(fp, "%6.6s = %d\n", "height", mri->height);
  fprintf(fp, "%6.6s = %d\n", "width", mri->width);
  fprintf(fp, "%6.6s = %d\n", "depth", mri->depth);
  fprintf(fp, "%6.6s = %d\n", "nframes", mri->nframes);
  fprintf(fp, "%6.6s = %d\n", "imnr0", mri->imnr0);
  fprintf(fp, "%6.6s = %d\n", "imnr1", mri->imnr1);
  fprintf(fp, "%6.6s = %d\n", "xnum", mri->width);
  fprintf(fp, "%6.6s = %d\n", "ynum", mri->height);
  fprintf(fp, "%6.6s = %f\n", "fov", mri->fov);
  fprintf(fp, "%6.6s = %f\n", "thick", mri->thick);
  fprintf(fp, "%6.6s = %f\n", "xstart", mri->xstart); /* strtx */
  fprintf(fp, "%6.6s = %f\n", "xend", mri->xend); /* endx */
  fprintf(fp, "%6.6s = %f\n", "ystart", mri->ystart); /* strty */
  fprintf(fp, "%6.6s = %f\n", "yend", mri->yend); /* endy */
  fprintf(fp, "%6.6s = %f\n", "zstart", mri->zstart); /* strtz */
  fprintf(fp, "%6.6s = %f\n", "zend", mri->zend); /* endz */
  fprintf(fp, "%6.6s = %d\n", "type", mri->type);
  fprintf(fp, "%6.6s = %d\n", "sl dir", mri->slice_direction);
  fprintf(fp, "%6.6s = %f\n", "xsize", mri->xsize);
  fprintf(fp, "%6.6s = %f\n", "ysize", mri->ysize);
  fprintf(fp, "%6.6s = %f\n", "zsize", mri->zsize);
  fprintf(fp, "%6.6s = %f %f %f\n", "x ras", mri->x_r, mri->x_a, mri->x_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "y ras", mri->y_r, mri->y_a, mri->y_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "z ras", mri->z_r, mri->z_a, mri->z_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "c ras", mri->c_r, mri->c_a, mri->c_s);
  fprintf(fp, "%s = %d\n", "ras_good_flag", mri->ras_good_flag);
  fprintf(fp, "%6.6s = %d\n", "xdir", mri->xdir);
  fprintf(fp, "%6.6s = %d\n", "ydir", mri->ydir);
  fprintf(fp, "%6.6s = %d\n", "zdir", mri->zdir);
  fprintf(fp, "%s = %d\n", "brightness", mri->brightness);
  fprintf(fp, "%s = %s\n", "subject_name", mri->subject_name);
  fprintf(fp, "%s = %s\n", "path_to_t1", mri->path_to_t1);
  fprintf(fp, "%s = %s\n", "fname_format", mri->fname_format);
  if(mri->register_mat != NULL)
  {
    fprintf(fp, "%s = \n", "register_mat");
    MatrixPrint(fp, mri->register_mat);
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Dump the non-zero elements of an MRI buffer to a file
------------------------------------------------------*/
int
MRIdumpBuffer(MRI *mri, FILE *fp)
{
  int  x, y, z ;

  for (z = 0 ; z < mri->depth ; z++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (x = 0 ; x < mri->width ; x++)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          if (!FZERO(mri->slices[z][y][x]))
            fprintf(fp, "[%d][%d][%d]: %d\n", x,y,z,(int)mri->slices[z][y][x]);
          break ;
        case MRI_FLOAT:
          /*          if (!FZERO(MRIFvox(mri,x,y,z)))*/
            fprintf(fp, "[%d][%d][%d]: %2.3f\n", x,y,z, MRIFvox(mri,x,y,z));
          break ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Find the peak intensity in an MRI image
------------------------------------------------------*/
int 
MRIpeak(MRI *mri, int *px, int *py, int *pz)
{
  int      max_row, max_col, max_slice, row, col, slice, width, height, depth ;
  BUFTYPE  val, max_val, *im ;
  long     lval, lmax_val, *lim ;
  int      ival, imax_val, *iim ;
  float    fval, fmax_val, *fim ;

  max_val = 0 ;
  lmax_val = 0L ;
  imax_val = 0 ;
  fmax_val = 0.0f ;
  max_row = max_col = max_slice = -1 ;   /* to prevent compiler warning */
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  for (slice = 0 ; slice < depth ; slice++)
  {
    for (row = 0 ; row < height ; row++)
    {
      switch (mri->type)
      {
      case MRI_UCHAR:
        im = mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          val = *im++ ;
          if (val > max_val)
          {
            max_val = val ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      case MRI_LONG:
        lim = (long *)mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          lval = *lim++ ;
          if (lval > lmax_val)
          {
            lmax_val = lval ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      case MRI_FLOAT:
        fim = (float *)mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          fval = *fim++ ;
          if (fval > fmax_val)
          {
            fmax_val = fval ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      case MRI_INT:
        iim = (int *)mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          ival = *iim++ ;
          if (ival > imax_val)
          {
            imax_val = ival ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      }
    }
  }

  *px = max_col ;
  *py = max_row ;
  *pz = max_slice ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Copy the header information from one MRI into
          another.
------------------------------------------------------*/
int
MRIcopyHeader(MRI *mri_src, MRI *mri_dst)
{
  mri_dst->dof = mri_src->dof ;
  mri_dst->mean = mri_src->mean ;
  mri_dst->xsize = mri_src->xsize ;
  mri_dst->ysize = mri_src->ysize ;
  mri_dst->zsize = mri_src->zsize ;
  if (mri_src->linear_transform)
  {
    copy_general_transform(&mri_src->transform, &mri_dst->transform) ;
    mri_dst->linear_transform = mri_src->linear_transform ;
    mri_dst->inverse_linear_transform = mri_src->inverse_linear_transform ;
    mri_dst->linear_transform = get_linear_transform_ptr(&mri_dst->transform) ;
    mri_dst->inverse_linear_transform = 
      get_inverse_linear_transform_ptr(&mri_dst->transform) ;
    mri_dst->free_transform = 1 ;
  }
  strcpy(mri_dst->transform_fname, mri_src->transform_fname) ;
  mri_dst->slice_direction = mri_src->slice_direction ;
  if (mri_dst->depth == mri_src->depth)
  {
    mri_dst->imnr0 = mri_src->imnr0 ;
    mri_dst->imnr1 = mri_src->imnr1 ;
  }
  mri_dst->ptype = mri_src->ptype ;
  mri_dst->fov = mri_src->fov ;
  mri_dst->thick = mri_src->thick ;
  mri_dst->ps = mri_src->ps ;
  mri_dst->location = mri_src->location ;
  mri_dst->xstart = mri_src->xstart ;
  mri_dst->xend = mri_src->xend ;
  mri_dst->ystart = mri_src->ystart ;
  mri_dst->yend = mri_src->yend ;
  mri_dst->zstart = mri_src->zstart ;
  mri_dst->zend = mri_src->zend ;
  mri_dst->flip_angle = mri_src->flip_angle ;
  mri_dst->tr = mri_src->tr ;
  mri_dst->te = mri_src->te ;
  mri_dst->ti = mri_src->ti ;
  strcpy(mri_dst->fname, mri_src->fname) ;
  mri_dst->xdir = mri_src->xdir;
  mri_dst->ydir = mri_src->ydir;
  mri_dst->zdir = mri_src->zdir;
  mri_dst->x_r = mri_src->x_r;
  mri_dst->x_a = mri_src->x_a;
  mri_dst->x_s = mri_src->x_s;
  mri_dst->y_r = mri_src->y_r;
  mri_dst->y_a = mri_src->y_a;
  mri_dst->y_s = mri_src->y_s;
  mri_dst->z_r = mri_src->z_r;
  mri_dst->z_a = mri_src->z_a;
  mri_dst->z_s = mri_src->z_s;
  mri_dst->c_r = mri_src->c_r;
  mri_dst->c_a = mri_src->c_a;
  mri_dst->c_s = mri_src->c_s;
  mri_dst->ras_good_flag = mri_src->ras_good_flag;

  mri_dst->brightness = mri_src->brightness;
  if(mri_src->register_mat != NULL)
    MatrixCopy(mri_src->register_mat, mri_dst->register_mat);
  else
    mri_dst->register_mat = NULL;
  strcpy(mri_dst->subject_name, mri_src->subject_name);
  strcpy(mri_dst->path_to_t1, mri_src->path_to_t1);
  strcpy(mri_dst->fname_format, mri_src->fname_format);

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Translate the MRI image mri_src by the vector
          dx, dy, dz  and store the result in mri_dst.
------------------------------------------------------*/
MRI *
MRItranslate(MRI *mri_src, MRI *mri_dst, double dx, double dy, double dz)
{
  int      y1, y2, y3, width, height, depth ;
  BUFTYPE *pdst ;
  Real     x1, x2, x3, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  for (y3 = 0 ; y3 < depth ; y3++)
  {
    x3 = (Real)y3 - dz ;
    if (x3 < 0 || x3 >= depth)
      continue ;

    for (y2 = 0 ; y2 < height ; y2++)
    {
      x2 = (Real)y2 - dy ;
      if (x2 < 0 || x2 >= height)
        continue ;

      pdst = &MRIvox(mri_dst, 0, y2, y3) ;
      for (y1 = 0 ; y1 < width ; y1++, pdst++)
      {
        x1 = (Real)y1 - dx ;
        if (x1 >= 0 && x1 < width)
        {
          MRIsampleVolume(mri_src, x1, x2, x3, &val) ;
          *pdst = (BUFTYPE)nint(val) ;
        }
      }
    }
  }

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Rotate mri_src around the Y axis and return the
           result in mri_dst
------------------------------------------------------*/
MRI *
MRIrotateX(MRI *mri_src, MRI *mri_dst, float x_angle)
{
  int    width, height, depth ;
  MATRIX *m, *mO ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, x_angle, X_ROTATION) ;

  mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
  mO->rptr[1][1] = mri_src->width / 2 ;
  mO->rptr[2][1] = mri_src->height / 2 ;
  mO->rptr[3][1] = mri_src->depth / 2 ;

  /* build rotation matrix */
  mri_dst = MRIrotate(mri_src, NULL, m, mO) ;
  MatrixFree(&m) ;
  MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Rotate mri_src around the Y axis and return the
           result in mri_dst
------------------------------------------------------*/
MRI *
MRIrotateY(MRI *mri_src, MRI *mri_dst, float y_angle)
{
  int    width, height, depth ;
  MATRIX *m, *mO ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  /* origin of coordinate system */
  mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
  mO->rptr[1][1] = mri_src->width / 2 ;
  mO->rptr[2][1] = mri_src->height / 2 ;
  mO->rptr[3][1] = mri_src->depth / 2 ;

  m = MatrixAllocRotation(3, y_angle, Y_ROTATION) ;

  /* build rotation matrix */
  mri_dst = MRIrotate(mri_src, NULL, m, mO) ;
  MatrixFree(&m) ;
  MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Rotate mri_src around the Z axis and return the
           result in mri_dst
------------------------------------------------------*/
MRI *
MRIrotateZ(MRI *mri_src, MRI *mri_dst, float z_angle)
{
  int    width, height, depth ;
  MATRIX *m, *mO ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, z_angle, Z_ROTATION) ;
  mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
  mO->rptr[1][1] = mri_src->width / 2 ;
  mO->rptr[2][1] = mri_src->height / 2 ;
  mO->rptr[3][1] = mri_src->depth / 2 ;

  mri_dst = MRIrotate(mri_src, NULL, m, mO) ;
  MatrixFree(&m) ;
  MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Rotate mri_src around the Y axis and return the
           result in mri_dst
------------------------------------------------------*/
MRI *
MRIrotateX_I(MRI *mri_src, MRI *mri_dst, float x_angle)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, x_angle, X_ROTATION) ;

  /* build rotation matrix */
  mri_dst = MRIrotate_I(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Rotate mri_src around the Y axis and return the
           result in mri_dst
------------------------------------------------------*/
MRI *
MRIrotateY_I(MRI *mri_src, MRI *mri_dst, float y_angle)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, y_angle, Y_ROTATION) ;

  /* build rotation matrix */
  mri_dst = MRIrotate_I(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Rotate mri_src around the Z axis and return the
           result in mri_dst
------------------------------------------------------*/
MRI *
MRIrotateZ_I(MRI *mri_src, MRI *mri_dst, float z_angle)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, z_angle, Z_ROTATION) ;
  mri_dst = MRIrotate_I(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Scale the MRI image mri_src by sx,sy,sz in the
          x, y, and z directions respectively.
------------------------------------------------------*/
MRI *
MRIscale(MRI *mri_src, MRI *mri_dst, float sx, float sy, float sz)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  m = MatrixAlloc(3, 3, MATRIX_REAL) ;

  /* build rotation matrix */
  m->rptr[1][1] = sx ;
  m->rptr[2][2] = sy ;
  m->rptr[3][3] = sz ;
  mri_dst = MRIaffine(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Rotate about the point mO
------------------------------------------------------*/
MRI *
MRIrotate(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO)
{
  int    x1, x2, x3, y1, y2, y3, width, height, depth, y1o, y2o, y3o, freeit ;
  MATRIX *mX, *mY ;   /* original and transformed coordinate systems */
  MATRIX *mRinv ;   /* inverse of R */

  mRinv = MatrixInverse(mR, NULL) ;
  if (!mRinv)
    ErrorReturn(NULL,(ERROR_BADPARM,"MRIrotate: rotation matrix is singular"));

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  if (!mO)
  {
    mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
    mO->rptr[1][1] = mri_src->width / 2 ;
    mO->rptr[2][1] = mri_src->height / 2 ;
    mO->rptr[3][1] = mri_src->depth / 2 ;
    freeit = 1 ;
  }
  else
    freeit = 0 ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* input coordinates */
  mY = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* transformed coordinates */


  y1o = mO->rptr[1][1] ;
  y2o = mO->rptr[2][1] ;
  y3o = mO->rptr[3][1] ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    mY->rptr[3][1] = y3 - y3o ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      mY->rptr[2][1] = y2 - y2o ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        mY->rptr[1][1] = y1 - y1o ;
        MatrixMultiply(mRinv, mY, mX) ;
        MatrixAdd(mX, mO, mX) ;
        
        /* should do bilinear interpolation here */
        x1 = nint(mX->rptr[1][1]) ;
        x2 = nint(mX->rptr[2][1]) ;
        x3 = nint(mX->rptr[3][1]) ;
        if (x1 >= 0 && x1 < width &&
            x2 >= 0 && x2 < height &&
            x3 >= 0 && x3 < depth)
          mri_dst->slices[y3][y2][y1] = mri_src->slices[x3][x2][x1] ;
      }
    }
  }

  MatrixFree(&mX) ;
  MatrixFree(&mRinv) ;
  MatrixFree(&mY) ;
  if (freeit)
    MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Rotate about the point mO using trilinear interpolation
------------------------------------------------------*/
MRI *
MRIrotate_I(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO)
{
  int    y1, y2, y3, width, height, depth, y1o, y2o, y3o, freeit ;
  MATRIX *mX, *mY ;   /* original and transformed coordinate systems */
  MATRIX *mRinv ;   /* inverse of R */
  float  x1, x2, x3 ;
  Real   val ;

  mRinv = MatrixInverse(mR, NULL) ;
  if (!mRinv)
    ErrorReturn(NULL,(ERROR_BADPARM,
                      "MRIrotate_I: rotation matrix is singular"));

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  if (!mO)
  {
    mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
    mO->rptr[1][1] = mri_src->width / 2 ;
    mO->rptr[2][1] = mri_src->height / 2 ;
    mO->rptr[3][1] = mri_src->depth / 2 ;
    freeit = 1 ;
  }
  else
    freeit = 0 ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* input coordinates */
  mY = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* transformed coordinates */


  y1o = mO->rptr[1][1] ;
  y2o = mO->rptr[2][1] ;
  y3o = mO->rptr[3][1] ;
  if (Gdiag == 99)
    MatrixPrint(stdout, mRinv) ;
  if (Gdiag == 99)
    MatrixPrint(stdout, mO) ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    mY->rptr[3][1] = y3 - y3o ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      mY->rptr[2][1] = y2 - y2o ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        mY->rptr[1][1] = y1 - y1o ;
        MatrixMultiply(mRinv, mY, mX) ;
        MatrixAdd(mX, mO, mX) ;
        
        /* do trilinear interpolation here */
        x1 = mX->rptr[1][1] ;
        x2 = mX->rptr[2][1] ;
        x3 = mX->rptr[3][1] ;

        MRIsampleVolume(mri_src, x1, x2, x3, &val) ;
        mri_dst->slices[y3][y2][y1] = (BUFTYPE)nint(val) ;
      }
    }
  }

  MatrixFree(&mX) ;
  MatrixFree(&mRinv) ;
  MatrixFree(&mY) ;
  if (freeit)
    MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Perform an affine coordinate transformation x' = Ax + B on
          the MRI image mri_src into mri_dst
------------------------------------------------------*/
MRI *
MRIaffine(MRI *mri_src, MRI *mri_dst, MATRIX *mA, MATRIX *mB)
{
  int    x1, x2, x3, y1, y2, y3, width, height, depth ;
  MATRIX *mX, *mY ;   /* original and transformed coordinate systems */

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* input coordinates */
  mY = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* transformed coordinates */

  for (x3 = 0 ; x3 < depth ; x3++)
  {
    mX->rptr[3][1] = x3 ;
    for (x2 = 0 ; x2 < height ; x2++)
    {
      mX->rptr[2][1] = x2 ;
      for (x1 = 0 ; x1 < width ; x1++)
      {
        mX->rptr[1][1] = x1 ;
        if (mA)
          MatrixMultiply(mA, mX, mY) ;
        else
          MatrixCopy(mX, mY) ;
        if (mB)
          MatrixAdd(mY, mB, mY) ;
        y1 = nint(mY->rptr[1][1]) ;
        y2 = nint(mY->rptr[2][1]) ;
        y3 = nint(mY->rptr[3][1]) ;
        if (y1 >= 0 && y1 < width &&
            y2 >= 0 && y2 < height &&
            y3 >= 0 && y3 < depth)
          mri_dst->slices[y3][y2][y1] = mri_src->slices[x3][x2][x1] ;
      }
    }
  }


  MatrixFree(&mX) ;
  MatrixFree(&mY) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Convert a slice of an MRI data structure into a HIPS image.
------------------------------------------------------*/
IMAGE *
MRItoImageView(MRI *mri, IMAGE *I, int slice, int view, int frame)
{
  int      width, height, depth, x, y, yp, w, h, d, *idst,
           xm, ym, zm, format ;
  float    *fdst ;
  BUFTYPE  *bdst ;

  d = w = h = xm = ym = zm = 0 ;  /* silly compiler warnings */
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (view == mri->slice_direction)
  {
    w = width ;
    h = height ;
    d = depth ;
  }
  else if (mri->slice_direction != MRI_CORONAL)
  {
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRItoImageView(%d, %d): unsupported view/slice direction %d",
                 slice, view, mri->slice_direction)) ;
  }
  else switch (view)
  {
  default:
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRItoImageView(%d, %d): unsupported view/slice direction %d",
                 slice, view, mri->slice_direction)) ;
    break ;
  case MRI_CORONAL:
    w = width ;
    h = height ;
    d = depth ;
    break ;
  case MRI_SAGITTAL:
    w = depth ;
    h = height ;
    d = width ;
    break ;
  case MRI_HORIZONTAL:
    w = width ;
    h = depth ;
    d = height ;
    break ;
  }

  if (slice < 0 || slice >= d)
    ErrorReturn(NULL, (ERROR_BADPARM, "MRItoImageView: bad slice %d\n",slice));

  format = (mri->type == MRI_UCHAR) ? PFBYTE :
    mri->type == MRI_INT   ? PFINT :
      mri->type == MRI_FLOAT ? PFFLOAT : 
        mri->type == MRI_SHORT ? PFFLOAT : PFBYTE ; 

  if (I && ((I->rows != h) || (I->cols != w) || (I->pixel_format != format)))
    I = NULL ;  /* must allocate a new one */

  if (!I)
    I = ImageAlloc(h, w, format, 1) ;

  switch (mri->type)
  {
  case MRI_SHORT:
    idst = (int *)IMAGEIpix(I, 0, 0) ;
    break ;
  case MRI_INT:
    idst = (int *)IMAGEIpix(I, 0, 0) ;
    break ;
  case MRI_FLOAT:
    fdst = IMAGEFpix(I, 0, 0) ;
    break ;
  case MRI_UCHAR:
    bdst = IMAGEpix(I, 0, 0) ;
    break ;
  }

  for (y = 0 ; y < h ; y++)
  {
    yp = h - (y+1) ;   /* hips coordinate system is inverted */
    for (x = 0 ; x < w ; x++)
    {
      if (view == mri->slice_direction)
      {
        xm = x ;
        ym = y ;
        zm = slice ;
      }
      else switch (view)   /* calculate coordinates in MR structure */
      {
      case MRI_CORONAL:
        xm = x ;
        ym = y ;
        zm = slice ;
        break ;
      case MRI_SAGITTAL:
        xm = slice ; ;
        ym = y ;
        zm = x ;
        break ;
      case MRI_HORIZONTAL:
        xm = x ;
        ym = slice ;
        zm = y ;
        break ;
      }

      switch (mri->type)
      {
#if 1
      case MRI_INT:
        *IMAGEIpix(I, x, yp) = MRIIseq_vox(mri, xm, ym, zm, frame) ;
        break ;
      case MRI_FLOAT:
        *IMAGEFpix(I, x, yp) = MRIFseq_vox(mri, xm, ym, zm, frame) ;
        break ;
      case MRI_UCHAR:
        *IMAGEpix(I, x, yp) = MRIseq_vox(mri, xm, ym, zm, frame) ;
        break ;
      case MRI_SHORT:
        *IMAGEFpix(I, x, yp) = MRISseq_vox(mri, xm, ym, zm, frame) ;
        break ;
#else
      case MRI_INT:
        *idst++ = MRIIvox(mri, xm, ym, zm) ;
        break ;
      case MRI_FLOAT:
        *fdst++ = MRIFvox(mri, xm, ym, zm) ;
        break ;
      case MRI_UCHAR:
        *bdst++ = MRIvox(mri, xm, ym, zm) ;
        break ;
#endif
      default:
      case MRI_LONG:
        ImageFree(&I) ;
        ErrorReturn(NULL, 
                    (ERROR_UNSUPPORTED, "MRItoImageView: unsupported type %d",
                     mri->type)) ;
        break ;
      }
    }
  }

  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Convert a slice of an MRI data structure into a HIPS image.
------------------------------------------------------*/
IMAGE *
MRItoImage(MRI *mri, IMAGE *I, int slice)
{
  int           width, height, y, yp ;
  
  width = mri->width ;
  height = mri->height ;
  if (slice < 0 || slice >= mri->depth)
    ErrorReturn(NULL, (ERROR_BADPARM, "MRItoImage: bad slice %d\n", slice)) ;

  if (!I)
  {
    I = ImageAlloc(height, width, 
                   mri->type == MRI_UCHAR ? PFBYTE :
                   mri->type == MRI_INT   ? PFINT :
                   mri->type == MRI_FLOAT ? PFFLOAT : PFBYTE, 1) ; 
  }

  for (y = 0 ; y < height ; y++)
  {
    yp = height - (y+1) ;

#if 0
    if (mri->type == MRI_UCHAR)
      image_to_buffer(I->image, mri, slice) ;
    else 
#endif
      switch (mri->type)
    {
    case MRI_INT:
      memcpy(IMAGEIpix(I, 0, yp), mri->slices[slice][y], width*sizeof(int)) ;
      break ;
    case MRI_FLOAT:
      memcpy(IMAGEFpix(I, 0, yp), mri->slices[slice][y], width*sizeof(float)) ;
      break ;
    case MRI_UCHAR:
      memcpy(IMAGEpix(I, 0, yp), mri->slices[slice][y], 
              width*sizeof(unsigned char)) ;
      break ;
    default:
    case MRI_LONG:
      ImageFree(&I) ;
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRItoImage: unsupported type %d",
                         mri->type)) ;
      break ;
    }
  }

  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Build an MRI from all values s.t. min_val <= val <= max_val
------------------------------------------------------*/
MRI *
MRIextractValues(MRI *mri_src, MRI *mri_dst, float min_val, float max_val)
{
  BUFTYPE  *psrc, *pdst, val ;
  float    *pfsrc, *pfdst, fval ;
  int      frame, x, y, z, width, height, depth ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_UCHAR:
          psrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            val = *psrc++ ;
            if ((val < min_val) || (val > max_val))
              val = 0 ;
            *pdst++ = val ;
          }
          break ;
        case MRI_FLOAT:
          pfsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
          pfdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            fval = *pfsrc++ ;
            if ((fval < min_val) || (fval > max_val))
              fval = 0.0f ;
            *pfdst++ = fval ;
          }
          break ;
        default:
          ErrorReturn(NULL, 
                      (ERROR_UNSUPPORTED, 
                       "MRIextractValues: unsupported type %d",mri_src->type));
        }
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIupsample2(MRI *mri_src, MRI *mri_dst)
{
  int     width, depth, height, x, y, z ;
  BUFTYPE *pdst ;
  
  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "MRIupsample2: source must be UCHAR"));

  width = 2*mri_src->width ;
  height = 2*mri_src->height ;
  depth = 2*mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
        *pdst++ = MRIvox(mri_src, x/2, y/2, z/2) ;
    }
  }
  
  mri_dst->imnr0 = mri_src->imnr0 ;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1 ;

  mri_dst->xsize = mri_src->xsize/2 ;
  mri_dst->ysize = mri_src->ysize/2 ;
  mri_dst->zsize = mri_src->zsize/2 ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
MRI *
MRIdownsample2LabeledVolume(MRI *mri_src, MRI *mri_dst)
{
  int     width, depth, height, x, y, z, x1, y1, z1, counts[256], label, 
          max_count, out_label ;
  BUFTYPE *psrc ;
  
  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRIdownsample2LabeledVolume: source must be UCHAR"));

  width = mri_src->width/2 ;
  height = mri_src->height/2 ;
  depth = mri_src->depth/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  MRIclear(mri_dst) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        memset(counts, 0, sizeof(counts)) ;
        if (x == 96 && y == 66 && z == 56)
          DiagBreak() ;
        for (z1 = 2*z ; z1 <= 2*z+1 ; z1++)
        {
          for (y1 = 2*y ; y1 <= 2*y+1 ; y1++)
          {
            psrc = &MRIvox(mri_src, 2*x, y1, z1) ;
            for (x1 = 2*x ; x1 <= 2*x+1 ; x1++)
            {
              label = *psrc++ ;
              counts[label]++ ;
            }
          }
        }
        for (out_label = label = 0, max_count = counts[0] ; label <= 255 ; label++)
        {
          if (counts[label] > max_count)
          {
            out_label = label ;
            max_count = counts[label] ;
          }
        }
        MRIvox(mri_dst, x, y, z) = out_label ;
      }
    }
  }

  mri_dst->imnr0 = mri_src->imnr0 ;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1 ;
  mri_dst->xsize = mri_src->xsize*2 ;
  mri_dst->ysize = mri_src->ysize*2 ;
  mri_dst->zsize = mri_src->zsize*2 ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIdownsample2(MRI *mri_src, MRI *mri_dst)
{
  int     width, depth, height, x, y, z, x1, y1, z1 ;
  BUFTYPE *psrc ;
  float   val ;
  
  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "MRIdownsample2: source must be UCHAR"));

  width = mri_src->width/2 ;
  height = mri_src->height/2 ;
  depth = mri_src->depth/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  MRIclear(mri_dst) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        for (val = 0.0f, z1 = 2*z ; z1 <= 2*z+1 ; z1++)
        {
          for (y1 = 2*y ; y1 <= 2*y+1 ; y1++)
          {
            psrc = &MRIvox(mri_src, 2*x, y1, z1) ;
            for (x1 = 2*x ; x1 <= 2*x+1 ; x1++)
              val += *psrc++ ;
          }
        }
        MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(val/8.0f) ;
      }
    }
  }

  mri_dst->imnr0 = mri_src->imnr0 ;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1 ;
  mri_dst->xsize = mri_src->xsize*2 ;
  mri_dst->ysize = mri_src->ysize*2 ;
  mri_dst->zsize = mri_src->zsize*2 ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Iteratively set all voxels in mri_dst that neighbor
          a voxel that has already been filled (starting with the seed),
          and for which the corresponding voxel in mri_src is
          below threshold.
------------------------------------------------------*/
MRI *
MRIfill(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y, 
        int seed_z, int threshold, int fill_val, int max_count)
{
  int     width, height, depth, x, y, z, nfilled, xmin, xmax, ymin, ymax,
          zmin, zmax, on, x0, x1, y0, y1, z0, z1 ;
  BUFTYPE *psrc, *pdst, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (seed_x < 0)
    seed_x = 0 ;
  else if (seed_x >= width)
    seed_x = width-1 ;
  if (seed_y < 0)
    seed_y = 0 ;
  else if (seed_y >= height)
    seed_y = height-1 ;
  if (seed_z < 0)
    seed_z = 0 ;
  else if (seed_z >= depth)
    seed_z = depth-1 ;

  xmin = MAX(seed_x-1,0) ; xmax = MIN(seed_x+1, width-1) ;
  ymin = MAX(seed_y-1,0) ; ymax = MIN(seed_y+1, height-1) ;
  zmin = MAX(seed_z-1,0) ; zmax = MIN(seed_z+1, depth-1) ;

  /* replace all occurrences of fill_val with fill_val-1 */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++, pdst++)
      {
        val = *pdst ;
        if (val == fill_val)
          *pdst = val-1 ;
      }
    }
  }

  MRIvox(mri_dst, seed_x, seed_y, seed_z) = fill_val ;
  do
  {
    z0 = zmin ; z1 = zmax ; y0 = ymin ; y1 = ymax ; x0 = xmin ; x1 = xmax ;
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        psrc = &MRIvox(mri_src, xmin, y, z) ;
        pdst = &MRIvox(mri_dst, xmin, y, z) ;
        for (x = xmin ; x <= xmax ; x++, psrc++, pdst++)
        {
          val = *psrc ;
          if ((val > threshold) || (*pdst == fill_val))
            continue ;

          on = 0 ;
          if (x > 0)
            on = (MRIvox(mri_dst, x-1, y, z) == fill_val) ;
          if (!on && (x < width-1))
            on = (MRIvox(mri_dst, x+1, y, z) == fill_val) ;
          if (!on && (y > 0))
            on = (MRIvox(mri_dst, x, y-1, z) == fill_val) ;
          if (!on && (y < height-1))
            on = (MRIvox(mri_dst, x, y+1, z) == fill_val) ;
          if (!on && (z > 0))
            on = (MRIvox(mri_dst, x, y, z-1) == fill_val) ;
          if (!on && (z < depth-1))
            on = (MRIvox(mri_dst, x, y, z+1) == fill_val) ;
          if (on)
          {
            if (x <= x0)
              x0 = MAX(x-1,0) ;
            if (x >= x1)
              x1 = MIN(x+1,width-1) ;
            if (y <= y0)
              y0 = MAX(y-1,0) ;
            if (y >= y1)
              y1 = MIN(y+1,height-1) ;
            if (z <= z0)
              z0 = MAX(z-1,0) ;
            if (z >= z1)
              z1 = MIN(z+1,depth-1) ;
            nfilled++ ;
            *pdst = fill_val ;
          }
        }
      }
    }
    zmin = z0 ; zmax = z1 ; ymin = y0 ; ymax = y1 ; xmin = x0 ; xmax = x1 ;
/*    fprintf(stderr, "# filled = %d\n", nfilled) ;*/
    if ((max_count > 0) && (nfilled > max_count))
      break ;
  } while (nfilled > 0) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIfillFG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y, 
        int seed_z, int threshold, int fill_val, int *npix)
{
  int     width, height, depth, x, y, z, nfilled, xmin, xmax, ymin, ymax,
          zmin, zmax, on, x0, x1, y0, y1, z0, z1, total_pix ;
  BUFTYPE *psrc, *pdst, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (seed_x < 0)
    seed_x = 0 ;
  else if (seed_x >= width)
    seed_x = width-1 ;
  if (seed_y < 0)
    seed_y = 0 ;
  else if (seed_y >= height)
    seed_y = width-1 ;
  if (seed_z < 0)
    seed_z = 0 ;
  else if (seed_z >= depth)
    seed_z = width-1 ;

  xmin = MAX(seed_x-1,0) ; xmax = MIN(seed_x+1, width-1) ;
  ymin = MAX(seed_y-1,0) ; ymax = MIN(seed_y+1, height-1) ;
  zmin = MAX(seed_z-1,0) ; zmax = MIN(seed_z+1, depth-1) ;

  /* replace all occurrences of fill_val with fill_val-1 */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++, pdst++)
      {
        val = *pdst ;
        if (val == fill_val)
          *pdst = val-1 ;
      }
    }
  }

  MRIvox(mri_dst, seed_x, seed_y, seed_z) = fill_val ;
  total_pix = 1 ;  /* include the seed point */
  do
  {
    z0 = zmin ; z1 = zmax ; y0 = ymin ; y1 = ymax ; x0 = xmin ; x1 = xmax ;
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        psrc = &MRIvox(mri_src, xmin, y, z) ;
        pdst = &MRIvox(mri_dst, xmin, y, z) ;
        for (x = xmin ; x <= xmax ; x++, psrc++, pdst++)
        {
          val = *psrc ;
          if ((val < threshold) || (*pdst == fill_val))
            continue ;

          on = 0 ;
          if (x > 0)
            on = (MRIvox(mri_dst, x-1, y, z) == fill_val) ;
          if (!on && (x < width-1))
            on = (MRIvox(mri_dst, x+1, y, z) == fill_val) ;
          if (!on && (y > 0))
            on = (MRIvox(mri_dst, x, y-1, z) == fill_val) ;
          if (!on && (y < height-1))
            on = (MRIvox(mri_dst, x, y+1, z) == fill_val) ;
          if (!on && (z > 0))
            on = (MRIvox(mri_dst, x, y, z-1) == fill_val) ;
          if (!on && (z < depth-1))
            on = (MRIvox(mri_dst, x, y, z+1) == fill_val) ;
          if (on)
          {
            if (x <= x0)
              x0 = MAX(x-1,0) ;
            if (x >= x1)
              x1 = MIN(x+1,width-1) ;
            if (y <= y0)
              y0 = MAX(y-1,0) ;
            if (y >= y1)
              y1 = MIN(y+1,height-1) ;
            if (z <= z0)
              z0 = MAX(z-1,0) ;
            if (z >= z1)
              z1 = MIN(z+1,depth-1) ;
            nfilled++ ;
            *pdst = fill_val ;
          }
        }
      }
    }
    zmin = z0 ; zmax = z1 ; ymin = y0 ; ymax = y1 ; xmin = x0 ; xmax = x1 ;
    total_pix += nfilled ;
/*    fprintf(stderr, "# filled = %d\n", nfilled) ;*/
  } while (nfilled > 0) ;

  if (npix)
    *npix = total_pix ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIfillBG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y, 
        int seed_z, int threshold, int fill_val, int *npix)
{
  int     width, height, depth, x, y, z, nfilled, xmin, xmax, ymin, ymax,
          zmin, zmax, on, x0, x1, y0, y1, z0, z1, total_pix ;
  BUFTYPE *psrc, *pdst, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (seed_x < 0)
    seed_x = 0 ;
  else if (seed_x >= width)
    seed_x = width-1 ;
  if (seed_y < 0)
    seed_y = 0 ;
  else if (seed_y >= height)
    seed_y = width-1 ;
  if (seed_z < 0)
    seed_z = 0 ;
  else if (seed_z >= depth)
    seed_z = width-1 ;

  xmin = MAX(seed_x-1,0) ; xmax = MIN(seed_x+1, width-1) ;
  ymin = MAX(seed_y-1,0) ; ymax = MIN(seed_y+1, height-1) ;
  zmin = MAX(seed_z-1,0) ; zmax = MIN(seed_z+1, depth-1) ;

  /* replace all occurrences of fill_val with fill_val-1 */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++, pdst++)
      {
        val = *pdst ;
        if (val == fill_val)
          *pdst = val-1 ;
      }
    }
  }

  MRIvox(mri_dst, seed_x, seed_y, seed_z) = fill_val ;
  total_pix = 0 ;
  do
  {
    z0 = zmin ; z1 = zmax ; y0 = ymin ; y1 = ymax ; x0 = xmin ; x1 = xmax ;
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        psrc = &MRIvox(mri_src, xmin, y, z) ;
        pdst = &MRIvox(mri_dst, xmin, y, z) ;
        for (x = xmin ; x <= xmax ; x++, psrc++, pdst++)
        {
          if (z == 130 && (y == 145 || y == 146) && x == 142)
            DiagBreak() ;
          val = *psrc ;
          if ((val > threshold) || (*pdst == fill_val))
            continue ;

          on = 0 ;
          if (x > 0)
            on = (MRIvox(mri_dst, x-1, y, z) == fill_val) ;
          if (!on && (x < width-1))
            on = (MRIvox(mri_dst, x+1, y, z) == fill_val) ;
          if (!on && (y > 0))
            on = (MRIvox(mri_dst, x, y-1, z) == fill_val) ;
          if (!on && (y < height-1))
            on = (MRIvox(mri_dst, x, y+1, z) == fill_val) ;
          if (!on && (z > 0))
            on = (MRIvox(mri_dst, x, y, z-1) == fill_val) ;
          if (!on && (z < depth-1))
            on = (MRIvox(mri_dst, x, y, z+1) == fill_val) ;
          if (on)
          {
            if (z == 130 && (y == 145 || y == 146) && x == 142)
              DiagBreak() ;
            if (x <= x0)
              x0 = MAX(x-1,0) ;
            if (x >= x1)
              x1 = MIN(x+1,width-1) ;
            if (y <= y0)
              y0 = MAX(y-1,0) ;
            if (y >= y1)
              y1 = MIN(y+1,height-1) ;
            if (z <= z0)
              z0 = MAX(z-1,0) ;
            if (z >= z1)
              z1 = MIN(z+1,depth-1) ;
            nfilled++ ;
            *pdst = fill_val ;
          }
        }
      }
    }
    zmin = z0 ; zmax = z1 ; ymin = y0 ; ymax = y1 ; xmin = x0 ; xmax = x1 ;
    total_pix += nfilled ;
/*    fprintf(stderr, "# filled = %d\n", nfilled) ;*/
  } while (nfilled > 0) ;

  if (npix)
    *npix = total_pix ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsetResolution(MRI *mri, float xres, float yres, float zres)
{
  mri->ps = (xres+yres+zres) / 3.0f ;
  mri->xsize = xres ;
  mri->ysize = yres ;
  mri->zsize = zres ;
  mri->xstart *= xres ; mri->xend *= xres ;
  mri->ystart *= yres ; mri->yend *= yres ;
  mri->zstart *= zres ; mri->zend *= zres ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsetTransform(MRI *mri,   General_transform *transform) 
{
  mri->transform = *transform ;
  mri->linear_transform = get_linear_transform_ptr(transform) ;
  mri->inverse_linear_transform = get_inverse_linear_transform_ptr(transform) ;
  
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIextractTalairachPlane(MRI *mri_src, MRI *mri_dst, int orientation, 
                                int x, int y, int z, int wsize)
{
  Real     e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi ;
  Real     ex, ey, ez, len ;

  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = x-whalf*mri_dst->xsize ;
    mri_dst->ystart = y-whalf*mri_dst->ysize ;
    mri_dst->zstart = z-whalf*mri_dst->zsize ;
    mri_dst->xend = mri_dst->xstart + wsize*mri_dst->xsize ;
    mri_dst->yend = mri_dst->ystart + wsize*mri_dst->ysize ;
    mri_dst->zend = mri_dst->zstart + wsize*mri_dst->zsize ;
    mri_dst->imnr0 = z + mri_src->imnr0 ;
    mri_dst->imnr1 = mri_dst->imnr0 ;
  }

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    /* the 'x' basis vector */
    ex = (Real)x+1 ; ey = (Real)y ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y+1 ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    /* the 'x' basis vector */
    ex = (Real)x+1 ; ey = (Real)y ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y ; ez = (Real)z+1 ;
    MRIvoxelToTalairachVoxel(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    /* the 'x' basis vector */
    ex = (Real)x ; ey = (Real)y ; ez = (Real)z+1.0 ;
    MRIvoxelToTalairachVoxel(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y+1.0 ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  }

  len = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
/*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  len = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
/*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      MRIvox(mri_dst, xk+whalf,yk+whalf,0) = MRIvox(mri_src, xi, yi, zi) ;
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIextractArbitraryPlane(MRI *mri_src, MRI *mri_dst, 
                         Real e1_x, Real e1_y, Real e1_z, 
                         Real e2_x, Real e2_y, Real e2_z, 
                         int x, int y, int z, int wsize)
{
  Real     xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi ;

  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = x-whalf*mri_dst->xsize ;
    mri_dst->ystart = y-whalf*mri_dst->ysize ;
    mri_dst->zstart = z-whalf*mri_dst->zsize ;
    mri_dst->xend = mri_dst->xstart + wsize*mri_dst->xsize ;
    mri_dst->yend = mri_dst->ystart + wsize*mri_dst->ysize ;
    mri_dst->zend = mri_dst->zstart + wsize*mri_dst->zsize ;
    mri_dst->imnr0 = z + mri_src->imnr0 ;
    mri_dst->imnr1 = mri_dst->imnr0 ;
  }



  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      MRIvox(mri_dst, xk+whalf,yk+whalf,0) = MRIvox(mri_src, xi, yi, zi) ;
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIeraseTalairachPlaneNew(MRI *mri, MRI *mri_mask, int orientation, int x, 
                          int y, int z, int wsize, int fill_val)
{
  Real     e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi, xki, yki, x0, y0 ;
  Real     ex, ey, ez, len ;

  whalf = (wsize-1)/2 ;

  x0 = mri_mask->width/2 ; y0 = mri_mask->height/2 ;
  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    /* the 'x' basis vector */
    ex = (Real)x+1 ; ey = (Real)y ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y+1 ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    /* the 'x' basis vector */
    ex = (Real)x+1 ; ey = (Real)y ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y ; ez = (Real)z+1 ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    /* the 'x' basis vector */
    ex = (Real)x ; ey = (Real)y ; ez = (Real)z+1.0 ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y+1.0 ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  }

/* 
   don't want to normalize basis - they are orthonormal in magnet space,
   not necessarily Talairach space.
*/
  len = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
/*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  len = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
/*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri->xi[nint(xbase + xk*e1_x)] ;
      yi = mri->yi[nint(ybase + xk*e1_y)] ;
      zi = mri->zi[nint(zbase + xk*e1_z)] ;
      xki = mri_mask->xi[xk+x0] ; yki = mri_mask->yi[yk+y0] ;
      if (MRIvox(mri_mask, xki, yki,0))
        MRIvox(mri, xi, yi, zi) = fill_val ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIeraseTalairachPlane(MRI *mri, MRI *mri_mask, int orientation, int x, int y, 
                       int z, int wsize, int fill_val)
{
  Real     e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi ;
  Real     ex, ey, ez, len ;

  whalf = (wsize-1)/2 ;

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    /* the 'x' basis vector */
    ex = (Real)x+1 ; ey = (Real)y ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y+1 ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    /* the 'x' basis vector */
    ex = (Real)x+1 ; ey = (Real)y ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y ; ez = (Real)z+1 ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    /* the 'x' basis vector */
    ex = (Real)x ; ey = (Real)y ; ez = (Real)z+1.0 ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (Real)x ; e1_y -= (Real)y ; e1_z -= (Real)z ; 

    /* the 'y' basis vector */
    ex = (Real)x ; ey = (Real)y+1.0 ; ez = (Real)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (Real)x ; e2_y -= (Real)y ; e2_z -= (Real)z ; 
    break ;
  }

/* 
   don't want to normalize basis - they are orthonormal in magnet space,
   not necessarily Talairach space.
*/
  len = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
/*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  len = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
/*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri->xi[nint(xbase + xk*e1_x)] ;
      yi = mri->yi[nint(ybase + xk*e1_y)] ;
      zi = mri->zi[nint(zbase + xk*e1_z)] ;
      if (MRIvox(mri_mask, xk+whalf,yk+whalf,0))
        MRIvox(mri, xi, yi, zi) = fill_val ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIextractPlane(MRI *mri_src, MRI *mri_dst, int orientation, int where)
{
  int      x, y, z, width, height ;

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    width = mri_src->width ; height = mri_src->height ;
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    width = mri_src->width ; height = mri_src->depth ;
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    width = mri_src->depth ; height = mri_src->height ;
    break ;
  }

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->zstart = where ; mri_dst->zend = where+1 ;
    mri_dst->imnr0 = 0 ; mri_dst->imnr1 = 1 ;
  }
  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
        MRIvox(mri_dst, x, y, 0) = MRIvox(mri_src, x, y, where) ;
    }
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
        MRIvox(mri_dst, x, z, 0) = MRIvox(mri_src, x, where, z) ;
    }
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    for (z = 0 ; z < mri_src->depth ; z++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
        MRIvox(mri_dst, z, y, 0) = MRIvox(mri_src, where, y, z) ;
    }
    break ;
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIerasePlane(MRI *mri, float x0, float y0, float z0,
              float dx, float dy, float dz, int fill_val)
{
  int      *pxi, *pyi, *pzi, xi, yi, zi, x, y, z ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, maxlen, l1, l2 ;

  maxlen = MAX(mri->width, mri->height) ; maxlen = MAX(maxlen, mri->depth) ; 

  /* don't care about sign (right-hand rule) */
  e1_x = dz*dz - dx*dy ; 
  e1_y = dx*dx - dy*dz ;
  e1_z = dy*dy - dz*dx ;
  l1 = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
  e1_x /= l1 ; e1_y /= l1 ; e1_z /= l1 ;

  e2_x = e1_y*dz - e1_z*dy ;  
  e2_y = e1_x*dz - e1_z*dx ;
  e2_z = e1_y*dx - e1_x*dy ;
  l2 = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
  e2_x /= l2 ; e2_y /= l2 ; e2_z /= l2 ;

  pxi = mri->xi ; pyi = mri->yi ; pzi = mri->zi ;
  maxlen *= 1.5 ;  /* make sure to get entire extent */
  for (l1 = -maxlen/2 ; l1 <= maxlen/2 ; l1 += 0.5f)
  {
    for (l2 = -maxlen/2 ; l2 <= maxlen/2 ; l2 += 0.5f)
    {
      x = nint(x0 + l1 * e1_x + l2 * e2_x) ; xi = pxi[x] ;
      y = nint(y0 + l1 * e1_y + l2 * e2_y) ; yi = pyi[y] ;
      z = nint(z0 + l1 * e1_z + l2 * e2_z) ; zi = pzi[z] ;
      MRIvox(mri, xi, yi, zi) = fill_val ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume to cubic voxels.
------------------------------------------------------*/
MRI *
MRIinterpolate(MRI *mri_src, MRI *mri_dst)
{
  int     xs, ys, zs, xd, yd, zd, max_dim, xorig, yorig, zorig, dorig ;
  float   sx, sy, sz, psize ;
  int     width, height, depth, i ;
  float   xsmd, ysmd, zsmd, xspd, yspd, zspd, weights[8], fout ;
  int     xsp, xsm, ysp, ysm, zsp, zsm ;  /* surrounding coordinates */
  float vals[8], outval ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth;
  if (width > height)
  {
    max_dim = width > depth ? width : depth ;
    psize = width > depth ? mri_src->xsize : mri_src->zsize ;
  }
  else
  {
    max_dim = height > depth ? height : depth ;
    psize = height > depth ? mri_src->ysize : mri_src->zsize ;
  }

  if (!mri_dst)
  {
    mri_dst = MRIalloc(max_dim, max_dim, max_dim, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  mri_dst->xsize = mri_dst->ysize = mri_dst->zsize = psize ;
  sx = (float)width / (float)max_dim  ;
  sy = (float)height / (float)max_dim  ;
  sz = (float)depth / (float)max_dim  ;

  mri_dst->imnr0 = 1;
  mri_dst->imnr1 = mri_dst->depth;
  mri_dst->xstart = -mri_dst->xsize * mri_dst->width / 2; 
  mri_dst->xend = -mri_dst->xstart;
  mri_dst->ystart = -mri_dst->ysize * mri_dst->height / 2;
  mri_dst->yend = -mri_dst->ystart;
  mri_dst->zstart = -mri_dst->zsize * mri_dst->depth / 2;
  mri_dst->zend = -mri_dst->zstart;

  xorig = (width-1)/2 ; yorig = (height-1)/2 ; zorig = (depth-1)/2 ;
  dorig = (max_dim-1)/2 ;

/*
  for each output voxel, find the 8 nearest input voxels and interpolate
  the output voxel from them.
  */
  for (zd = 0 ; zd < max_dim ; zd++)
  {
printf("interpolate: %d/%d\n",zd+1,max_dim);
    for (yd = 0 ; yd < max_dim ; yd++)
    {
      for (xd = 0 ; xd < max_dim ; xd++)
      {
        /* do trilinear interpolation here */
        xs = sx*(xd-dorig) + xorig ;
        ys = sy*(yd-dorig) + yorig ;
        zs = sz*(zd-dorig) + zorig ;

        /* 
           these boundary conditions will cause reflection across the border
           for the 1st negative pixel.
           */
        if (xs > -1 && xs < width &&
            ys > -1 && ys < height &&
            zs > -1 && zs < depth)
        {
          xsm = (int)xs ;
          xsp = MIN(width-1, xsm+1) ;
          ysm = (int)ys ;
          ysp = MIN(height-1, ysm+1) ;
          zsm = (int)zs ;
          zsp = MIN(depth-1, zsm+1) ;
          xsmd = xs - (float)xsm ;
          ysmd = ys - (float)ysm ;
          zsmd = zs - (float)zsm ;
          xspd = (1.0f - xsmd) ;
          yspd = (1.0f - ysmd) ;
          zspd = (1.0f - zsmd) ;

/*          vals[0] = mri_src->slices[zsm][ysm][xsm] ;
          vals[1] = mri_src->slices[zsm][ysm][xsp] ;
          vals[2] = mri_src->slices[zsm][ysp][xsm] ;
          vals[3] = mri_src->slices[zsm][ysp][xsp] ;
          vals[4] = mri_src->slices[zsp][ysm][xsm] ;
          vals[5] = mri_src->slices[zsp][ysm][xsp] ;
          vals[6] = mri_src->slices[zsp][ysp][xsm] ;
          vals[7] = mri_src->slices[zsp][ysp][xsp] ;
*/
/* different vox types... */
          vals[0] = MRIFvox(mri_src, xsm, ysm, zsm);
          vals[1] = MRIFvox(mri_src, xsp, ysm, zsm);
          vals[2] = MRIFvox(mri_src, xsm, ysp, zsm);
          vals[3] = MRIFvox(mri_src, xsp, ysp, zsm);
          vals[4] = MRIFvox(mri_src, xsm, ysm, zsp);
          vals[5] = MRIFvox(mri_src, xsp, ysm, zsp);
          vals[6] = MRIFvox(mri_src, xsm, ysp, zsp);
          vals[7] = MRIFvox(mri_src, xsp, ysp, zsp);

          weights[0] = zsmd * ysmd * xsmd ;
          weights[1] = zsmd * ysmd * xspd ;
          weights[2] = zsmd * yspd * xsmd ;
          weights[3] = zsmd * yspd * xspd ;
          weights[4] = zspd * ysmd * xsmd ;
          weights[5] = zspd * ysmd * xspd ;
          weights[6] = zspd * yspd * xsmd ;
          weights[7] = zspd * yspd * xspd ;
/*
for(i = 0;i < 8;i++)
  printf("%d, %f, %f\n",i,vals[i], weights[i]);
*/
          for (fout = 0.0f, i = 0 ; i < 8 ; i++)
            fout += (float)vals[i] * weights[i] ;
          outval = (float)nint(fout) ;
          MRIvox(mri_dst, xd, yd, zd) = (BUFTYPE)nint(fout) ;

        }
      }
    }
  }

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsampleVolumeFrame(MRI *mri,Real x,Real y,Real z,int frame,Real *pval)
{
  int  OutOfBounds;
  int  xm, xp, ym, yp, zm, zp, width, height, depth ;
  Real val, xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  if (frame >= mri->nframes)
  {
    *pval = 1.0 ;
    return(NO_ERROR) ;
  }

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if(OutOfBounds == 1){
    /* unambiguoulsy out of bounds */
    *pval = 0.0;
    return(NO_ERROR) ;
  }

  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRIseq_vox(mri, xm, ym, zm, frame) +
      xpd * ypd * zmd * (Real)MRIseq_vox(mri, xm, ym, zp, frame) +
      xpd * ymd * zpd * (Real)MRIseq_vox(mri, xm, yp, zm, frame) +
      xpd * ymd * zmd * (Real)MRIseq_vox(mri, xm, yp, zp, frame) +
      xmd * ypd * zpd * (Real)MRIseq_vox(mri, xp, ym, zm, frame) +
      xmd * ypd * zmd * (Real)MRIseq_vox(mri, xp, ym, zp, frame) +
      xmd * ymd * zpd * (Real)MRIseq_vox(mri, xp, yp, zm, frame) +
      xmd * ymd * zmd * (Real)MRIseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_FLOAT:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRIFseq_vox(mri, xm, ym, zm, frame) +
      xpd * ypd * zmd * (Real)MRIFseq_vox(mri, xm, ym, zp, frame) +
      xpd * ymd * zpd * (Real)MRIFseq_vox(mri, xm, yp, zm, frame) +
      xpd * ymd * zmd * (Real)MRIFseq_vox(mri, xm, yp, zp, frame) +
      xmd * ypd * zpd * (Real)MRIFseq_vox(mri, xp, ym, zm, frame) +
      xmd * ypd * zmd * (Real)MRIFseq_vox(mri, xp, ym, zp, frame) +
      xmd * ymd * zpd * (Real)MRIFseq_vox(mri, xp, yp, zm, frame) +
      xmd * ymd * zmd * (Real)MRIFseq_vox(mri, xp, yp, zp, frame) ;
    break ;
   case MRI_SHORT:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRISseq_vox(mri, xm, ym, zm, frame) +
      xpd * ypd * zmd * (Real)MRISseq_vox(mri, xm, ym, zp, frame) +
      xpd * ymd * zpd * (Real)MRISseq_vox(mri, xm, yp, zm, frame) +
      xpd * ymd * zmd * (Real)MRISseq_vox(mri, xm, yp, zp, frame) +
      xmd * ypd * zpd * (Real)MRISseq_vox(mri, xp, ym, zm, frame) +
      xmd * ypd * zmd * (Real)MRISseq_vox(mri, xp, ym, zp, frame) +
      xmd * ymd * zpd * (Real)MRISseq_vox(mri, xp, yp, zm, frame) +
      xmd * ymd * zmd * (Real)MRISseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_INT:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRIIseq_vox(mri, xm, ym, zm, frame) +
      xpd * ypd * zmd * (Real)MRIIseq_vox(mri, xm, ym, zp, frame) +
      xpd * ymd * zpd * (Real)MRIIseq_vox(mri, xm, yp, zm, frame) +
      xpd * ymd * zmd * (Real)MRIIseq_vox(mri, xm, yp, zp, frame) +
      xmd * ypd * zpd * (Real)MRIIseq_vox(mri, xp, ym, zm, frame) +
      xmd * ypd * zmd * (Real)MRIIseq_vox(mri, xp, ym, zp, frame) +
      xmd * ymd * zpd * (Real)MRIIseq_vox(mri, xp, yp, zm, frame) +
      xmd * ymd * zmd * (Real)MRIIseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_LONG:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRILseq_vox(mri, xm, ym, zm, frame) +
      xpd * ypd * zmd * (Real)MRILseq_vox(mri, xm, ym, zp, frame) +
      xpd * ymd * zpd * (Real)MRILseq_vox(mri, xm, yp, zm, frame) +
      xpd * ymd * zmd * (Real)MRILseq_vox(mri, xm, yp, zp, frame) +
      xmd * ypd * zpd * (Real)MRILseq_vox(mri, xp, ym, zm, frame) +
      xmd * ypd * zmd * (Real)MRILseq_vox(mri, xp, ym, zp, frame) +
      xmd * ymd * zpd * (Real)MRILseq_vox(mri, xp, yp, zm, frame) +
      xmd * ymd * zmd * (Real)MRILseq_vox(mri, xp, yp, zp, frame) ;
    break ;
 default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIsampleVolumeFrame: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}



/*---------------------------------------------------------------------------
   Purpose: to return the approximate fraction of a voxel centered 
     at the given point
     is labeled with a given label by the labeled volume: mri

   Input: mri  is the labeled volume. Ea voxel contains the uchar label index
          x,y,z is the floating point location of the center of a voxel 
          whose labeling is to be determined. The voxel is examined to see 
          how much of it is labeled with the label, ucharLabel

   Output: pval is the fraction which the given voxel location is labeled 
           by ucharLabel 
           returns NO_ERROR or ERROR_UNSUPPORTED if an unsupported 
           (non-uchar) mri labeled volume is passed in
   AAM: 7/26/00
--------------------------------------------------------------------------*/

#ifndef uchar
#define uchar unsigned char
#endif
int
MRIsampleLabeledVolume(MRI *mri, Real x, Real y, Real z, Real *pval, unsigned char ucharLabel)
{
  /* m's are the mri grid locations less than x (or y or z)
     (i.e. floor(x), p's are essentially rounding up  to the next 
     grid location  greater than x  */
  int  xm, xp, ym, yp, zm, zp;  
  int width, height, depth ;
  Real xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */
  uchar ucharDmmm;
  uchar ucharDmmp;
  uchar ucharDmpm;
  uchar ucharDmpp;
  uchar ucharDpmm;
  uchar ucharDpmp;
  uchar ucharDppm;
  uchar ucharDppp;


  *pval=0.0;



  width = mri->width ; 
  height = mri->height ; 
  depth = mri->depth ; 
  /*  if (x >= width) 
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;
  */

  /* if the x,y,z point is out of range then return that none of the given voxel was labeled by ucharLabel */
  if (x >= width) return(NO_ERROR) ;
  if (y >= height) return(NO_ERROR) ;
  if (z >= depth) return(NO_ERROR) ;
  if (x < 0.0) return(NO_ERROR) ;
  if (y < 0.0) return(NO_ERROR) ;
  if (z < 0.0) return(NO_ERROR) ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  ucharDmmm = MRIvox(mri, xm, ym, zm) == ucharLabel ? 1 : 0;
  ucharDmmp = MRIvox(mri, xm, ym, zp) == ucharLabel ? 1 : 0;
  ucharDmpm = MRIvox(mri, xm, yp, zm) == ucharLabel ? 1 : 0;
  ucharDmpp = MRIvox(mri, xm, yp, zp) == ucharLabel ? 1 : 0;
  ucharDpmm = MRIvox(mri, xp, ym, zm) == ucharLabel ? 1 : 0;
  ucharDpmp = MRIvox(mri, xp, ym, zp) == ucharLabel ? 1 : 0;
  ucharDppm = MRIvox(mri, xp, yp, zm) == ucharLabel ? 1 : 0;
  ucharDppp = MRIvox(mri, xp, yp, zp) == ucharLabel ? 1 : 0;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval =  
      xpd * ypd * zpd * (Real)ucharDmmm +
      xpd * ypd * zmd * (Real)ucharDmmp +
      xpd * ymd * zpd * (Real)ucharDmpm +
      xpd * ymd * zmd * (Real)ucharDmpp +
      xmd * ypd * zpd * (Real)ucharDpmm +
      xmd * ypd * zmd * (Real)ucharDpmp +
      xmd * ymd * zpd * (Real)ucharDppm +
      xmd * ymd * zmd * (Real)ucharDppp ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIsampleVolume: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsampleVolumeType(MRI *mri, Real x, Real y, Real z, Real *pval, int type)
{
  int   xv, yv, zv ;
  int OutOfBounds;

  switch (type)
  {
  default:
  case SAMPLE_NEAREST:
    break ;
  case SAMPLE_TRILINEAR:
    return(MRIsampleVolume(mri, x, y, z, pval)) ;
  case SAMPLE_SINC:
    return(MRIsincSampleVolume(mri, x, y, z, 5, pval)) ;
  }

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if(OutOfBounds == 1){
    /* unambiguoulsy out of bounds */
    *pval = 0.0;
    return(NO_ERROR) ;
  }

  xv = nint(x) ; yv = nint(y) ; zv = nint(z) ; 
  if (xv < 0)
    xv = 0 ;
  if (xv >= mri->width)
    xv = mri->width-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= mri->height)
    yv = mri->height-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= mri->depth)
    zv = mri->depth-1 ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = (float)MRIvox(mri, xv, yv, zv) ;
    break ;
  case MRI_SHORT:
    *pval = (float)MRISvox(mri, xv, yv, zv) ;
    break ;
  case MRI_INT:
    *pval = (float)MRIIvox(mri, xv, yv, zv) ;
    break ;
  case MRI_FLOAT:
    *pval = MRIFvox(mri, xv, yv, zv) ;
    break ;
  default:
    *pval = 0 ;
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, "MRIsampleVolumeType: unsupported volume type %d",
                 mri->type)) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIsampleVolume(MRI *mri, Real x, Real y, Real z, Real *pval)
{
  int  OutOfBounds;
  int  xm, xp, ym, yp, zm, zp, width, height, depth ;
  Real val, xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if(OutOfBounds == 1){
    /* unambiguoulsy out of bounds */
    *pval = 0.0;
    return(NO_ERROR) ;
  }

  width = mri->width ; height = mri->height ; depth = mri->depth ; 

  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRIvox(mri, xm, ym, zm) +
      xpd * ypd * zmd * (Real)MRIvox(mri, xm, ym, zp) +
      xpd * ymd * zpd * (Real)MRIvox(mri, xm, yp, zm) +
      xpd * ymd * zmd * (Real)MRIvox(mri, xm, yp, zp) +
      xmd * ypd * zpd * (Real)MRIvox(mri, xp, ym, zm) +
      xmd * ypd * zmd * (Real)MRIvox(mri, xp, ym, zp) +
      xmd * ymd * zpd * (Real)MRIvox(mri, xp, yp, zm) +
      xmd * ymd * zmd * (Real)MRIvox(mri, xp, yp, zp) ;
    break ;
  case MRI_FLOAT:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRIFvox(mri, xm, ym, zm) +
      xpd * ypd * zmd * (Real)MRIFvox(mri, xm, ym, zp) +
      xpd * ymd * zpd * (Real)MRIFvox(mri, xm, yp, zm) +
      xpd * ymd * zmd * (Real)MRIFvox(mri, xm, yp, zp) +
      xmd * ypd * zpd * (Real)MRIFvox(mri, xp, ym, zm) +
      xmd * ypd * zmd * (Real)MRIFvox(mri, xp, ym, zp) +
      xmd * ymd * zpd * (Real)MRIFvox(mri, xp, yp, zm) +
      xmd * ymd * zmd * (Real)MRIFvox(mri, xp, yp, zp) ;
    break ;
  case MRI_SHORT:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRISvox(mri, xm, ym, zm) +
      xpd * ypd * zmd * (Real)MRISvox(mri, xm, ym, zp) +
      xpd * ymd * zpd * (Real)MRISvox(mri, xm, yp, zm) +
      xpd * ymd * zmd * (Real)MRISvox(mri, xm, yp, zp) +
      xmd * ypd * zpd * (Real)MRISvox(mri, xp, ym, zm) +
      xmd * ypd * zmd * (Real)MRISvox(mri, xp, ym, zp) +
      xmd * ymd * zpd * (Real)MRISvox(mri, xp, yp, zm) +
      xmd * ymd * zmd * (Real)MRISvox(mri, xp, yp, zp) ;
    break ;
  case MRI_INT:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRIIvox(mri, xm, ym, zm) +
      xpd * ypd * zmd * (Real)MRIIvox(mri, xm, ym, zp) +
      xpd * ymd * zpd * (Real)MRIIvox(mri, xm, yp, zm) +
      xpd * ymd * zmd * (Real)MRIIvox(mri, xm, yp, zp) +
      xmd * ypd * zpd * (Real)MRIIvox(mri, xp, ym, zm) +
      xmd * ypd * zmd * (Real)MRIIvox(mri, xp, ym, zp) +
      xmd * ymd * zpd * (Real)MRIIvox(mri, xp, yp, zm) +
      xmd * ymd * zmd * (Real)MRIIvox(mri, xp, yp, zp) ;
    break ;
  case MRI_LONG:
    *pval = val = 
      xpd * ypd * zpd * (Real)MRILvox(mri, xm, ym, zm) +
      xpd * ypd * zmd * (Real)MRILvox(mri, xm, ym, zp) +
      xpd * ymd * zpd * (Real)MRILvox(mri, xm, yp, zm) +
      xpd * ymd * zmd * (Real)MRILvox(mri, xm, yp, zp) +
      xmd * ypd * zpd * (Real)MRILvox(mri, xp, ym, zm) +
      xmd * ypd * zmd * (Real)MRILvox(mri, xp, ym, zp) +
      xmd * ymd * zpd * (Real)MRILvox(mri, xp, yp, zm) +
      xmd * ymd * zmd * (Real)MRILvox(mri, xp, yp, zp) ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED, 
                (ERROR_UNSUPPORTED, 
                 "MRIsampleVolume: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define IMIN(a,b) (a < b ? a : b)
#define IMAX(a,b) (a > b ? a : b)
double ham_sinc(double x,double fullwidth)
{
  double ham;
  if( fabs(x) < 1.0e-5)
    ham = 1.0;
  else { 
    ham = sin(PI*x)/(PI*x);
    ham *= 0.54 + 0.46 * cos(2.0*PI*x/fullwidth);
  }
  return ham;
}
/*-------------------------------------------------------------------------*/
int
MRIsincSampleVolume(MRI *mri, Real x, Real y, Real z, int hw, Real *pval)
{

  int  OutOfBounds;
  int  width, height, depth ;
  int nwidth; 
  int ix_low,ix_high,iy_low,iy_high,iz_low,iz_high;
  int jx1,jy1,jz1,jx_rel,jy_rel,jz_rel;
  double coeff_x[128],coeff_y[128],coeff_z[128];
  double coeff_x_sum,coeff_y_sum,coeff_z_sum;
  double sum_x,sum_y,sum_z;
  double xsize,ysize,zsize;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if(OutOfBounds == 1){
    /* unambiguoulsy out of bounds */
    *pval = 0.0;
    return(NO_ERROR) ;
  }

  xsize = mri->xsize; ysize=mri->ysize; zsize=mri->zsize;
  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  nwidth = hw;
  ix_low = floor((double)x);
  ix_high = ceil((double)x);
  iy_low = floor((double)y);
  iy_high = ceil((double)y);
  iz_low = floor((double)z);
  iz_high = ceil((double)z);

  coeff_x_sum = coeff_y_sum = coeff_z_sum = 0; 
  if(iz_low>=0 && iz_high < depth) {
    for (jx1=IMAX(ix_high-nwidth,0), jx_rel=0;
   jx1<IMIN(ix_low+nwidth,width-1); 
   jx1++,jx_rel++) {
      coeff_x[jx_rel] = ham_sinc((double)(x-jx1)*xsize,2*xsize*nwidth);
      coeff_x_sum += coeff_x[jx_rel];
    }
    for (jy1=IMAX(iy_high-nwidth,0), jy_rel=0;
   jy1<IMIN(iy_low+nwidth,height-1); 
   jy1++,jy_rel++) {
      coeff_y[jy_rel] = ham_sinc((double)(y-jy1)*ysize,2*nwidth*ysize);
      coeff_y_sum += coeff_y[jy_rel];
    }
    for (jz1=IMAX(iz_high-nwidth,0), jz_rel=0;
   jz1<IMIN(iz_low+nwidth,depth-1); 
   jz1++,jz_rel++) {
      coeff_z[jz_rel] = ham_sinc((double)(z-jz1)*zsize,2*nwidth*zsize);
      coeff_z_sum += coeff_z[jz_rel];
    }
    
    for(sum_z=0., jz1=IMAX(iz_high-nwidth,0), jz_rel = 0;
  jz1 < IMIN(iz_low+nwidth,depth-1);
  jz1++, jz_rel++) {
      
      for(sum_y=0., jy1=IMAX(iy_high-nwidth,0), jy_rel = 0;
    jy1 < IMIN(iy_low+nwidth,height-1);
    jy1++, jy_rel++) {
  for(sum_x=0., jx1=IMAX(ix_high-nwidth,0), jx_rel = 0;
      jx1 < IMIN(ix_low+nwidth,width-1);
      jx1++, jx_rel++) {
    
    switch(mri->type)
      {
      case MRI_UCHAR:
        sum_x += (coeff_x[jx_rel]/coeff_x_sum) 
    * (double)MRIvox(mri,jx1,jy1,jz1);
        break;
      case MRI_FLOAT:
        sum_x += (coeff_x[jx_rel]/coeff_x_sum) 
    * (double)MRIFvox(mri,jx1,jy1,jz1);
        break; 
      case MRI_SHORT:
        sum_x += (coeff_x[jx_rel]/coeff_x_sum) 
    * (double)MRISvox(mri,jx1,jy1,jz1);
        break;
      default:
        ErrorReturn(ERROR_UNSUPPORTED, 
        (ERROR_UNSUPPORTED, 
         "MRIsincSampleVolume: unsupported type %d", 
         mri->type)) ;
        break;
      }
  }
  sum_y += sum_x * (coeff_y[jy_rel]/coeff_y_sum);
      }
      sum_z += sum_y * (coeff_z[jz_rel]/coeff_z_sum); 
    }
    if((mri->type == MRI_UCHAR || mri->type == MRI_SHORT) && sum_z<0.0)
      *pval = 0.0;
    else if(mri->type == MRI_UCHAR && sum_z >255.0)
      *pval = 255.0;
    else if(mri->type == MRI_SHORT && sum_z > 65535.0)
      *pval = 65535.0;
    else
      *pval = sum_z;
  } else 
    *pval = 0.0;
  
  return(NO_ERROR);
}
/*-----------------------------------------------------------------
  MRIindexNotInVolume() - determines whether a col, row, slice point is
  in the mri volume. If it is unambiguously in the volume, then 0
  is returned. If it is within 0.5 of the edge of the volume, -1
  is returned. Otherwise 1 is returned. Flagging the case where
  the point is within 0.5 of the edge can be used for assigning
  a nearest neighbor when the point is outside but close to the 
  volume. In this case the index of the nearest neighbor can safely 
  be computed as the nearest integer to col, row, and slice.
  -----------------------------------------------------------------*/
int MRIindexNotInVolume(MRI *mri, Real col, Real row, Real slice) 
{
  float nicol, nirow, nislice;

  /* unambiguously in the volume */
  if(col   >= 0 && col   <= mri->width-1  &&
     row   >= 0 && row   <= mri->height-1 &&
     slice >= 0 && slice <= mri->depth-1  ) 
    return(0);

  /* within 0.5 of the edge of the volume */
  nicol   = rint(col);
  nirow   = rint(row);
  nislice = rint(slice);
  if(nicol   >= 0 && nicol   < mri->width  &&
     nirow   >= 0 && nirow   < mri->height &&
     nislice >= 0 && nislice < mri->depth  ) 
    return(-1);

  /* unambiguously NOT in the volume */
  return(1);

}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume directional derivative using
          trilinear interpolation.
------------------------------------------------------*/
float
MRIsampleCardinalDerivative(MRI *mri, int x, int y, int z,
                            int xk, int yk, int zk)
{
  float d ;

  if (xk)
    d = MRIsampleXDerivative(mri, x, y, z, xk) ;
  else if (yk)
    d = MRIsampleYDerivative(mri, x, y, z, yk) ;
  else
    d = MRIsampleZDerivative(mri, x, y, z, zk) ;
  return(d) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume directional derivative using
          trilinear interpolation.
------------------------------------------------------*/
float 
MRIsampleXDerivative(MRI *mri, int x, int y, int z, int dir)
{
  float dx ;
  int   yk, zk, xi, yi, zi, nvox ;

  dx = 0.0 ;

  xi = mri->xi[x+dir] ;
  for (nvox = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dx += dir*MRIvox(mri, xi, yi, zi) ;  /* x+dir */
      dx -= dir*MRIvox(mri, x, yi, zi) ;   /* - x */
      nvox += 2 ;
    }
  }
  dx /= ((float)nvox*mri->xsize) ;
  return(dx) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume directional derivative using
          trilinear interpolation.
------------------------------------------------------*/
float 
MRIsampleYDerivative(MRI *mri, int x, int y, int z, int dir)
{
  float dy ;
  int   xk, zk, xi, yi, zi, nvox ;

  dy = 0.0 ;

  yi = mri->yi[y+dir] ;
  for (nvox = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (xk = -1 ; xk <= 1 ; xk++)
    {
      xi = mri->xi[x+xk] ;
      dy += dir*MRIvox(mri, xi, yi, zi) ;  /* x+dir */
      dy -= dir*MRIvox(mri, x, yi, zi) ;   /* - x */
      nvox += 2 ;
    }
  }
  dy /= ((float)nvox*mri->ysize) ;
  return(dy) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume directional derivative using
          trilinear interpolation.
------------------------------------------------------*/
float 
MRIsampleZDerivative(MRI *mri, int x, int y, int z, int dir)
{
  float dz ;
  int   xk, yk, xi, yi, zi, nvox ;

  dz = 0.0 ;

  zi = mri->zi[z+dir] ;
  for (nvox = 0, xk = -1 ; xk <= 1 ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dz += dir*MRIvox(mri, xi, yi, zi) ;  /* x+dir */
      dz -= dir*MRIvox(mri, x, yi, zi) ;   /* - x */
      nvox += 2 ;
    }
  }
  dz /= ((float)nvox*mri->zsize) ;
  return(dz) ;
}
int   
MRIsampleVolumeDirectionScale(MRI *mri, Real x, Real y, Real z,
                              Real dx, Real dy, Real dz, Real *pmag,
                              double sigma)
{
  int  width, height, depth ;
  Real xp1, xm1, yp1, ym1, zp1, zm1, len ;
  Real dist, val, k, ktotal, step_size, total_val ;
  int  n ;

  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  step_size = MAX(.5,sigma/2) ;
  for (total_val = ktotal = 0.0,n = 0, len = 0.0, dist = step_size ; 
       dist <= MAX(2*sigma,step_size); 
       dist += step_size, n++)
  {
    if (FZERO(sigma))
      k = 1.0 ;
    else
      k = exp(-dist*dist/(2*sigma*sigma)) ; 
    ktotal += k ;
    len += dist ;
    xp1 = x + dist*dx ; yp1 = y + dist*dy ; zp1 = z + dist*dz ;
    MRIsampleVolume(mri, xp1, yp1, zp1, &val) ;
    total_val += k*val ;
    
    xm1 = x - dist*dx ; ym1 = y - dist*dy ; zm1 = z - dist*dz ;
    MRIsampleVolume(mri, xm1, ym1, zm1, &val) ;
    total_val += k*val ;
    if (FZERO(step_size))
      break ;
  }
  total_val /= (double)2.0*ktotal ; 

  *pmag = total_val ;
  return(NO_ERROR) ;
}

int
MRIsampleVolumeDerivativeScale(MRI *mri, Real x, Real y, Real z, Real dx, 
                               Real dy, Real dz, Real *pmag, double sigma)
{
  int  width, height, depth ;
  Real xp1, xm1, yp1, ym1, zp1, zm1, vp1, vm1, len ;
  Real dist, val, k, ktotal, step_size ;
  int  n ;

  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  step_size = MAX(.5,sigma/2) ;
  for (ktotal = 0.0,n = 0, len = vp1 = vm1 = 0.0, dist = step_size ; 
       dist <= MAX(2*sigma,step_size); 
       dist += step_size, n++)
  {
    if (FZERO(sigma))
      k = 1.0 ;
    else
      k = exp(-dist*dist/(2*sigma*sigma)) ; 
    ktotal += k ;
    len += dist ;
    xp1 = x + dist*dx ; yp1 = y + dist*dy ; zp1 = z + dist*dz ;
    MRIsampleVolume(mri, xp1, yp1, zp1, &val) ;
    vp1 += k*val ;
    
    xm1 = x - dist*dx ; ym1 = y - dist*dy ; zm1 = z - dist*dz ;
    MRIsampleVolume(mri, xm1, ym1, zm1, &val) ;
    vm1 += k*val ;
    if (FZERO(step_size))
      break ;
  }
  vm1 /= (double)ktotal ; vp1 /= (double)ktotal ; len /= (double)ktotal ;

  *pmag = (vp1-vm1) / (2.0*len) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume directional derivative using
          trilinear interpolation.
------------------------------------------------------*/
int
MRIsampleVolumeDerivative(MRI *mri, Real x, Real y, Real z,
                          Real dx, Real dy, Real dz, Real *pmag)
{
  int  width, height, depth ;
  Real xp1, xm1, yp1, ym1, zp1, zm1, vp1, vm1, len ;

  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;
#if 1
    {
      Real dist, val ;
      int  n ;

      for (n = 0, len = vp1 = vm1 = 0.0, dist = .5 ; dist <= 2 ; 
           dist += 0.5, n++)
      {
        len += dist ;
        xp1 = x + dist*dx ; yp1 = y + dist*dy ; zp1 = z + dist*dz ;
        MRIsampleVolume(mri, xp1, yp1, zp1, &val) ;
        vp1 += val ;

        xm1 = x - dist*dx ; ym1 = y - dist*dy ; zm1 = z - dist*dz ;
        MRIsampleVolume(mri, xm1, ym1, zm1, &val) ;
        vm1 += val ;
      }
      vm1 /= (double)n ; vp1 /= (double)n ; len /= (double)n ;
    }
#else
  xp1 = x + dx ; xm1 = x - dx ;
  yp1 = y + dy ; ym1 = y - dy ;
  zp1 = z + dz ; zm1 = z - dz ;
  len = sqrt(dx*dx+dy*dy+dz*dz) ;
  MRIsampleVolume(mri, xp1, yp1, zp1, &vp1) ;
  MRIsampleVolume(mri, xm1, ym1, zm1, &vm1) ;
#endif

  *pmag = (vp1-vm1) / (2.0*len) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Interpolate the volume gradient to cubic voxels.
------------------------------------------------------*/
int
MRIsampleVolumeGradient(MRI *mri, Real x, Real y, Real z, 
                        Real *pdx, Real *pdy, Real *pdz)
{
  int  width, height, depth ;
  Real xp1, xm1, yp1, ym1, zp1, zm1 ;

  width = mri->width ; height = mri->height ; depth = mri->depth ; 
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;
  MRIsampleVolume(mri, x+1.0, y, z, &xp1) ;
  MRIsampleVolume(mri, x-1.0, y, z, &xm1) ;

  MRIsampleVolume(mri, x, y+1.0, z, &yp1) ;
  MRIsampleVolume(mri, x, y-1.0, z, &ym1) ;

  MRIsampleVolume(mri, x, y, z+1.0, &zp1) ;
  MRIsampleVolume(mri, x, y, z-1.0, &zm1) ;

  *pdx = (xp1-xm1)/(2.0*mri->xsize) ; 
  *pdy = (yp1-ym1)/(2.0*mri->ysize) ; 
  *pdz = (zp1-zm1)/(2.0*mri->zsize) ; 
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIneighborsOn(MRI *mri, int x0, int y0, int z0, int min_val)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,mri->xi[x0-1],y0,z0) >= min_val)
    nbrs++ ;
  if (MRIvox(mri,mri->xi[x0+1],y0,z0) >= min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,mri->yi[y0+1],z0) >= min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,mri->yi[y0-1],z0) >= min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,mri->zi[z0+1]) >= min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,mri->zi[z0-1]) >= min_val)
    nbrs++ ;
  return(nbrs) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
MRIneighborsOn3x3(MRI *mri, int x, int y, int z, int min_val)
{
  int   xk, yk, zk, xi, yi, zi, nbrs ;

  for (nbrs = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (!zk && !yk && !xk)
          continue ;
        if (MRIvox(mri, xi, yi, zi) > min_val)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIneighborsOff(MRI *mri, int x0, int y0, int z0, int min_val)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,x0-1,y0,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0+1,y0,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0+1,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0-1,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0+1) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0-1) < min_val)
    nbrs++ ;
  return(nbrs) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIneighborsOff3x3(MRI *mri, int x, int y, int z, int min_val)
{
  int   xk, yk, zk, xi, yi, zi, nbrs ;

  for (nbrs = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (!zk && !yk && !xk)
          continue ;
        if (MRIvox(mri, xi, yi, zi) < min_val)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
/*-----------------------------------------------------
  Perform an linear coordinate transformation x' = Ax on
  the MRI image mri_src into mri_dst
------------------------------------------------------*/
MRI *
MRIinverseLinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA)
{
  MATRIX  *m_inv ;

  m_inv = MatrixInverse(mA, NULL) ;
  if   (!m_inv)
    ErrorReturn(NULL, (ERROR_BADPARM, 
                       "MRIinverseLinearTransform: xform is singular!")) ;
  mri_dst = MRIlinearTransform(mri_src, mri_dst, m_inv) ;
  MatrixFree(&m_inv) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Convert a transform from RAS to voxel coordinates, then apply
  it to an MRI.
------------------------------------------------------*/
MRI *
MRIapplyRASlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *m_ras_xform)
{
  MATRIX   *m_voxel_xform ;

  m_voxel_xform = MRIrasXformToVoxelXform(mri_src, mri_dst, m_ras_xform, NULL);
  mri_dst = MRIlinearTransform(mri_src, mri_dst, m_voxel_xform) ;
  MatrixFree(&m_voxel_xform) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Convert a transform from RAS to voxel coordinates, then apply
  it to an MRI.
------------------------------------------------------*/
MRI *
MRIapplyRASinverseLinearTransform(MRI *mri_src, MRI *mri_dst, 
                                  MATRIX *m_ras_xform)
{
  MATRIX   *m_voxel_xform ;

  m_voxel_xform = MRIrasXformToVoxelXform(mri_src, mri_dst, m_ras_xform, NULL);
  mri_dst = MRIinverseLinearTransform(mri_src, mri_dst, m_voxel_xform) ;
  MatrixFree(&m_voxel_xform) ;
  return(mri_dst) ;
}

/*-----------------------------------------------------
  Perform an linear coordinate transformation x' = Ax on
  the MRI image mri_src into mri_dst using sinc interp.
------------------------------------------------------*/
MRI *
MRIsincTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA, int hw)
{
  int    y1, y2, y3, width, height, depth ;
  VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  MATRIX *mAinv ;     /* inverse of mA */
  Real   val, x1, x2, x3 ;

  mAinv = MatrixInverse(mA, NULL) ;      /* will sample from dst back to src */
  if (!mAinv)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRIsincTransform: xform is singular")) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  v_Y->rptr[4][1] = 1.0f ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    V3_Z(v_Y) = y3 ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        V3_X(v_Y) = y1 ;
        MatrixMultiply(mAinv, v_Y, v_X) ;
        
        x1 = V3_X(v_X) ; x2 = V3_Y(v_X) ; x3 = V3_Z(v_X) ;

        if (nint(y1) == 13 && nint(y2) == 10 && nint(y3) == 7)
          DiagBreak() ;
        if (nint(x1) == 13 && nint(x2) == 10 && nint(x3) == 7)
        {
#if 0
          fprintf(stderr, "(%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
                  (float)x1, (float)x2, (float)x3, 
                  (float)y1, (float)y2, (float)y3) ;
#endif
          DiagBreak() ;
        }

        if (x1 > -1 && x1 < width &&
            x2 > -1 && x2 < height &&
            x3 > -1 && x3 < depth)
        {
          MRIsincSampleVolume(mri_src, x1, x2, x3, hw, &val);
          MRIvox(mri_dst,y1,y2,y3) = (BUFTYPE)nint(val) ;
        }
      }
    }
  }

  MatrixFree(&v_X) ;
  MatrixFree(&mAinv) ;
  MatrixFree(&v_Y) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------------------
  MRIlinearTransform() - for historical reasons, this uses trilinear
  interpolation. This the operations under this function name can
  now (2/20/02) be found under MRIlinearTransformInterp().
  -----------------------------------------------------------------*/
MRI *
MRIlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA)
{
  MRIlinearTransformInterp(mri_src, mri_dst, mA, SAMPLE_TRILINEAR);
  return(mri_dst);
}
/*-------------------------------------------------------------------
  MRIlinearTransformInterp() Perform linear coordinate transformation 
  x' = Ax on the MRI image mri_src into mri_dst using the specified
  interpolation method. A is a voxel-to-voxel transform.
  ------------------------------------------------------------------*/
MRI *
MRIlinearTransformInterp(MRI *mri_src, MRI *mri_dst, MATRIX *mA,
       int InterpMethod)
{
  int    y1, y2, y3, width, height, depth ;
  VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  MATRIX *mAinv ;     /* inverse of mA */
  Real   val, x1, x2, x3 ;

  if(InterpMethod != SAMPLE_NEAREST &&
     InterpMethod != SAMPLE_TRILINEAR &&
     InterpMethod != SAMPLE_SINC ){
    printf("ERROR: MRIlinearTransformInterp: unrecoginzed interpolation "
     "method %d\n",InterpMethod);
  }

  mAinv = MatrixInverse(mA, NULL) ;      /* will sample from dst back to src */
  if (!mAinv)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRIlinearTransform: xform is singular")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  width = mri_dst->width ; height = mri_dst->height ; depth = mri_dst->depth ;
  v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  v_Y->rptr[4][1] = 1.0f ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    V3_Z(v_Y) = y3 ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        V3_X(v_Y) = y1 ;
        MatrixMultiply(mAinv, v_Y, v_X) ;
        
        x1 = V3_X(v_X) ; x2 = V3_Y(v_X) ; x3 = V3_Z(v_X) ;

        if (nint(y1) == 13 && nint(y2) == 10 && nint(y3) == 7)
          DiagBreak() ;
        if (nint(x1) == 13 && nint(x2) == 10 && nint(x3) == 7)
        {
#if 0
          fprintf(stderr, "(%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
                  (float)x1, (float)x2, (float)x3, 
                  (float)y1, (float)y2, (float)y3) ;
#endif
          DiagBreak() ;
        }

        //MRIsampleVolume(mri_src, x1, x2, x3, &val);
        MRIsampleVolumeType(mri_src, x1, x2, x3, &val, InterpMethod);

        switch (mri_dst->type)
        {
        case MRI_UCHAR:
          MRIvox(mri_dst,y1,y2,y3) = (BUFTYPE)nint(val) ;
          break ;
        case MRI_SHORT:
          MRISvox(mri_dst,y1,y2,y3) = (short)nint(val) ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri_dst,y1,y2,y3) = (float)(val) ;
          break ;
        case MRI_INT:
          MRIIvox(mri_dst,y1,y2,y3) = nint(val) ;
          break ;
        default:
          ErrorReturn(NULL, 
                      (ERROR_UNSUPPORTED, 
                       "MRIlinearTransform: unsupported dst type %d",
                       mri_dst->type)) ;
          break ;
        }
      }
    }
  }

  MatrixFree(&v_X) ;
  MatrixFree(&mAinv) ;
  MatrixFree(&v_Y) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
MRI *
MRIconcatenateFrames(MRI *mri_frame1, MRI *mri_frame2, MRI *mri_dst)
{
  int       width, height, depth, x, y, z ;
  BUFTYPE   *pf1, *pf2, *pdst1, *pdst2 ;

  if (mri_frame1->type != MRI_UCHAR || mri_frame1->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,"MRIconcatenateFrames: src must be UCHAR"));

  width = mri_frame1->width ; 
  height = mri_frame1->height ;
  depth = mri_frame1->depth ;
  if (mri_dst == NULL)
  {
    mri_dst = MRIallocSequence(width, height, depth, mri_frame1->type, 2) ;
    MRIcopyHeader(mri_frame1, mri_dst) ;
  }
  if (!mri_dst)
    ErrorExit(ERROR_NOMEMORY, "MRIconcatenateFrames: could not alloc dst") ;


  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst1 = &MRIvox(mri_dst, 0, y, z) ;
      pdst2 = &MRIseq_vox(mri_dst, 0, y, z, 1) ;
      pf1 = &MRIvox(mri_frame1, 0, y, z) ;
      pf2 = &MRIvox(mri_frame2, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        *pdst1++ = *pf1++ ;
        *pdst2++ = *pf2++ ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *
MRIcopyFrame(MRI *mri_src, MRI *mri_dst, int src_frame, int dst_frame)
{
  int       width, height, depth, y, z, bytes ;
  BUFTYPE   *psrc, *pdst ;

  width = mri_src->width ; 
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = 
      MRIallocSequence(width, height, depth, mri_src->type, dst_frame+1) ;
  if (!mri_dst)
    ErrorExit(ERROR_NOMEMORY, "MRIcopyFrame: could not alloc dst") ;

  if (mri_src->type != mri_dst->type)
    ErrorReturn(NULL,(ERROR_UNSUPPORTED,
                      "MRIcopyFrame: src and dst must be same type"));


  switch (mri_src->type)
  {
  case MRI_UCHAR:
    bytes = sizeof(unsigned char) ;
    break ;
  case MRI_SHORT:
    bytes = sizeof(short) ;
    break ;
  case MRI_INT:
    bytes = sizeof(int) ;
    break ;
  case MRI_FLOAT:
    bytes = sizeof(float) ;
    break ;
  default:
    ErrorReturn(NULL, (ERROR_BADPARM, 
                       "MRIcopyFrame: unsupported src format %d",
                       mri_src->type));
    break ;
  }
  bytes *= width ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      switch (mri_src->type)
      {
      default:  /* already handled above */
      case MRI_UCHAR:
        psrc = &MRIseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = &MRIseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      case MRI_SHORT:
        psrc = (BUFTYPE *)&MRISseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = (BUFTYPE *)&MRISseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      case MRI_INT:
        psrc = (BUFTYPE *)&MRIIseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = (BUFTYPE *)&MRIIseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      case MRI_FLOAT:
        psrc = (BUFTYPE *)&MRIFseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = (BUFTYPE *)&MRIFseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      }
      memmove(pdst, psrc, bytes) ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
double
MRImeanFrame(MRI *mri, int frame)
{
  int       width, height, depth, x, y, z ;
  double    mean ;
  BUFTYPE   *psrc ;

  width = mri->width ; 
  height = mri->height ;
  depth = mri->depth ;
  if (mri->type != MRI_UCHAR)
    ErrorReturn(0.0,(ERROR_UNSUPPORTED,"MRImeanFrame: src must be UCHAR"));
  if (mri->nframes <= frame)
    ErrorReturn(0.0,(ERROR_BADPARM,
                     "MRImeanFrame: frame %d out of bounds (%d)",
                     frame, mri->nframes));

  for (mean = 0.0, z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIseq_vox(mri, 0, y, z, frame) ;
      for (x = 0 ; x < width ; x++)
      {
        mean += (double)*psrc++ ;
      }
    }
  }
  mean /= (double)(width*height*depth) ;
  return(mean) ;
}

#define UCHAR_MIN  0.0
#define UCHAR_MAX  255.0
#define SHORT_MIN  -32768.0
#define SHORT_MAX  32767.0
#define INT_MIN    -2147483648.0
#define INT_MAX    2147483647.0
#define LONG_MIN   -2147483648.0
#define LONG_MAX   2147483647.0

#define N_HIST_BINS 1000

/*--------------------------------------------------------------
  MRISeqchangeType() - changes the data type for a 3D or 4D volume.
  This simply changes the volume dimensions so that it appears to be a
  3D volume, then calls MRIchangeType(), and then resets the
  dimensions to their original values. The values of the volume can be
  rescaled between f_low and f_high.
  ------------------------------------------------------------*/
MRI *MRISeqchangeType(MRI *vol, int dest_type, float f_low, 
          float f_high, int no_scale_option_flag)
{
  int nslices, nframes;
  MRI *mri;

  /* Change vol dimensions to make it look like a single frame */
  nslices = vol->depth;
  nframes = vol->nframes;
  vol->depth = nslices*nframes;
  vol->nframes = 1;

  /* Change the type */
  mri = MRIchangeType(vol,dest_type,f_low,f_high,no_scale_option_flag);

  /* Change vol dimensions back to original */
  vol->depth = nslices;
  vol->nframes = nframes;

  /* Check for error */
  if(mri == NULL) {
    fprintf(stderr,"ERROR: MRISeqchangeType: MRIchangeType\n");
    return(NULL);
  }

  /* Change mri dimensions back to original */
  mri->depth = nslices;
  mri->nframes = nframes;

  return(mri);
}
/*-----------------------------------------------------------
  MRIchangeType() - changes the data type of a 3D MRI volume,
  with optional rescaling. Use MRISeqchangeType() for 3D or
  4D volumes.
  ---------------------------------------------------------*/
MRI *MRIchangeType(MRI *src, int dest_type, float f_low, 
       float f_high, int no_scale_option_flag)
{

  MRI *dest = NULL;
  int i, j, k;
  float val;
  int no_scale_flag = FALSE;
  float scale, dest_min, dest_max; /* new = scale * (val - min) */
  float src_min, src_max;
  int hist_bins[N_HIST_BINS];
  float bin_size;
  int bin;
  int nth, n_passed;

  /* ----- shut the compiler up ----- */
  val = 0.0;
  dest_min = dest_max = 0.0;

  if(src->type == dest_type)
  {
    dest = MRIcopy(src, NULL);
    return(dest);
  }

  MRIlimits(src, &src_min, &src_max);

  if(src->type == MRI_UCHAR && (dest_type == MRI_SHORT || dest_type == MRI_INT || dest_type == MRI_LONG || dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else if(src->type == MRI_SHORT && (dest_type == MRI_INT || dest_type == MRI_LONG || dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else if(src->type == MRI_LONG && (dest_type == MRI_INT || dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else if(src->type == MRI_INT && (dest_type == MRI_LONG || dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else
  {
    if(no_scale_option_flag)
    {
      if(dest_type == MRI_UCHAR && src_min >= UCHAR_MIN && src_max <= UCHAR_MAX)
        no_scale_flag = TRUE;
      if(dest_type == MRI_SHORT && src_min >= SHORT_MIN && src_max <= SHORT_MAX)
        no_scale_flag = TRUE;
      if(dest_type == MRI_INT && src_min >= INT_MIN && src_max <= INT_MAX)
        no_scale_flag = TRUE;
      if(dest_type == MRI_LONG && src_min >= LONG_MIN && src_max <= LONG_MAX)
        no_scale_flag = TRUE;
    }
  }

  if(no_scale_flag)
  {

    dest = MRIalloc(src->width, src->height, src->depth, dest_type);
    MRIcopyHeader(src, dest);
    dest->type = dest_type;

    for(i = 0;i < src->width;i++)
      for(j = 0;j < src->height;j++)
        for(k = 0;k < src->depth;k++)
        {

          if(src->type == MRI_UCHAR)
            val = (float)MRIvox(src, i, j, k);
          if(src->type == MRI_SHORT)
            val = (float)MRISvox(src, i, j, k);
          if(src->type == MRI_INT)
            val = (float)MRIIvox(src, i, j, k);
          if(src->type == MRI_LONG)
            val = (float)MRILvox(src, i, j, k);
          if(src->type == MRI_FLOAT)
            val = (float)MRIFvox(src, i, j, k);

          if(dest_type == MRI_UCHAR)
            MRIvox(dest, i, j, k) = (unsigned char)val;
          if(dest_type == MRI_SHORT)
            MRISvox(dest, i, j, k) = (short)val;
          if(dest_type == MRI_INT)
            MRIIvox(dest, i, j, k) = (int)val;
          if(dest_type == MRI_LONG)
            MRILvox(dest, i, j, k) = (long)val;
          if(dest_type == MRI_FLOAT)
            MRIFvox(dest, i, j, k) = (float)val;

        }

  }
  else
  {

    /* ----- build a histogram ----- */

    bin_size = (src_max - src_min) / (float)N_HIST_BINS;

    for(i = 0;i < N_HIST_BINS;i++)
      hist_bins[i] = 0;

    for(i = 0;i < src->width;i++)
      for(j = 0;j < src->height;j++)
        for(k = 0;k < src->depth;k++)
        {

          if(src->type == MRI_UCHAR)
            val = (float)MRIvox(src, i, j, k);
          if(src->type == MRI_SHORT)
            val = (float)MRISvox(src, i, j, k);
          if(src->type == MRI_INT)
            val = (float)MRIIvox(src, i, j, k);
          if(src->type == MRI_LONG)
            val = (float)MRILvox(src, i, j, k);
          if(src->type == MRI_FLOAT)
            val = (float)MRIFvox(src, i, j, k);

          bin = (int)((val - src_min) / bin_size);

          if(bin < 0)
            bin = 0;
          if(bin >= N_HIST_BINS)
            bin = N_HIST_BINS-1;

          hist_bins[bin]++;

        }

    nth = (int)(f_low * src->width * src->height * src->depth);
    for(n_passed = 0,bin = 0;n_passed < nth && bin < N_HIST_BINS;bin++)
      n_passed += hist_bins[bin];
    src_min = (float)bin * bin_size + src_min;

    nth = (int)((1.0-f_high) * src->width * src->height * src->depth);
    for(n_passed = 0,bin = N_HIST_BINS-1;n_passed < nth && bin > 0;bin--)
      n_passed += hist_bins[bin];    
    src_max = (float)bin * bin_size + src_min;

    if(src_min >= src_max)
    {
      ErrorReturn(NULL, (ERROR_BADPARM, "MRIchangeType(): after hist: src_min = %g, src_max = %g (f_low = %g, f_high = %g)", src_min, src_max, f_low, f_high));
    }

    /* ----- continue ----- */

    if(dest_type == MRI_UCHAR)
    {
      dest_min = UCHAR_MIN;
      dest_max = UCHAR_MAX;
    }
    if(dest_type == MRI_SHORT)
    {
      dest_min = SHORT_MIN;
      dest_max = SHORT_MAX;
    }
    if(dest_type == MRI_INT)
    {
      dest_min = INT_MIN;
      dest_max = INT_MAX;
    }
    if(dest_type == MRI_LONG)
    {
      dest_min = LONG_MIN;
      dest_max = LONG_MAX;
    }

    scale = (dest_max - dest_min) / (src_max - src_min);

    dest = MRIalloc(src->width, src->height, src->depth, dest_type);
    MRIcopyHeader(src, dest);
    dest->type = dest_type;

    for(i = 0;i < src->width;i++)
      for(j = 0;j < src->height;j++)
        for(k = 0;k < src->depth;k++)
        {

          if(src->type == MRI_SHORT)
            val = MRISvox(src, i, j, k);
          if(src->type == MRI_INT)
            val = MRIIvox(src, i, j, k);
          if(src->type == MRI_LONG)
            val = MRILvox(src, i, j, k);
          if(src->type == MRI_FLOAT)
            val = MRIFvox(src, i, j, k);

          val = dest_min + scale * (val - src_min);

          if(dest->type == MRI_UCHAR)
          {
            if(val < UCHAR_MIN)
              val = UCHAR_MIN;
            if(val > UCHAR_MAX)
              val = UCHAR_MAX;
            MRIvox(dest, i, j, k) = (unsigned char)val;
          }
          if(dest->type == MRI_SHORT)
          {
            if(val < SHORT_MIN)
              val = SHORT_MIN;
            if(val > SHORT_MAX)
              val = SHORT_MAX;
            MRISvox(dest, i, j, k) = (short)val;
          }
          if(dest->type == MRI_INT)
          {
            if(val < INT_MIN)
              val = INT_MIN;
            if(val > INT_MAX)
              val = INT_MAX;
            MRIIvox(dest, i, j, k) = (int)val;
          }
          if(dest->type == MRI_LONG)
          {
            if(val < LONG_MIN)
              val = LONG_MIN;
            if(val > LONG_MAX)
              val = LONG_MAX;
            MRILvox(dest, i, j, k) = (long)val;
          }

        }

  }

  return(dest);

} /* end MRIchangeType() */

/*-----------------------------------------------------*/
MATRIX *MRIgetResampleMatrix(MRI *src, MRI *template_vol)
{

  MATRIX *src_mat, *dest_mat; /* from i to ras */
  float src_center_x, src_center_y, src_center_z;
  float dest_center_x, dest_center_y, dest_center_z;
  float s14, s24, s34;
  float d14, d24, d34;
  float src_det, dest_det;
  MATRIX *src_inv, *m;

  /* ----- fake the ras values if ras_good_flag is not set ----- */
  if(!src->ras_good_flag)
  {
    if(src->slice_direction == MRI_CORONAL)
    {
      src->x_r = -1; src->x_a =  0; src->x_s =  0;
      src->y_r =  0; src->y_a =  0; src->y_s = -1;
      src->z_r =  0; src->z_a =  1; src->z_s =  0;
    }
    else if(src->slice_direction == MRI_SAGITTAL)
    {
      src->x_r =  0; src->x_a =  1; src->x_s =  0;
      src->y_r =  0; src->y_a =  0; src->y_s = -1;
      src->z_r =  1; src->z_a =  0; src->z_s =  0;
    }
    else if(src->slice_direction == MRI_HORIZONTAL)
    {
      src->x_r = -1; src->x_a =  0; src->x_s =  0;
      src->y_r =  0; src->y_a = -1; src->y_s =  0;
      src->z_r =  0; src->z_a =  0; src->z_s = -1;
    }
    else
    {
      ErrorReturn(NULL, (ERROR_BADPARM, "MRIresample(): source volume orientation is unknown"));
    }
  }

  src_mat = MatrixAlloc(4, 4, MATRIX_REAL);
  if(src_mat == NULL)
  {
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIresample(): can't allocate source volume matrix"));
  }

  dest_mat = MatrixAlloc(4, 4, MATRIX_REAL);
  if(dest_mat == NULL)
  {
    MatrixFree(&src_mat);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIresample(): can't allocate destination volume matrix"));
  }

/*

below: solve each row of src_mat * [midx;midy;midz;1] = centerr, centera, centers
and for dest

  src->x_r * src->xsize * src_center_x + src->y_r * src->ysize * src_center_y + src->z_r * src->zsize * src_center_z + s14 = src->c_r;
  src->x_a * src->xsize * src_center_x + src->y_a * src->ysize * src_center_y + src->z_a * src->zsize * src_center_z + s24 = src->c_a;
  src->x_s * src->xsize * src_center_x + src->y_s * src->ysize * src_center_y + src->z_s * src->zsize * src_center_z + s34 = src->c_s;

  dest->x_r * dest->xsize * dest_center_x + dest->y_r * dest->ysize * dest_center_y + dest->z_r * dest->zsize * dest_center_z + d14 = dest->c_r;
  dest->x_a * dest->xsize * dest_center_x + dest->y_a * dest->ysize * dest_center_y + dest->z_a * dest->zsize * dest_center_z + d24 = dest->c_a;
  dest->x_s * dest->xsize * dest_center_x + dest->y_s * dest->ysize * dest_center_y + dest->z_s * dest->zsize * dest_center_z + d34 = dest->c_s;

*/

  src_center_x = (float)(src->width - 1) / 2.0;
  src_center_y = (float)(src->height - 1) / 2.0;
  src_center_z = (float)(src->depth - 1) / 2.0;

  dest_center_x = (float)(template_vol->width - 1) / 2.0;
  dest_center_y = (float)(template_vol->height - 1) / 2.0;
  dest_center_z = (float)(template_vol->depth - 1) / 2.0;

  s14 = src->c_r - (src->x_r * src->xsize * src_center_x + src->y_r * src->ysize * src_center_y + src->z_r * src->zsize * src_center_z);
  s24 = src->c_a - (src->x_a * src->xsize * src_center_x + src->y_a * src->ysize * src_center_y + src->z_a * src->zsize * src_center_z);
  s34 = src->c_s - (src->x_s * src->xsize * src_center_x + src->y_s * src->ysize * src_center_y + src->z_s * src->zsize * src_center_z);

  d14 = template_vol->c_r - (template_vol->x_r * template_vol->xsize * dest_center_x + template_vol->y_r * template_vol->ysize * dest_center_y + template_vol->z_r * template_vol->zsize * dest_center_z);
  d24 = template_vol->c_a - (template_vol->x_a * template_vol->xsize * dest_center_x + template_vol->y_a * template_vol->ysize * dest_center_y + template_vol->z_a * template_vol->zsize * dest_center_z);
  d34 = template_vol->c_s - (template_vol->x_s * template_vol->xsize * dest_center_x + template_vol->y_s * template_vol->ysize * dest_center_y + template_vol->z_s * template_vol->zsize * dest_center_z);

  stuff_four_by_four(src_mat, src->x_r * src->xsize, src->y_r * src->ysize, src->z_r * src->zsize, s14, 
                              src->x_a * src->xsize, src->y_a * src->ysize, src->z_a * src->zsize, s24, 
                              src->x_s * src->xsize, src->y_s * src->ysize, src->z_s * src->zsize, s34, 
                                                0.0,                   0.0,                   0.0, 1.0);

  stuff_four_by_four(dest_mat, template_vol->x_r * template_vol->xsize, template_vol->y_r * template_vol->ysize, template_vol->z_r * template_vol->zsize, d14, 
                               template_vol->x_a * template_vol->xsize, template_vol->y_a * template_vol->ysize, template_vol->z_a * template_vol->zsize, d24, 
                               template_vol->x_s * template_vol->xsize, template_vol->y_s * template_vol->ysize, template_vol->z_s * template_vol->zsize, d34, 
                                                                   0.0,                                     0.0,                                     0.0, 1.0);

  src_det = MatrixDeterminant(src_mat);
  dest_det = MatrixDeterminant(dest_mat);

  if(src_det == 0.0)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "MRIresample(): source matrix has zero determinant; matrix is:");
    MatrixPrint(stderr, src_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  if(dest_det == 0.0)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "MRIresample(): destination matrix has zero determinant; matrix is:");
    MatrixPrint(stderr, dest_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  src_inv = MatrixInverse(src_mat, NULL);

  if(src_inv == NULL)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "MRIresample(): error inverting matrix; determinant is %g, matrix is:", src_det);
    MatrixPrint(stderr, src_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  m = MatrixMultiply(src_inv, dest_mat, NULL);
  if(m == NULL)
    return(NULL);

  MatrixFree(&src_inv);
  MatrixFree(&src_mat);
  MatrixFree(&dest_mat);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("MRIresample() matrix is:\n");
    MatrixPrint(stdout, m);
  }
  return(m) ;

} /* end MRIreslice() */

MRI *MRIresample(MRI *src, MRI *template_vol, int resample_type)
{

  MRI *dest = NULL;
  MATRIX *src_mat, *dest_mat; /* from i to ras */
  float src_center_x, src_center_y, src_center_z;
  float dest_center_x, dest_center_y, dest_center_z;
  float s14, s24, s34;
  float d14, d24, d34;
  float src_det, dest_det;
  MATRIX *src_inv, *m;
  int di, dj, dk;
  int si, sj, sk;
  float si_f, sj_f, sk_f;
  float si_ff, sj_ff, sk_ff;
  MATRIX *sp, *dp;
  float val, val000, val001, val010, val011, val100, val101, val110, val111;
  float w000, w001, w010, w011, w100, w101, w110, w111;
  float si_f2, sj_f2, sk_f2;
  float ii2, ij2, ik2, isi_f2, isj_f2, isk_f2;
  float w[8];
  int wi[8];
  int mwi;
  Real pval;

  /* ----- keep the compiler quiet ----- */
  val = 0.0;
  val000 = val001 = val010 = val011 = val100 = val101 = val110 = val111 = 0.0;

#if 0
  if(src->type != template_vol->type)
  {
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIresample(): source and destination types must be identical"));
  }
#endif

  /* ----- fake the ras values if ras_good_flag is not set ----- */
  if(!src->ras_good_flag)
  {
    if(src->slice_direction == MRI_CORONAL)
    {
      src->x_r = -1; src->x_a =  0; src->x_s =  0;
      src->y_r =  0; src->y_a =  0; src->y_s = -1;
      src->z_r =  0; src->z_a =  1; src->z_s =  0;
    }
    else if(src->slice_direction == MRI_SAGITTAL)
    {
      src->x_r =  0; src->x_a =  1; src->x_s =  0;
      src->y_r =  0; src->y_a =  0; src->y_s = -1;
      src->z_r =  1; src->z_a =  0; src->z_s =  0;
    }
    else if(src->slice_direction == MRI_HORIZONTAL)
    {
      src->x_r = -1; src->x_a =  0; src->x_s =  0;
      src->y_r =  0; src->y_a = -1; src->y_s =  0;
      src->z_r =  0; src->z_a =  0; src->z_s = -1;
    }
    else
    {
      ErrorReturn(NULL, (ERROR_BADPARM, "MRIresample(): source volume orientation is unknown"));
    }
  }

  src_mat = MatrixAlloc(4, 4, MATRIX_REAL);
  if(src_mat == NULL)
  {
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIresample(): can't allocate source volume matrix"));
  }

  dest_mat = MatrixAlloc(4, 4, MATRIX_REAL);
  if(dest_mat == NULL)
  {
    MatrixFree(&src_mat);
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIresample(): can't allocate destination volume matrix"));
  }

/*

below: solve each row of src_mat * [midx;midy;midz;1] = centerr, centera, centers
and for dest

  src->x_r * src->xsize * src_center_x + src->y_r * src->ysize * src_center_y + src->z_r * src->zsize * src_center_z + s14 = src->c_r;
  src->x_a * src->xsize * src_center_x + src->y_a * src->ysize * src_center_y + src->z_a * src->zsize * src_center_z + s24 = src->c_a;
  src->x_s * src->xsize * src_center_x + src->y_s * src->ysize * src_center_y + src->z_s * src->zsize * src_center_z + s34 = src->c_s;

  dest->x_r * dest->xsize * dest_center_x + dest->y_r * dest->ysize * dest_center_y + dest->z_r * dest->zsize * dest_center_z + d14 = dest->c_r;
  dest->x_a * dest->xsize * dest_center_x + dest->y_a * dest->ysize * dest_center_y + dest->z_a * dest->zsize * dest_center_z + d24 = dest->c_a;
  dest->x_s * dest->xsize * dest_center_x + dest->y_s * dest->ysize * dest_center_y + dest->z_s * dest->zsize * dest_center_z + d34 = dest->c_s;

*/

  src_center_x = (float)(src->width - 1) / 2.0;
  src_center_y = (float)(src->height - 1) / 2.0;
  src_center_z = (float)(src->depth - 1) / 2.0;

  dest_center_x = (float)(template_vol->width - 1) / 2.0;
  dest_center_y = (float)(template_vol->height - 1) / 2.0;
  dest_center_z = (float)(template_vol->depth - 1) / 2.0;

  s14 = src->c_r - (src->x_r * src->xsize * src_center_x + src->y_r * src->ysize * src_center_y + src->z_r * src->zsize * src_center_z);
  s24 = src->c_a - (src->x_a * src->xsize * src_center_x + src->y_a * src->ysize * src_center_y + src->z_a * src->zsize * src_center_z);
  s34 = src->c_s - (src->x_s * src->xsize * src_center_x + src->y_s * src->ysize * src_center_y + src->z_s * src->zsize * src_center_z);

  d14 = template_vol->c_r - (template_vol->x_r * template_vol->xsize * dest_center_x + template_vol->y_r * template_vol->ysize * dest_center_y + template_vol->z_r * template_vol->zsize * dest_center_z);
  d24 = template_vol->c_a - (template_vol->x_a * template_vol->xsize * dest_center_x + template_vol->y_a * template_vol->ysize * dest_center_y + template_vol->z_a * template_vol->zsize * dest_center_z);
  d34 = template_vol->c_s - (template_vol->x_s * template_vol->xsize * dest_center_x + template_vol->y_s * template_vol->ysize * dest_center_y + template_vol->z_s * template_vol->zsize * dest_center_z);

  stuff_four_by_four(src_mat, src->x_r * src->xsize, src->y_r * src->ysize, src->z_r * src->zsize, s14, 
                              src->x_a * src->xsize, src->y_a * src->ysize, src->z_a * src->zsize, s24, 
                              src->x_s * src->xsize, src->y_s * src->ysize, src->z_s * src->zsize, s34, 
                                                0.0,                   0.0,                   0.0, 1.0);

  stuff_four_by_four(dest_mat, template_vol->x_r * template_vol->xsize, template_vol->y_r * template_vol->ysize, template_vol->z_r * template_vol->zsize, d14, 
                               template_vol->x_a * template_vol->xsize, template_vol->y_a * template_vol->ysize, template_vol->z_a * template_vol->zsize, d24, 
                               template_vol->x_s * template_vol->xsize, template_vol->y_s * template_vol->ysize, template_vol->z_s * template_vol->zsize, d34, 
                                                                   0.0,                                     0.0,                                     0.0, 1.0);

  src_det = MatrixDeterminant(src_mat);
  dest_det = MatrixDeterminant(dest_mat);

  if(src_det == 0.0)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "MRIresample(): source matrix has zero determinant; matrix is:");
    MatrixPrint(stderr, src_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  if(dest_det == 0.0)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "MRIresample(): destination matrix has zero determinant; matrix is:");
    MatrixPrint(stderr, dest_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  src_inv = MatrixInverse(src_mat, NULL);

  if(src_inv == NULL)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "MRIresample(): error inverting matrix; determinant is %g, matrix is:", src_det);
    MatrixPrint(stderr, src_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  m = MatrixMultiply(src_inv, dest_mat, NULL);
  if(m == NULL)
  {
    MatrixFree(&src_inv);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("MRIresample() matrix is:\n");
    MatrixPrint(stdout, m);
  }

  dest = MRIalloc(template_vol->width, template_vol->height, template_vol->depth, src->type);
  if(dest == NULL)
    return(NULL);
  MRIcopyHeader(template_vol, dest);

  sp = MatrixAlloc(4, 1, MATRIX_REAL);
  dp = MatrixAlloc(4, 1, MATRIX_REAL);

  *MATRIX_RELT(dp, 4, 1) = 1.0;

  for(di = 0;di < template_vol->width;di++)
  {
    for(dj = 0;dj < template_vol->height;dj++)
    {
      for(dk = 0;dk < template_vol->depth;dk++)
      {

        *MATRIX_RELT(dp, 1, 1) = (float)di;
        *MATRIX_RELT(dp, 2, 1) = (float)dj;
        *MATRIX_RELT(dp, 3, 1) = (float)dk;

        MatrixMultiply(m, dp, sp);

        si_ff = *MATRIX_RELT(sp, 1, 1);
        sj_ff = *MATRIX_RELT(sp, 2, 1);
        sk_ff = *MATRIX_RELT(sp, 3, 1);

        si = (int)floor(si_ff);
        sj = (int)floor(sj_ff);
        sk = (int)floor(sk_ff);

        if (si == 147 && sj == 91 && sk == 86)
          DiagBreak() ;
        if (di == 129 && dj == 164 && dk == 147)
          DiagBreak() ;
        si_f = si_ff - si;
        sj_f = sj_ff - sj;
        sk_f = sk_ff - sk;

        if(resample_type == RESAMPLE_SINC)
        {
          MRIsincSampleVolume(src, si_ff, sj_ff, sk_ff, 5, &pval);
          val = (float)pval;
        }
        else if(si < 0.0 || si >= (float)(src->width - 2) ||
                sj < 0.0 || sj >= (float)(src->height - 2) ||
                sk < 0.0 || sk >= (float)(src->depth - 2))
        {
          val = 0.0;
        }
        else
        {

          if(src->type == MRI_UCHAR)
          {
            val000 = (float)MRIvox(src, si    , sj    , sk    );
            val001 = (float)MRIvox(src, si    , sj    , sk + 1);
            val010 = (float)MRIvox(src, si    , sj + 1, sk    );
            val011 = (float)MRIvox(src, si    , sj + 1, sk + 1);
            val100 = (float)MRIvox(src, si + 1, sj    , sk    );
            val101 = (float)MRIvox(src, si + 1, sj    , sk + 1);
            val110 = (float)MRIvox(src, si + 1, sj + 1, sk    );
            val111 = (float)MRIvox(src, si + 1, sj + 1, sk + 1);
          }

          if (si == 154 && sj == 134 && sk == 136)
            DiagBreak() ;
          if(src->type == MRI_SHORT)
          {
            val000 = (float)MRISvox(src, si    , sj    , sk    );
            val001 = (float)MRISvox(src, si    , sj    , sk + 1);
            val010 = (float)MRISvox(src, si    , sj + 1, sk    );
            val011 = (float)MRISvox(src, si    , sj + 1, sk + 1);
            val100 = (float)MRISvox(src, si + 1, sj    , sk    );
            val101 = (float)MRISvox(src, si + 1, sj    , sk + 1);
            val110 = (float)MRISvox(src, si + 1, sj + 1, sk    );
            val111 = (float)MRISvox(src, si + 1, sj + 1, sk + 1);
          }

          if(src->type == MRI_INT)
          {
            val000 = (float)MRIIvox(src, si    , sj    , sk    );
            val001 = (float)MRIIvox(src, si    , sj    , sk + 1);
            val010 = (float)MRIIvox(src, si    , sj + 1, sk    );
            val011 = (float)MRIIvox(src, si    , sj + 1, sk + 1);
            val100 = (float)MRIIvox(src, si + 1, sj    , sk    );
            val101 = (float)MRIIvox(src, si + 1, sj    , sk + 1);
            val110 = (float)MRIIvox(src, si + 1, sj + 1, sk    );
            val111 = (float)MRIIvox(src, si + 1, sj + 1, sk + 1);
          }

          if(src->type == MRI_LONG)
          {
            val000 = (float)MRILvox(src, si    , sj    , sk    );
            val001 = (float)MRILvox(src, si    , sj    , sk + 1);
            val010 = (float)MRILvox(src, si    , sj + 1, sk    );
            val011 = (float)MRILvox(src, si    , sj + 1, sk + 1);
            val100 = (float)MRILvox(src, si + 1, sj    , sk    );
            val101 = (float)MRILvox(src, si + 1, sj    , sk + 1);
            val110 = (float)MRILvox(src, si + 1, sj + 1, sk    );
            val111 = (float)MRILvox(src, si + 1, sj + 1, sk + 1);
          }

          if(src->type == MRI_FLOAT)
          {
            val000 = (float)MRIFvox(src, si    , sj    , sk    );
            val001 = (float)MRIFvox(src, si    , sj    , sk + 1);
            val010 = (float)MRIFvox(src, si    , sj + 1, sk    );
            val011 = (float)MRIFvox(src, si    , sj + 1, sk + 1);
            val100 = (float)MRIFvox(src, si + 1, sj    , sk    );
            val101 = (float)MRIFvox(src, si + 1, sj    , sk + 1);
            val110 = (float)MRIFvox(src, si + 1, sj + 1, sk    );
            val111 = (float)MRIFvox(src, si + 1, sj + 1, sk + 1);
          }

          if(resample_type == RESAMPLE_INTERPOLATE)
          {
            val = (1.0-si_f) * (1.0-sj_f) * (1.0-sk_f) * val000 + 
                  (1.0-si_f) * (1.0-sj_f) * (    sk_f) * val001 + 
                  (1.0-si_f) * (    sj_f) * (1.0-sk_f) * val010 + 
                  (1.0-si_f) * (    sj_f) * (    sk_f) * val011 + 
                  (    si_f) * (1.0-sj_f) * (1.0-sk_f) * val100 + 
                  (    si_f) * (1.0-sj_f) * (    sk_f) * val101 + 
                  (    si_f) * (    sj_f) * (1.0-sk_f) * val110 + 
                  (    si_f) * (    sj_f) * (    sk_f) * val111;
          }

          if(resample_type == RESAMPLE_NEAREST)
          {
            if(si_f < 0.5)
            {
              if(sj_f < 0.5)
              {
                if(sk_f < 0.5)
                  val = val000;
                else
                  val = val001;
              }
              else
              {
                if(sk_f < 0.5)
                  val = val010;
                else
                  val = val011;
              }
            }
            else
            {
              if(sj_f < 0.5)
              {
                if(sk_f < 0.5)
                  val = val100;
                else
                  val = val101;
              }
              else
              {
                if(sk_f < 0.5)
                  val = val110;
                else
                  val = val111;
              }
            }
          }

          if(resample_type == RESAMPLE_WEIGHTED)
          {
/* unfinished */
            si_f2 = si_f * si_f;
            sj_f2 = sj_f * sj_f;
            sk_f2 = sk_f * sk_f;

            ii2 = 1. / (1 - 2*si_f + si_f2);
            ij2 = 1. / (1 - 2*sj_f + sj_f2);
            ik2 = 1. / (1 - 2*sk_f + sk_f2);

            isi_f2 = 1 / si_f2;
            isj_f2 = 1 / sj_f2;
            isk_f2 = 1 / sk_f2;

            w000 =   ii2    + ij2    + ik2;
            w001 =   ii2    + ij2    + isk_f2;
            w010 =   ii2    + isj_f2 + ik2;
            w011 =   ii2    + isj_f2 + isk_f2;
            w100 =   isi_f2 + ij2    + ik2;
            w101 =   isi_f2 + ij2    + isk_f2;
            w110 =   isi_f2 + isj_f2 + ik2;
            w111 =   isi_f2 + isj_f2 + isk_f2;

            w[0] = w[1] = w[2] = w[3] = w[4] = w[5] = w[6] = w[7] = 0.0;

            wi[0] = 0;
            wi[1] = 1;
            wi[2] = 2;
            wi[3] = 3;
            wi[4] = 4;
            wi[5] = 5;
            wi[6] = 6;
            wi[7] = 7;

            if(val001 == val000)
              wi[1] = 0;

            if(val010 == val001)
              wi[2] = 1;
            if(val010 == val000)
              wi[2] = 0;

            if(val011 == val010)
              wi[3] = 2;
            if(val011 == val001)
              wi[3] = 1;
            if(val011 == val000)
              wi[3] = 0;

            if(val100 == val011)
              wi[4] = 3;
            if(val100 == val010)
              wi[4] = 2;
            if(val100 == val001)
              wi[4] = 1;
            if(val100 == val000)
              wi[4] = 0;

            if(val101 == val100)
              wi[5] = 4;
            if(val101 == val011)
              wi[5] = 3;
            if(val101 == val010)
              wi[5] = 2;
            if(val101 == val001)
              wi[5] = 1;
            if(val101 == val000)
              wi[5] = 0;

            if(val110 == val101)
              wi[6] = 5;
            if(val110 == val100)
              wi[6] = 4;
            if(val110 == val011)
              wi[6] = 3;
            if(val110 == val010)
              wi[6] = 2;
            if(val110 == val001)
              wi[6] = 1;
            if(val110 == val000)
              wi[6] = 0;

            if(val111 == val110)
              wi[7] = 6;
            if(val111 == val101)
              wi[7] = 5;
            if(val111 == val100)
              wi[7] = 4;
            if(val111 == val011)
              wi[7] = 3;
            if(val111 == val010)
              wi[7] = 2;
            if(val111 == val001)
              wi[7] = 1;
            if(val111 == val000)
              wi[7] = 0;

            w[wi[0]] += w000;
            w[wi[1]] += w001;
            w[wi[2]] += w010;
            w[wi[3]] += w011;
            w[wi[4]] += w100;
            w[wi[5]] += w101;
            w[wi[6]] += w110;
            w[wi[7]] += w111;

            mwi = 0;

            if(w[1] > w[mwi])
              mwi = 1;
            if(w[2] > w[mwi])
              mwi = 2;
            if(w[3] > w[mwi])
              mwi = 3;
            if(w[4] > w[mwi])
              mwi = 4;
            if(w[5] > w[mwi])
              mwi = 5;
            if(w[6] > w[mwi])
              mwi = 6;
            if(w[7] > w[mwi])
              mwi = 7;

            if(mwi == 0)
              val = val000;
            if(mwi == 1)
              val = val001;
            if(mwi == 2)
              val = val010;
            if(mwi == 3)
              val = val011;
            if(mwi == 4)
              val = val100;
            if(mwi == 5)
              val = val101;
            if(mwi == 6)
              val = val110;
            if(mwi == 7)
              val = val111;

          }

        }

        if(dest->type == MRI_UCHAR)
          MRIvox(dest, di, dj, dk) = (unsigned char)val;
        if(dest->type == MRI_SHORT)
          MRISvox(dest, di, dj, dk) = (short)val;
        if(dest->type == MRI_INT)
          MRIIvox(dest, di, dj, dk) = (int)val;
        if(dest->type == MRI_LONG)
          MRILvox(dest, di, dj, dk) = (long)val;
        if(dest->type == MRI_FLOAT)
          MRIFvox(dest, di, dj, dk) = (float)val;

      }
    }
  }

  {
    MATRIX  *m_old_voxel_to_ras, *m_voxel_to_ras,
            *m_old_ras_to_voxel, *v_ras, *v_vox, *m_new_ras_to_voxel, *v_vox2 ;

    m_old_voxel_to_ras = MRIgetVoxelToRasXform(src) ;

    m_voxel_to_ras = MatrixMultiply(m_old_voxel_to_ras, m, NULL) ;
    dest->x_r = *MATRIX_RELT(m_voxel_to_ras,1,1)/dest->xsize;
    dest->x_a = *MATRIX_RELT(m_voxel_to_ras,2,1)/dest->xsize;
    dest->x_s = *MATRIX_RELT(m_voxel_to_ras,3,1)/dest->xsize;

    dest->y_r = *MATRIX_RELT(m_voxel_to_ras,1,2)/dest->ysize;
    dest->y_a = *MATRIX_RELT(m_voxel_to_ras,2,2)/dest->ysize;
    dest->y_s = *MATRIX_RELT(m_voxel_to_ras,3,2)/dest->ysize;

    dest->z_r = *MATRIX_RELT(m_voxel_to_ras,1,3)/dest->zsize;
    dest->z_a = *MATRIX_RELT(m_voxel_to_ras,2,3)/dest->zsize;
    dest->z_s = *MATRIX_RELT(m_voxel_to_ras,3,3)/dest->zsize;

    /* compute the RAS coordinates of the center of the dest. image
       and put them in c_r, c_a, and c_s.
    */
    v_vox = VectorAlloc(4, MATRIX_REAL) ;

    /* voxel coords of center of dest image */
    VECTOR_ELT(v_vox,4) = 1.0 ;
    VECTOR_ELT(v_vox,1) = (dest->width-1)/2.0 ;
    VECTOR_ELT(v_vox,2) = (dest->height-1)/2.0 ;
    VECTOR_ELT(v_vox,3) = (dest->depth-1)/2.0 ;

    v_vox2 = MatrixMultiply(m, v_vox, NULL) ; /* voxel coords in source image */
    v_ras = MatrixMultiply(m_old_voxel_to_ras, v_vox2, NULL) ; /* ras cntr of dest */

    dest->c_r = VECTOR_ELT(v_ras, 1);
    dest->c_a = VECTOR_ELT(v_ras, 2);
    dest->c_s = VECTOR_ELT(v_ras, 3);
    dest->ras_good_flag = 1 ;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      m_new_ras_to_voxel = MRIgetRasToVoxelXform(dest) ;
      m_old_ras_to_voxel = MRIgetRasToVoxelXform(src) ;
      V3_X(v_ras) = V3_Y(v_ras) = V3_Z(v_ras) = 0.0 ;
      MatrixMultiply(m_old_ras_to_voxel, v_ras, v_vox) ;
      printf("old RAS (0,0,0) -> (%2.0f, %2.0f, %2.0f)\n",
             VECTOR_ELT(v_vox,1), VECTOR_ELT(v_vox,2), VECTOR_ELT(v_vox,3)) ;
      MatrixMultiply(m_new_ras_to_voxel, v_ras, v_vox) ;
      printf("new RAS (0,0,0) -> (%2.0f, %2.0f, %2.0f)\n",
             VECTOR_ELT(v_vox,1), VECTOR_ELT(v_vox,2), VECTOR_ELT(v_vox,3)) ;
      MatrixFree(&m_new_ras_to_voxel) ; MatrixFree(&m_old_ras_to_voxel) ;
    }


    MatrixFree(&v_vox) ; MatrixFree(&v_vox2) ; MatrixFree(&v_ras) ;
    MatrixFree(&m_voxel_to_ras) ; MatrixFree(&m_old_voxel_to_ras) ;
  }


  MatrixFree(&dp);
  MatrixFree(&sp);

  MatrixFree(&src_inv);
  MatrixFree(&m);
  MatrixFree(&src_mat);
  MatrixFree(&dest_mat);

  return(dest);

} /* end MRIreslice() */

int MRIlimits(MRI *mri, float *min, float *max)
{

  float val;
  int i, j, k;

  if(mri == NULL)
    return(NO_ERROR);

  if(mri->slices == NULL)
    return(NO_ERROR);

  if(mri->type == MRI_UCHAR)
  {
    *min = *max = (float)MRIvox(mri, 0, 0, 0);
    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
        {
          val = (float)MRIvox(mri, i, j, k);
          if(val < *min)
            *min = val;
          if(val > *max)
            *max = val;
        }

  }

  if(mri->type == MRI_SHORT)
  {

    *min = *max = (float)MRISvox(mri, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
        {
          val = (float)MRISvox(mri, i, j, k);
          if(val < *min)
            *min = val;
          if(val > *max)
            *max = val;
        }

  }

  if(mri->type == MRI_INT)
  {

    *min = *max = (float)MRILvox(mri, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
        {
          val = (float)MRILvox(mri, i, j, k);
          if(val < *min)
            *min = val;
          if(val > *max)
            *max = val;
        }

  }

  if(mri->type == MRI_LONG)
  {

    *min = *max = (float)MRILvox(mri, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
        {
          val = (float)MRILvox(mri, i, j, k);
          if(val < *min)
            *min = val;
          if(val > *max)
            *max = val;
        }

  }

  if(mri->type == MRI_FLOAT)
  {

    *min = *max = (float)MRIFvox(mri, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
        {
          val = (float)MRIFvox(mri, i, j, k);
          if(val < *min)
            *min = val;
          if(val > *max)
            *max = val;
        }

  }

  return(NO_ERROR);

} /* end MRIlimits() */

int MRIprintStats(MRI *mri, FILE *stream)
{

  float min, max, mean, std;
  int n;

  MRIstats(mri, &min, &max, &n, &mean, &std);

  fprintf(stream, "%d values\n", n);
  fprintf(stream, "min = %g\n", min);
  fprintf(stream, "max = %g\n", max);
  fprintf(stream, "mean = %g\n", mean);
  fprintf(stream, "std = %g\n", std);

  return(NO_ERROR);

} /* end MRIprintStats() */

int MRIstats(MRI *mri, float *min, float *max, int *n_voxels, float *mean, float *std)
{

  float val;
  float sum, sq_sum;
  int i, j, k, t;
  float n;

  if(mri == NULL)
    return(NO_ERROR);

  if(mri->slices == NULL)
    return(NO_ERROR);

  sum = sq_sum = 0.0;
  *n_voxels = mri->width * mri->height * mri->depth * mri->nframes;
  n = (float)(*n_voxels);

  if(mri->type == MRI_UCHAR)
  {

    *max = *min = MRIseq_vox(mri, 0, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
          for(t = 0;t < mri->nframes;t++)
          {

            val = (float)MRIseq_vox(mri, i, j, k, t);

            if(val < *min)
              *min = val;
            if(val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if(mri->type == MRI_SHORT)
  {

    *min = *max = (float)MRISseq_vox(mri, 0, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
          for(t = 0;t < mri->nframes;t++)
          {

            val = (float)MRISseq_vox(mri, i, j, k, t);

            if(val < *min)
              *min = val;
            if(val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if(mri->type == MRI_INT)
  {

    *min = *max = (float)MRILseq_vox(mri, 0, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
          for(t = 0;t < mri->nframes;t++)
          {

            val = (float)MRILseq_vox(mri, i, j, k, t);

            if(val < *min)
              *min = val;
            if(val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if(mri->type == MRI_LONG)
  {

    *min = *max = (float)MRILseq_vox(mri, 0, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
          for(t = 0;t < mri->nframes;t++)
          {

            val = (float)MRILseq_vox(mri, i, j, k, t);

            if(val < *min)
              *min = val;
            if(val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if(mri->type == MRI_FLOAT)
  {

    *min = *max = (float)MRIFseq_vox(mri, 0, 0, 0, 0);

    for(i = 0;i < mri->width;i++)
      for(j = 0;j < mri->height;j++)
        for(k = 0;k < mri->depth;k++)
          for(t = 0;t < mri->nframes;t++)
          {

            val = (float)MRIFseq_vox(mri, i, j, k, t);

            if(val < *min)
              *min = val;
            if(val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  *mean = sum / n;
  *std = (sq_sum - n * (*mean) * (*mean)) / (n - 1.0);

  *std = sqrt(*std);

  return(NO_ERROR);

} /* end MRIstats() */

float MRIvolumeDeterminant(MRI *mri)
{

  MATRIX *m;
  float det;

  m = MatrixAlloc(4, 4, MATRIX_REAL);
  if(m == NULL)
    return(0.0);

  stuff_four_by_four(m, mri->x_r, mri->y_r, mri->z_r, 0.0, 
                        mri->x_a, mri->y_a, mri->z_a, 0.0, 
                        mri->x_s, mri->y_s, mri->z_s, 0.0, 
                             0.0,      0.0,      0.0, 1.0);

  det = MatrixDeterminant(m);

  MatrixFree(&m);

  return(det);

} /* end MRIvolumeDeterminant() */

int stuff_four_by_four(MATRIX *m, float m11, float m12, float m13, float m14, 
                                  float m21, float m22, float m23, float m24, 
                                  float m31, float m32, float m33, float m34, 
                                  float m41, float m42, float m43, float m44)
{

  if(m == NULL)
  {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "stuff_four_by_four(): matrix is NULL"));
  }

  if(m->rows != 4 || m->cols != 4)
  {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "stuff_four_by_four(): matrix is not four-by-four"));
  }

  *MATRIX_RELT(m, 1, 1) = m11;  *MATRIX_RELT(m, 1, 2) = m12;  *MATRIX_RELT(m, 1, 3) = m13;  *MATRIX_RELT(m, 1, 4) = m14;
  *MATRIX_RELT(m, 2, 1) = m21;  *MATRIX_RELT(m, 2, 2) = m22;  *MATRIX_RELT(m, 2, 3) = m23;  *MATRIX_RELT(m, 2, 4) = m24;
  *MATRIX_RELT(m, 3, 1) = m31;  *MATRIX_RELT(m, 3, 2) = m32;  *MATRIX_RELT(m, 3, 3) = m33;  *MATRIX_RELT(m, 3, 4) = m34;
  *MATRIX_RELT(m, 4, 1) = m41;  *MATRIX_RELT(m, 4, 2) = m42;  *MATRIX_RELT(m, 4, 3) = m43;  *MATRIX_RELT(m, 4, 4) = m44;

  return(NO_ERROR);

} /* end stuff_four_by_four() */

int apply_i_to_r(MRI *mri, MATRIX *m)
{

  float x_r, x_a, x_s;
  float y_r, y_a, y_s;
  float z_r, z_a, z_s;
  float mag;
  MATRIX *origin, *c;

  x_r = *MATRIX_RELT(m, 1, 1);  x_a = *MATRIX_RELT(m, 2, 1);  x_s = *MATRIX_RELT(m, 3, 1);
  mag = sqrt(x_r*x_r + x_a*x_a + x_s*x_s);
  mri->x_r = x_r / mag;  mri->x_a = x_a / mag;  mri->x_s = x_s / mag;
  mri->xsize = mag;

  y_r = *MATRIX_RELT(m, 1, 2);  y_a = *MATRIX_RELT(m, 2, 2);  y_s = *MATRIX_RELT(m, 3, 2);
  mag = sqrt(y_r*y_r + y_a*y_a + y_s*y_s);
  mri->y_r = y_r / mag;  mri->y_a = y_a / mag;  mri->y_s = y_s / mag;
  mri->ysize = mag;

  z_r = *MATRIX_RELT(m, 1, 3);  z_a = *MATRIX_RELT(m, 2, 3);  z_s = *MATRIX_RELT(m, 3, 3);
  mag = sqrt(z_r*z_r + z_a*z_a + z_s*z_s);
  mri->z_r = z_r / mag;  mri->z_a = z_a / mag;  mri->z_s = z_s / mag;
  mri->zsize = mag;

  origin = MatrixAlloc(4, 1, MATRIX_REAL);
  if(origin == NULL)
  {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "apply_i_to_r(): error allocating matrix"));
  }
  *MATRIX_RELT(origin, 1, 1) = (mri->width  - 1.0) / 2.0;
  *MATRIX_RELT(origin, 2, 1) = (mri->height - 1.0) / 2.0;
  *MATRIX_RELT(origin, 3, 1) = (mri->depth  - 1.0) / 2.0;
  *MATRIX_RELT(origin, 4, 1) = 1.0;

  c = MatrixMultiply(m, origin, NULL);
  if(c == NULL)
  {
    MatrixFree(&origin);
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "apply_i_to_r(): error multiplying matrices"));
  }

  mri->c_r = *MATRIX_RELT(c, 1, 1);
  mri->c_a = *MATRIX_RELT(c, 2, 1);
  mri->c_s = *MATRIX_RELT(c, 3, 1);

  MatrixFree(&origin);
  MatrixFree(&c);

  mri->ras_good_flag = 1;

  return(NO_ERROR);

} /* end apply_i_to_r() */

MATRIX *
MRIrasXformToVoxelXform(MRI *mri_src, MRI *mri_dst, MATRIX *m_ras_xform, 
                        MATRIX *m_voxel_xform)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras, *m_tmp ;

  if (!mri_dst)
    mri_dst = mri_src ;   /* assume they will be in the same space */

  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_src) ;
  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_dst) ;

  m_tmp = MatrixMultiply(m_ras_xform, m_voxel_to_ras, NULL) ;
  m_voxel_xform = MatrixMultiply(m_ras_to_voxel, m_tmp, m_voxel_xform) ;

  MatrixFree(&m_voxel_to_ras); MatrixFree(&m_ras_to_voxel); MatrixFree(&m_tmp);

  return(m_voxel_xform) ;
}

MATRIX *
MRIvoxelXformToRasXform(MRI *mri_src, MRI *mri_dst, MATRIX *m_voxel_xform, 
                        MATRIX *m_ras_xform)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras, *m_tmp ;

  if (!mri_dst)
    mri_dst = mri_src ;  /* assume they will be in the same space */

  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_src) ;
  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_dst) ;

  m_tmp = MatrixMultiply(m_voxel_xform, m_ras_to_voxel, NULL) ;
  m_ras_xform = MatrixMultiply(m_voxel_to_ras, m_tmp, m_ras_xform) ;

  MatrixFree(&m_voxel_to_ras); MatrixFree(&m_ras_to_voxel); MatrixFree(&m_tmp);

  return(m_ras_xform) ;
}

MATRIX *extract_r_to_i(MRI *mri)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras ;

  m_voxel_to_ras = extract_i_to_r(mri) ;
  m_ras_to_voxel = MatrixInverse(m_voxel_to_ras, NULL) ;
  MatrixFree(&m_voxel_to_ras) ;
  return(m_ras_to_voxel) ;
}
  
MATRIX *extract_i_to_r(MRI *mri)
{

  MATRIX *m;
  float m11, m12, m13, m14;
  float m21, m22, m23, m24;
  float m31, m32, m33, m34;
  float ci, cj, ck;

  m = MatrixAlloc(4, 4, MATRIX_REAL);
  if(m == NULL)
  {
    ErrorReturn(NULL, (ERROR_BADPARM, "extract_i_to_r(): error allocating matrix"));
  }

  if(mri->ras_good_flag)
  {
    m11 = mri->xsize * mri->x_r;  m12 = mri->ysize * mri->y_r;  m13 = mri->zsize * mri->z_r;
    m21 = mri->xsize * mri->x_a;  m22 = mri->ysize * mri->y_a;  m23 = mri->zsize * mri->z_a;
    m31 = mri->xsize * mri->x_s;  m32 = mri->ysize * mri->y_s;  m33 = mri->zsize * mri->z_s;

    ci = (mri->width - 1.0) / 2.0;
    cj = (mri->height - 1.0) / 2.0;
    ck = (mri->depth - 1.0) / 2.0;

    m14 = mri->c_r - (m11 * ci + m12 * cj + m13 * ck);
    m24 = mri->c_a - (m21 * ci + m22 * cj + m23 * ck);
    m34 = mri->c_s - (m31 * ci + m32 * cj + m33 * ck);

  }
  else if(mri->slice_direction == MRI_CORONAL)
  {
    m11 = -mri->xsize;  m12 =  0.0;         m13 = 0.0;         m14 =  mri->xsize * mri->width / 2.0;
    m21 =  0.0;         m22 =  0.0;         m23 = mri->zsize;  m24 = -mri->zsize * mri->depth / 2.0;
    m31 =  0.0;         m32 = -mri->ysize;  m33 = 0.0;         m34 =  mri->ysize * mri->height / 2.0;
  }
  else if(mri->slice_direction == MRI_SAGITTAL)
  {
    m11 =  0.0;         m12 =  0.0;         m13 = mri->zsize;  m14 = -mri->zsize * mri->depth / 2.0;
    m21 =  mri->xsize;  m22 =  0.0;         m23 = 0.0;         m24 = -mri->xsize * mri->width / 2.0;
    m31 =  0.0;         m32 = -mri->ysize;  m33 = 0.0;         m34 =  mri->ysize * mri->height / 2.0;
  }
  else if(mri->slice_direction == MRI_HORIZONTAL)
  {
    m11 = -mri->xsize;  m12 =  0.0;         m13 = 0.0;         m14 =  mri->xsize * mri->width / 2.0;
    m21 =  0.0;         m22 = -mri->ysize;  m23 = 0.0;         m24 =  mri->ysize * mri->height / 2.0;
    m31 =  0.0;         m32 =  0.0;         m33 = mri->zsize;  m34 = -mri->zsize * mri->depth / 2.0;
  }
  else
  {
    m11 = -mri->xsize;  m12 =  0.0;         m13 = 0.0;         
    m14 =  mri->xsize * mri->width / 2.0;
    m21 =  0.0;         m22 = -mri->ysize;  m23 = 0.0;         
    m24 =  mri->ysize * mri->height / 2.0;
    m31 =  0.0;         m32 =  0.0;         m33 = mri->zsize;  
    m34 = -mri->zsize * mri->depth / 2.0;
  }
  
  stuff_four_by_four(m, m11, m12, m13, m14, 
                        m21, m22, m23, m24, 
                        m31, m32, m33, m34, 
                        0.0, 0.0, 0.0, 1.0);


  return(m);

} /* end extract_i_to_r() */

/* eof */
MRI *
MRIscaleMeanIntensities(MRI *mri_src, MRI *mri_ref, MRI *mri_dst)
{
  int    width, height, depth, x, y, z, val ;
  double ref_mean, src_mean, nref_vox, nsrc_vox, scale ;
  
  mri_dst = MRIcopy(mri_src, mri_dst) ;

  width = mri_dst->width ; height = mri_dst->height ; depth = mri_dst->depth;

  nref_vox = nsrc_vox = src_mean = ref_mean = 0.0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_ref,x,y,z) > 10)
        {
          nref_vox++ ;
          ref_mean += (double)MRIvox(mri_ref, x, y, z) ;
        }
        if (MRIvox(mri_src,x,y,z) > 10)
        {
          src_mean += (double)MRIvox(mri_src, x, y, z) ;
          nsrc_vox++ ;
        }
      }
    }
  }

  ref_mean /= nref_vox ; src_mean /= nsrc_vox ; 
  fprintf(stderr, "mean brightnesses: ref = %2.1f, in = %2.1f\n",
          ref_mean, src_mean) ;
  scale = ref_mean / src_mean ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIvox(mri_src, x, y, z) ;
        val = nint(val*scale) ;
        if (val > 255)
          val = 255 ;
        MRIvox(mri_src, x, y, z) = val ;
      }
    }
  }

  return(mri_dst) ;
}

MRI *MRIsmoothParcellation(MRI *mri, int smooth_parcellation_count)
{

  MRI *mri2;
  int i, j, k;
  short vals[26];
  int counts[32768];
  int c;

  if(mri->type != MRI_SHORT)
  {
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIsmoothParcellation(): only supported for shorts data"));
  }

  mri2 = MRIcopy(mri, NULL);
  if(mri2 == NULL)
  {
    ErrorReturn(NULL, (ERROR_NOMEMORY, "MRIsmoothParcellation(): error copying structre"));
  }

  for(i = 1;i < mri->width-1;i++)
  {
    for(j = 1;j < mri->height-1;j++)
    {
      for(k = 1;k < mri->depth-1;k++)
      {

        memset(counts, 0x00, 32768 * sizeof(int));

        vals[ 0] = MRISvox(mri, i+1, j+1, k+1);
        vals[ 1] = MRISvox(mri, i+1, j+1, k  );
        vals[ 2] = MRISvox(mri, i+1, j+1, k-1);
        vals[ 3] = MRISvox(mri, i+1, j  , k+1);
        vals[ 4] = MRISvox(mri, i+1, j  , k  );
        vals[ 5] = MRISvox(mri, i+1, j  , k-1);
        vals[ 6] = MRISvox(mri, i+1, j-1, k+1);
        vals[ 7] = MRISvox(mri, i+1, j-1, k  );
        vals[ 8] = MRISvox(mri, i+1, j-1, k-1);

        vals[ 9] = MRISvox(mri, i  , j+1, k+1);
        vals[10] = MRISvox(mri, i  , j+1, k  );
        vals[11] = MRISvox(mri, i  , j+1, k-1);
        vals[12] = MRISvox(mri, i  , j  , k+1);
        /* --- ignore the voxel itself --- */
        vals[13] = MRISvox(mri, i  , j  , k-1);
        vals[14] = MRISvox(mri, i  , j-1, k+1);
        vals[15] = MRISvox(mri, i  , j-1, k  );
        vals[16] = MRISvox(mri, i  , j-1, k-1);

        vals[17] = MRISvox(mri, i-1, j+1, k+1);
        vals[18] = MRISvox(mri, i-1, j+1, k  );
        vals[19] = MRISvox(mri, i-1, j+1, k-1);
        vals[20] = MRISvox(mri, i-1, j  , k+1);
        vals[21] = MRISvox(mri, i-1, j  , k  );
        vals[22] = MRISvox(mri, i-1, j  , k-1);
        vals[23] = MRISvox(mri, i-1, j-1, k+1);
        vals[24] = MRISvox(mri, i-1, j-1, k  );
        vals[25] = MRISvox(mri, i-1, j-1, k-1);

        for(c = 0;c < 26;c++)
          counts[vals[c]]++;

        for(c = 0;c < 26;c++)
          if(counts[vals[c]] >= smooth_parcellation_count)
            MRISvox(mri2, i, j, k) = vals[c];

      }
    }
  }

  return(mri2);

} /* end MRIsmoothParcellation() */

int
MRIeraseBorderPlanes(MRI *mri)
{
  int  x, y, z ;

  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
    {
      MRIvox(mri, x, y, 0) = MRIvox(mri, x, y, mri->depth-1) = 0 ;
    }

  for (y = 0 ; y < mri->height ; y++)
    for (z = 0 ; z < mri->depth ; z++)
    {
      MRIvox(mri, 0, y, z) = MRIvox(mri, mri->width-1, y, z) = 0 ;
    }

  for (x = 0 ; x < mri->width ; x++)
    for (z = 0 ; z < mri->depth ; z++)
    {
      MRIvox(mri, x, 0, z) = MRIvox(mri, x, mri->height-1, z) = 0 ;
    }

  return(NO_ERROR) ;
}
int
MRIcopyPulseParameters(MRI *mri_src, MRI *mri_dst)
{
  mri_dst->flip_angle = mri_src->flip_angle ;
  mri_dst->tr = mri_src->tr ;
  mri_dst->te = mri_src->te ;
  mri_dst->ti = mri_src->ti ;
  return(NO_ERROR) ;
}
float
MRIfindNearestNonzero(MRI *mri, int wsize, Real xr, Real yr, Real zr)
{
  int   xk, yk, zk, xi, yi, zi, whalf, x, y, z ;
  float dist, min_dist, min_val, dx, dy, dz ;

  x = nint(xr) ; y = nint(yr) ; z = nint(zr) ;
  if (MRIvox(mri, x, y, z) > 0)
    return((float)MRIvox(mri, x, y, z)) ;

  min_dist = 100000 ; min_val = 0 ;
  whalf = (wsize-1)/2 ;
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    dz = zi-zr ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dy = yi-yr ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        dx = xi-xr ;
        if (MRIvox(mri, xi, yi, zi) > 0)
        {
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            min_val = MRIvox(mri, xi, yi, zi) ;
          }
        }
      }
    }
  }
  return(min_val) ;
}
float
MRIfindNearestNonzeroLocation(MRI *mri, int wsize, Real xr, Real yr, Real zr,
                              int *pxv, int *pyv, int *pzv)
{
  int   xk, yk, zk, xi, yi, zi, whalf, x, y, z ;
  float dist, min_dist, min_val, dx, dy, dz ;

  x = nint(xr) ; y = nint(yr) ; z = nint(zr) ;
  if (MRIvox(mri, x, y, z) > 0)
    return((float)MRIvox(mri, x, y, z)) ;

  min_dist = 100000 ; min_val = 0 ;
  whalf = (wsize-1)/2 ;
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    dz = zi-zr ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dy = yi-yr ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        dx = xi-xr ;
        if (MRIvox(mri, xi, yi, zi) > 0)
        {
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            if (pxv)
            { *pxv = xi ; *pyv = yi ; *pzv = zi ;}
            min_dist = dist ;
            min_val = MRIvox(mri, xi, yi, zi) ;
          }
        }
      }
    }
  }
  return(min_val) ;
}

