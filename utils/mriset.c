/*
 *       FILE NAME:   mriset.c
 *
 *       DESCRIPTION: utilities for MRI data structure
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        1/8/97
 *
 */

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>

#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "filter.h"
#include "box.h"
#include "region.h"

/*-----------------------------------------------------
  MACROS AND CONSTANTS
  -------------------------------------------------------*/

#define DEBUG_POINT(x,y,z)  (((x==72) && (y==142)) &&((z)==127))

/*-----------------------------------------------------
  STATIC DATA
  -------------------------------------------------------*/

/*-----------------------------------------------------
  STATIC PROTOTYPES
  -------------------------------------------------------*/
/*-----------------------------------------------------
  GLOBAL FUNCTIONS
  -------------------------------------------------------*/

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description: A sanity-check: make sure the volumes being
  operated upon are the same dimensions
  ------------------------------------------------------*/
static void MRIcheckVolDims(MRI *mri1, MRI *mri2)
{
  if (mri1 == NULL) return;
  if (mri2 == NULL) return;

  if (mri1->height != mri2->height)
    ErrorExit(ERROR_BADPARM, 
              "ERROR: %s-MRIcheckVolDims: volume1 height=%d != "
              "volume2 height=%d.\n",
              Progname, mri1->height, mri2->height);
  if (mri1->width != mri2->width)
    ErrorExit(ERROR_BADPARM, 
              "ERROR: %s-MRIcheckVolDims: volume1 width=%d != "
              "volume2 width=%d.\n",
              Progname, mri1->width, mri2->width);
  if (mri1->depth != mri2->depth)
    ErrorExit(ERROR_BADPARM, 
              "ERROR: %s-MRIcheckVolDims: volume1 depth=%d != "
              "volume2 depth=%d.\n",
              Progname, mri1->depth, mri2->depth);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIunion(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst, v1, v2 ;

  MRIcheckVolDims(mri1, mri2);

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri1, NULL) ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          p1 = &MRIvox(mri1, 0, y, z) ;
          p2 = &MRIvox(mri2, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              v1 = *p1++ ;
              v2 = *p2++ ;
              *pdst++ = MAX(v1, v2) ;
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
MRIintersect(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst, v1, v2 ;

  MRIcheckVolDims(mri1, mri2);

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri1, NULL) ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          p1 = &MRIvox(mri1, 0, y, z) ;
          p2 = &MRIvox(mri2, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              v1 = *p1++ ;
              v2 = *p2++ ;
              *pdst++ = MIN(v1, v2) ;
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
MRIcomplement(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, b ;

  MRIcheckVolDims(mri_src, mri_dst);

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
              *pdst++ = !b ;
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
MRIand(MRI *mri1, MRI *mri2, MRI *mri_dst, int thresh)
{
  int     width, height, depth, x, y, z, f ;
  Real    val1, val2 ;

  MRIcheckVolDims(mri1, mri2);

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri1, NULL) ;

  for (f = 0 ; f < mri1->nframes ; f++)
    {
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++)
                {
                  MRIsampleVolumeFrame(mri1, x, y, z, f, &val1) ;
                  if (val1 < thresh)
                    val1 = 0 ;
                  MRIsampleVolumeFrame(mri2, x, y, z, f, &val2) ;
                  if (val2 < thresh)
                    val2 = 0 ;
                  MRIsetVoxVal(mri_dst, x, y, z, f, val1 && val2) ;
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
MRIxor(MRI *mri1, MRI *mri2, MRI *mri_dst, int t1, int t2)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst, v1, v2 ;

  MRIcheckVolDims(mri1, mri2);

  if ((mri1->type != MRI_UCHAR) || (mri2->type != MRI_UCHAR))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIxor: inputs must be UCHAR")) ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri1, NULL) ;
  else if (mri_dst->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIxor: destination must be UCHAR")) ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          p1 = &MRIvox(mri1, 0, y, z) ;
          p2 = &MRIvox(mri2, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              v1 = *p1++ ;
              v2 = *p2++ ;
              v1 = ((v1 >= t1) && (v1 <= t2)) ;
              v2 = ((v2 >= t1) && (v2 <= t2)) ;
              *pdst++ = v1 ^ v2 ;
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
#define KLEN 3
MRI *
MRImorph(MRI *mri_src, MRI *mri_dst, int which)
{
  BUFTYPE kernel[KLEN][KLEN][KLEN] ;

  switch (which)
    {
    case FILTER_OPEN:
      memset(kernel, 1, sizeof(BUFTYPE)*KLEN*KLEN*KLEN) ;
      break ;
    case FILTER_CLOSE:
      break ;
    case FILTER_DILATE:
      memset(kernel, 1, sizeof(BUFTYPE)*KLEN*KLEN*KLEN) ;
      break ;
    case FILTER_ERODE:
      memset(kernel, 1, sizeof(BUFTYPE)*KLEN*KLEN*KLEN) ;
      break ;
    default:
      ErrorReturn(NULL,
                  (ERROR_UNSUPPORTED, "MRImorph: unknown type %d", which)) ;
    }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIerode(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi, same ;
  BUFTYPE *pdst, min_val, val ;

  MRIcheckVolDims(mri_src, mri_dst);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_dst == mri_src)
    {
      same = 1 ;
      mri_dst = MRIclone(mri_src, NULL) ;
    }
  else
    same = 0 ;

  if (mri_src->type != MRI_UCHAR || mri_dst->type != MRI_UCHAR)
    {
      Real fmin_val, fval ;

      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++)
                {
                  fmin_val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
                  for (z0 = -1 ; z0 <= 1 ; z0++)
                    {
                      zi = mri_src->zi[z+z0] ;
                      for (y0 = -1 ; y0 <= 1 ; y0++)
                        {
                          yi = mri_src->yi[y+y0] ;
                          for (x0 = -1 ; x0 <= 1 ; x0++)
                            {
                              xi = mri_src->xi[x+x0] ;
                              fval = MRIgetVoxVal(mri_src, xi,yi,zi, 0) ;
                              if (fval < fmin_val)
                                fmin_val = fval ;
                            }
                        }
                    }
                  MRIsetVoxVal(mri_dst, x, y, z, 0, fmin_val) ;
                }
            }
        }
    }
  else
    {
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              pdst = &MRIvox(mri_dst, 0, y, z) ;
              for (x = 0 ; x < width ; x++)
                {
                  min_val = 255 ;
                  for (z0 = -1 ; z0 <= 1 ; z0++)
                    {
                      zi = mri_src->zi[z+z0] ;
                      for (y0 = -1 ; y0 <= 1 ; y0++)
                        {
                          yi = mri_src->yi[y+y0] ;
                          for (x0 = -1 ; x0 <= 1 ; x0++)
                            {
                              xi = mri_src->xi[x+x0] ;
                              val = MRIvox(mri_src, xi,yi,zi) ;
                              if (val < min_val)
                                min_val = val ;
                            }
                        }
                    }
                  *pdst++ = min_val ;
                }
            }
        }
    }
  if (same)
    {
      MRIcopy(mri_dst, mri_src) ;
      MRIfree(&mri_dst) ;
      mri_dst = mri_src ;
    }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdilateThreshLabel(MRI *mri_src, MRI *mri_val, MRI *mri_dst, int label,
                     int niter, int thresh)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi,
    xmin, xmax, ymin, ymax, zmin, zmax, i, nadded ;
  BUFTYPE *psrc, out_val, val ;
  MRI     *mri_tmp = NULL ;

  MRIcheckVolDims(mri_src, mri_val);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  /* get everything outside of bounding box */
  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (i = 0 ; i < niter ; i++)
    {
      nadded = 0 ;
      mri_tmp = MRIcopy(mri_dst, mri_tmp) ; /* will allocate first time */
      xmax = 0 ; xmin = width-1 ;
      ymax = 0 ; ymin = height-1 ;
      zmax = 0 ; zmin = depth-1 ;
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              psrc = &MRIvox(mri_src, 0, y, z) ;
              for (x = 0 ; x < width ; x++)
                {
                  if (*psrc++ == label)
                    {
                      if (x-1 < xmin)
                        xmin = x-1 ;
                      if (x+1 > xmax)
                        xmax = x+1 ;
                      if (y-1 < ymin)
                        ymin = y-1 ;
                      if (y+1 > ymax)
                        ymax = y+1 ;
                      if (z-1 < zmin)
                        zmin = z-1 ;
                      if (z+1 > zmax)
                        zmax = z+1 ;
                    }
                }
            }
        }
      xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
      ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
      zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;
      for (z = zmin ; z <= zmax ; z++)
        {
          for (y = ymin ; y <= ymax ; y++)
            {
              for (x = xmin ; x <= xmax ; x++)
                {
                  out_val = MRIvox(mri_tmp, x, y, z) ;
                  if (MRIvox(mri_val, x, y, z) >= thresh)
                    {
                      for (z0 = -1 ; z0 <= 1 ; z0++)
                        {
                          zi = mri_tmp->zi[z+z0] ;
                          for (y0 = -1 ; y0 <= 1 ; y0++)
                            {
                              yi = mri_tmp->yi[y+y0] ;
                              for (x0 = -1 ; x0 <= 1 ; x0++)
                                {
                                  xi = mri_tmp->xi[x+x0] ;
                                  val = MRIvox(mri_tmp, xi,yi,zi) ;
                                  if (val == label)
                                    out_val = label ;
                                }
                            }
                        }
                      if (out_val == label)
                        {
                          if (xi < xmin)
                            xmin = xi ;
                          if (xi > xmax)
                            xmax = xi ;
                          if (yi < ymin)
                            ymin = yi ;
                          if (yi > ymax)
                            ymax = yi ;
                          if (zi < zmin)
                            zmin = zi ;
                          if (zi > zmax)
                            zmax = zi ;
                          if (MRIvox(mri_dst, x,y,z) != label)
                            nadded++ ;
                        }
                    }
                  MRIvox(mri_dst, x,y,z) = out_val ;
                }
            }
        }
      if (nadded == 0)
        break ;
    }
  MRIfree(&mri_tmp) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdilateInvThreshLabel(MRI *mri_src, MRI *mri_val, MRI *mri_dst, int label,
                        int niter, int thresh)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi,
    xmin, xmax, ymin, ymax, zmin, zmax, i, nadded ;
  BUFTYPE *psrc, out_val, val ;
  MRI     *mri_tmp = NULL ;

  MRIcheckVolDims(mri_src, mri_val);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  /* get everything outside of bounding box */
  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (i = 0 ; i < niter ; i++)
    {
      nadded = 0 ;
      mri_tmp = MRIcopy(mri_dst, mri_tmp) ; /* will allocate first time */
      xmax = 0 ; xmin = width-1 ;
      ymax = 0 ; ymin = height-1 ;
      zmax = 0 ; zmin = depth-1 ;
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              psrc = &MRIvox(mri_src, 0, y, z) ;
              for (x = 0 ; x < width ; x++)
                {
                  if (*psrc++ == label)
                    {
                      if (x-1 < xmin)
                        xmin = x-1 ;
                      if (x+1 > xmax)
                        xmax = x+1 ;
                      if (y-1 < ymin)
                        ymin = y-1 ;
                      if (y+1 > ymax)
                        ymax = y+1 ;
                      if (z-1 < zmin)
                        zmin = z-1 ;
                      if (z+1 > zmax)
                        zmax = z+1 ;
                    }
                }
            }
        }
      xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
      ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
      zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;
      for (z = zmin ; z <= zmax ; z++)
        {
          for (y = ymin ; y <= ymax ; y++)
            {
              for (x = xmin ; x <= xmax ; x++)
                {
                  out_val = MRIvox(mri_tmp, x, y, z) ;
                  if (MRIvox(mri_val, x, y, z) <= thresh)
                    {
                      for (z0 = -1 ; z0 <= 1 ; z0++)
                        {
                          zi = mri_tmp->zi[z+z0] ;
                          for (y0 = -1 ; y0 <= 1 ; y0++)
                            {
                              yi = mri_tmp->yi[y+y0] ;
                              for (x0 = -1 ; x0 <= 1 ; x0++)
                                {
                                  xi = mri_tmp->xi[x+x0] ;
                                  val = MRIvox(mri_tmp, xi,yi,zi) ;
                                  if (val == label)
                                    out_val = label ;
                                }
                            }
                        }
                      if (out_val == label)
                        {
                          if (xi < xmin)
                            xmin = xi ;
                          if (xi > xmax)
                            xmax = xi ;
                          if (yi < ymin)
                            ymin = yi ;
                          if (yi > ymax)
                            ymax = yi ;
                          if (zi < zmin)
                            zmin = zi ;
                          if (zi > zmax)
                            zmax = zi ;
                          if (MRIvox(mri_dst, x,y,z) != label)
                            nadded++ ;
                        }
                    }
                  MRIvox(mri_dst, x,y,z) = out_val ;
                }
            }
        }
      if (nadded == 0)
        break ;
    }
  MRIfree(&mri_tmp) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdilateLabelUchar(MRI *mri_src, MRI *mri_dst, int label, int niter)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi,
    xmin, xmax, ymin, ymax, zmin, zmax, i;
  BUFTYPE *psrc, out_val, val ;
  MRI     *mri_tmp = NULL ;

  MRIcheckVolDims(mri_src, mri_dst);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  /* get everything outside of bounding box */
  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (i = 0 ; i < niter ; i++)
    {
      mri_tmp = MRIcopy(mri_dst, mri_tmp) ; /* will allocate first time */
      xmax = 0 ; xmin = width-1 ;
      ymax = 0 ; ymin = height-1 ;
      zmax = 0 ; zmin = depth-1 ;
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              psrc = &MRIvox(mri_src, 0, y, z) ;
              for (x = 0 ; x < width ; x++)
                {
                  if (*psrc++ == label)
                    {
                      if (x-1 < xmin)
                        xmin = x-1 ;
                      if (x+1 > xmax)
                        xmax = x+1 ;
                      if (y-1 < ymin)
                        ymin = y-1 ;
                      if (y+1 > ymax)
                        ymax = y+1 ;
                      if (z-1 < zmin)
                        zmin = z-1 ;
                      if (z+1 > zmax)
                        zmax = z+1 ;
                    }
                }
            }
        }
      xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
      ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
      zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;
      for (z = zmin ; z <= zmax ; z++)
        {
          for (y = ymin ; y <= ymax ; y++)
            {
              for (x = xmin ; x <= xmax ; x++)
                {
                  out_val = MRIvox(mri_tmp, x, y, z) ;
                  for (z0 = -1 ; z0 <= 1 ; z0++)
                    {
                      if (out_val == label)
                        break ;
                      zi = mri_tmp->zi[z+z0] ;
                      for (y0 = -1 ; y0 <= 1 ; y0++)
                        {
                          if (out_val == label)
                            break ;
                          yi = mri_tmp->yi[y+y0] ;
                          for (x0 = -1 ; x0 <= 1 ; x0++)
                            {
                              if (out_val == label)
                                break ;
                              xi = mri_tmp->xi[x+x0] ;
                              val = MRIvox(mri_tmp, xi,yi,zi) ;
                              if (val == label)
                                out_val = label ;
                            }
                        }
                    }
                  MRIvox(mri_dst, x,y,z) = out_val ;
                }
            }
        }
    }
  MRIfree(&mri_tmp) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdilateLabel(MRI *mri_src, MRI *mri_dst, int label, int niter)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi,
    xmin, xmax, ymin, ymax, zmin, zmax, i, f ;
  MRI     *mri_tmp = NULL ;
  Real    out_val, val ;

  MRIcheckVolDims(mri_src, mri_dst);

  if (mri_src->type == MRI_UCHAR)
    return(MRIdilateLabelUchar(mri_src, mri_dst, label, niter)) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  /* get everything outside of bounding box */
  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (i = 0 ; i < niter ; i++)
    {
      mri_tmp = MRIcopy(mri_dst, mri_tmp) ; /* will allocate first time */
      xmax = 0 ; xmin = width-1 ;
      ymax = 0 ; ymin = height-1 ;
      zmax = 0 ; zmin = depth-1 ;
      for (f = 0 ; f < mri_src->nframes ; f++)
        {
          for (z = 0 ; z < depth ; z++)
            {
              for (y = 0 ; y < height ; y++)
                {
                  for (x = 0 ; x < width ; x++)
                    {
                      MRIsampleVolumeFrameType
                        (mri_src, x, y, z, f, SAMPLE_NEAREST, &val) ;
                      if (val == label)
                        {
                          if (x-1 < xmin)
                            xmin = x-1 ;
                          if (x+1 > xmax)
                            xmax = x+1 ;
                          if (y-1 < ymin)
                            ymin = y-1 ;
                          if (y+1 > ymax)
                            ymax = y+1 ;
                          if (z-1 < zmin)
                            zmin = z-1 ;
                          if (z+1 > zmax)
                            zmax = z+1 ;
                        }
                    }
                }
            }
        }
      xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
      ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
      zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;
      for (f = 0 ; f < mri_src->nframes ; f++)
        {
          for (z = zmin ; z <= zmax ; z++)
            {
              for (y = ymin ; y <= ymax ; y++)
                {
                  for (x = xmin ; x <= xmax ; x++)
                    {
                      MRIsampleVolumeFrameType
                        (mri_src, x, y, z, f, SAMPLE_NEAREST, &out_val) ;
                      for (z0 = -1 ; z0 <= 1 ; z0++)
                        {
                          if (out_val == label)
                            break ;
                          zi = mri_tmp->zi[z+z0] ;
                          for (y0 = -1 ; y0 <= 1 ; y0++)
                            {
                              if (out_val == label)
                                break ;
                              yi = mri_tmp->yi[y+y0] ;
                              for (x0 = -1 ; x0 <= 1 ; x0++)
                                {
                                  if (out_val == label)
                                    break ;
                                  xi = mri_tmp->xi[x+x0] ;
                                  MRIsampleVolumeFrameType
                                    (mri_src,x,y,z,f,SAMPLE_NEAREST, &val) ;
                                  if (val == label)
                                    {
                                      break ;
                                      out_val = label ;
                                    }
                                }
                            }
                        }
                    }
                  MRIsetVoxVal(mri_dst, x, y, z, f, out_val) ;
                }
            }
        }
    }
  MRIfree(&mri_tmp) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdilateUchar(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi, same,
    xmin, xmax, ymin, ymax, zmin, zmax;
  BUFTYPE *psrc, max_val, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIcheckVolDims(mri_src, mri_dst);

  if (mri_dst == mri_src)
	{
		same = 1 ;
		mri_dst = MRIclone(mri_src, NULL) ;
	}
  else
    same = 0 ;

  xmax = 0 ; xmin = width-1 ;
  ymax = 0 ; ymin = height-1 ;
  zmax = 0 ; zmin = depth-1 ;
  for (z = 0 ; z < depth ; z++)
	{
		for (y = 0 ; y < height ; y++)
		{
			psrc = &MRIvox(mri_src, 0, y, z) ;
			for (x = 0 ; x < width ; x++)
			{
				if (*psrc++ > 0)
				{
					if (x-1 < xmin)
						xmin = x-1 ;
					if (x+1 > xmax)
						xmax = x+1 ;
					if (y-1 < ymin)
						ymin = y-1 ;
					if (y+1 > ymax)
						ymax = y+1 ;
					if (z-1 < zmin)
						zmin = z-1 ;
					if (z+1 > zmax)
						zmax = z+1 ;
				}
			}
		}
	}
  xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
  ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
  zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;
  for (z = zmin ; z <= zmax ; z++)
	{
		for (y = ymin ; y <= ymax ; y++)
		{
			for (x = xmin ; x <= xmax ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)
					DiagBreak() ;
				max_val = 0 ;
				for (z0 = -1 ; z0 <= 1 ; z0++)
				{
					zi = mri_src->zi[z+z0] ;
					for (y0 = -1 ; y0 <= 1 ; y0++)
					{
						yi = mri_src->yi[y+y0] ;
						for (x0 = -1 ; x0 <= 1 ; x0++)
						{
							xi = mri_src->xi[x+x0] ;
							val = MRIvox(mri_src, xi,yi,zi) ;
							if (val > max_val)
								max_val = val ;
						}
					}
				}
				MRIvox(mri_dst, x,y,z) = max_val ;
			}
		}
	}
  if (same)
	{
		MRIcopy(mri_dst, mri_src) ;
		MRIfree(&mri_dst) ;
		mri_dst = mri_src ;
	}
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdilate(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi, same,
    xmin, xmax, ymin, ymax, zmin, zmax, f;
  Real    val, max_val ;

  if (mri_src->type == MRI_UCHAR)
    return(MRIdilateUchar(mri_src, mri_dst)) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIcheckVolDims(mri_src, mri_dst);

  if (mri_dst == mri_src)
    {
      same = 1 ;
      mri_dst = MRIclone(mri_src, NULL) ;
    }
  else
    same = 0 ;

#if 1
  xmax = 0 ; xmin = width-1 ;
  ymax = 0 ; ymin = height-1 ;
  zmax = 0 ; zmin = depth-1 ;
  for (f = 0 ; f < mri_src->nframes ; f++)
    {
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++)
                {
                  MRIsampleVolumeFrameType
                    (mri_src, x, y, z, f, SAMPLE_NEAREST, &val) ;
                  if (val > 0)
                    {
                      if (x-1 < xmin)
                        xmin = x-1 ;
                      if (x+1 > xmax)
                        xmax = x+1 ;
                      if (y-1 < ymin)
                        ymin = y-1 ;
                      if (y+1 > ymax)
                        ymax = y+1 ;
                      if (z-1 < zmin)
                        zmin = z-1 ;
                      if (z+1 > zmax)
                        zmax = z+1 ;
                    }
                }
            }
        }
    }
  xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
  ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
  zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;
#else
  xmin = 0 ; xmax = width-1 ;
  ymin = 0 ; ymax = height-1 ;
  zmin = 0 ; zmax = depth-1 ;
#endif
  for (f = 0 ; f < mri_src->nframes ; f++)
    {
      for (z = zmin ; z <= zmax ; z++)
        {
          for (y = ymin ; y <= ymax ; y++)
            {
              for (x = xmin ; x <= xmax ; x++)
                {
                  max_val = 0 ;
                  for (z0 = -1 ; z0 <= 1 ; z0++)
                    {
                      zi = mri_src->zi[z+z0] ;
                      for (y0 = -1 ; y0 <= 1 ; y0++)
                        {
                          yi = mri_src->yi[y+y0] ;
                          for (x0 = -1 ; x0 <= 1 ; x0++)
                            {
                              xi = mri_src->xi[x+x0] ;
                              MRIsampleVolumeFrameType
                                (mri_src,xi,yi,zi,f,SAMPLE_NEAREST,&val)  ;
                              if (val > max_val)
                                max_val = val ;
                            }
                        }
                    }
                  MRIsetVoxVal(mri_dst, x,y,z, f, max_val) ;
                }
            }
        }
    }
  if (same)
    {
      MRIcopy(mri_dst, mri_src) ;
      MRIfree(&mri_dst) ;
      mri_dst = mri_src ;
    }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIopen(MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_tmp ;

  mri_tmp = MRIerode(mri_src, NULL) ;
  mri_dst = MRIdilate(mri_tmp, mri_dst) ;
  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIclose(MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_tmp ;

  mri_tmp = MRIdilate(mri_src, NULL) ;
  mri_dst = MRIerode(mri_tmp, mri_dst) ;
  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIerode6(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi;
  BUFTYPE *pdst, min_val, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIcheckVolDims(mri_src, mri_dst);

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              xi = mri_src->xi[x] ;
              yi = mri_src->yi[y] ;
              min_val = 255 ;
              for (z1 = -1 ; z1 <= 1 ; z1++)
                {
                  zi = mri_src->zi[z+z1] ;
                  val = MRIvox(mri_src, xi, yi, zi) ;
                  if (val < min_val)
                    min_val = val ;
                }
              zi = mri_src->zi[z] ;
              for (y1 = -1 ; y1 <= 1 ; y1++)
                {
                  if (!y1)    /* already done */
                    continue ;
                  yi = mri_src->yi[y+y1] ;
                  val = MRIvox(mri_src, xi, yi, zi) ;
                  if (val < min_val)
                    min_val = val ;
                }
              yi = mri_src->yi[y] ;
              for (x1 = -1 ; x1 <= 1 ; x1++)
                {
                  if (!y1)    /* already done */
                    continue ;
                  xi = mri_src->xi[x+x1] ;
                  val = MRIvox(mri_src, xi, yi, zi) ;
                  if (val < min_val)
                    min_val = val ;
                }
              *pdst++ = min_val ;
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
MRIdilate6(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi;
  BUFTYPE *pdst, max_val, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIcheckVolDims(mri_src, mri_dst);

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              xi = mri_src->xi[x] ;
              yi = mri_src->yi[y] ;
              max_val = 0 ;
              for (z1 = -1 ; z1 <= 1 ; z1++)
                {
                  zi = mri_src->zi[z+z1] ;
                  val = MRIvox(mri_src, xi, yi, zi) ;
                  if (val > max_val)
                    max_val = val ;
                }
              zi = mri_src->zi[z] ;
              for (y1 = -1 ; y1 <= 1 ; y1++)
                {
                  if (!y1)    /* already done */
                    continue ;
                  yi = mri_src->yi[y+y1] ;
                  val = MRIvox(mri_src, xi, yi, zi) ;
                  if (val > max_val)
                    max_val = val ;
                }
              yi = mri_src->yi[y] ;
              for (x1 = -1 ; x1 <= 1 ; x1++)
                {
                  if (!y1)    /* already done */
                    continue ;
                  xi = mri_src->xi[x+x1] ;
                  val = MRIvox(mri_src, xi, yi, zi) ;
                  if (val > max_val)
                    max_val = val ;
                }
              *pdst++ = max_val ;
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
MRIopen6(MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_tmp ;

  mri_tmp = MRIerode6(mri_src, NULL) ;
  mri_dst = MRIdilate6(mri_tmp, mri_dst) ;
  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIclose6(MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_tmp ;

  mri_tmp = MRIdilate6(mri_src, NULL) ;
  mri_dst = MRIerode6(mri_tmp, mri_dst) ;
  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIerodeRegion(MRI *mri_src, MRI *mri_dst, int wsize, MRI_REGION *region)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0, val,
    xk, yk, zk, xi, yi, zi, min_val ;
  BUFTYPE *pdst ;

  MRIcheckVolDims(mri_src, mri_dst);

  whalf = wsize/2 ;
  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

  if (!mri_dst)
    {
      int  w, h, d ;

      w = width - region->x ;
      h = height - region->y ;
      d = depth - region->z ;
      mri_dst = MRIalloc(w, h, d, MRI_UCHAR) ;
      MRIcopyHeader(mri_src, mri_dst) ;
      mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize ;
      mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize ;
      mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize ;
      mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
      mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
      mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
    }

  for (z = z0 ; z < depth ; z++)
    {
      for (y = y0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y-y0, z-z0) ;
          for (x = x0 ; x < width ; x++)
            {
              min_val = 255 ;
              for (zk = -whalf ; zk <= whalf ; zk++)
                {
                  zi = mri_src->zi[z+zk] ;
                  for (yk = -whalf ; yk <= whalf ; yk++)
                    {
                      yi = mri_src->yi[y+yk] ;
                      for (xk = -whalf ; xk <= whalf ; xk++)
                        {
                          xi = mri_src->xi[x+xk] ;
                          val = (int)MRIvox(mri_src, xi, yi, zi) ;
                          if (val < min_val)
                            min_val = val ;
                        }
                    }
                }
              *pdst++ = min_val ;
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
MRIcomputeResidual(MRI *mri1, MRI *mri2, MRI *mri_dst, int t1, int t2)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst, v1, v2, out ;

  if ((mri1->type != MRI_UCHAR) || (mri2->type != MRI_UCHAR))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIxor: inputs must be UCHAR")) ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri1, NULL) ;
  else if (mri_dst->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "MRIcomputeResidual: destination must be UCHAR")) ;

  MRIcheckVolDims(mri1, mri2);

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          p1 = &MRIvox(mri1, 0, y, z) ;
          p2 = &MRIvox(mri2, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              v1 = *p1++ ;
              v2 = *p2++ ;
              v1 = ((v1 >= t1) && (v1 <= t2)) ;
              v2 = ((v2 >= t1) && (v2 <= t2)) ;
              if (v1 && !v2)
                out = 255 ;   /* v2 off and v1 on - make it 'positive' */
              else if (!v1 && v2)
                out = 0 ;     /* v2 on and v1 off - make it 'negative' */
              else
                out = 128 ;   /* both on or both off */
              *pdst++ = out ;
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
MRIminmax(MRI *mri_src, MRI *mri_dst, MRI *mri_dir, int wsize)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi, whalf, offset;
  float   max_val, val, min_val ;
  float   fmin, fmax ;

  MRIcheckVolDims(mri_src, mri_dst);

  MRIvalRange(mri_src, &fmin, &fmax) ;
  whalf = (wsize-1) / 2 ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              offset = MRIgetVoxVal(mri_dir, x, y, z, 0);
              max_val = fmin ;
              min_val = fmax ;
              if (offset != OFFSET_ZERO) for (z0 = -whalf ; z0 <= whalf ; z0++)
                {
                  zi = mri_src->zi[z+z0] ;
                  for (y0 = -whalf ; y0 <= whalf ; y0++)
                    {
                      yi = mri_src->yi[y+y0] ;
                      for (x0 = -whalf ; x0 <= whalf ; x0++)
                        {
                          xi = mri_src->xi[x+x0] ;
                          val = MRIgetVoxVal(mri_src, xi,yi,zi,0) ;
                          if (val > max_val)
                            max_val = val ;
                          if (val < min_val)
                            min_val = val ;
                        }
                    }
                }
              switch (offset)
                {
                case OFFSET_GRADIENT_DIRECTION:
                  val = max_val ;
                  break ;
                case OFFSET_NEGATIVE_GRADIENT_DIRECTION:
                  val = min_val ;
                  break ;
                default:
                case OFFSET_ZERO:
                  val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
                  break ;
                }
              MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
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
MRIreplaceValues(MRI *mri_src, MRI *mri_dst, float in_val, float out_val)
{
  int     width, height, depth, x, y, z;
  float   val ;

  MRIcheckVolDims(mri_src, mri_dst);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_src->type == MRI_UCHAR && mri_dst->type == MRI_UCHAR)
    return(MRIreplaceValuesUchar
           (mri_src, mri_dst, (BUFTYPE)nint(in_val), (BUFTYPE)nint(out_val))) ;
  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
              if (FEQUAL(val, in_val))
                val = out_val ;
              MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
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
MRIreplaceValuesUchar(MRI *mri_src, MRI *mri_dst,
                      BUFTYPE in_val, BUFTYPE out_val)
{
  int     width, height, depth, x, y, z;
  BUFTYPE  val ;
  BUFTYPE *psrc, *pdst ;

  MRIcheckVolDims(mri_src, mri_dst);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_src->type != MRI_UCHAR || mri_dst->type != MRI_UCHAR)
    return(MRIreplaceValues(mri_src, mri_dst, (float)in_val, (float)out_val)) ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          psrc = &MRIvox(mri_src, 0, y, z) ;
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              val = *psrc++ ;
              if (val == in_val)
                val = out_val ;
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
MRImeanMask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,int mask, int wsize)
{
  int     width, height, depth, x, y, z, mask_val;
  BUFTYPE *pmask ;
  float   val ;

  MRIcheckVolDims(mri_src, mri_mask);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_src->type != mri_dst->type)
    ErrorReturn
      (NULL,
       (ERROR_UNSUPPORTED, "MRImeanMask: src and dst must be same type")) ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pmask = &MRIvox(mri_mask, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              if (x == Gx && y == Gy && z == Gz)
                DiagBreak() ;
              mask_val = MRIgetVoxVal(mri_mask, x, y, z, 0) ;
              val = MRIgetVoxVal(mri_src, x, y, z, 0) ;

              if (mask_val == mask)
                {
                  val = MRIvoxelMean(mri_src, x, y, z, wsize) ;
                }
              MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
            }
        }
    }
  return(mri_dst) ;
}
MRI   *MRImaskDifferentGeometry(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, 
																int mask, float out_val)
{
	MATRIX *m_vox2vox ;
	VECTOR *v1, *v2 ;
	int    x, y, z, xd, yd, zd, f, mask_val, val ;

	v1 = VectorAlloc(4, MATRIX_REAL) ;
	v2 = VectorAlloc(4, MATRIX_REAL) ;
	VECTOR_ELT(v1, 4) = 1.0 ;
	VECTOR_ELT(v2, 4) = 1.0 ;
	m_vox2vox = MRIgetVoxelToVoxelXform(mri_src, mri_mask) ;

  if(!mri_dst) 
		mri_dst = MRIclone(mri_src, NULL) ;
	for (x = 0 ; x < mri_dst->width ; x++)
	{
		V3_X(v1) = x ;
		for (y = 0 ; y < mri_dst->height ; y++)
		{
			V3_Y(v1) = y ;
			for (z = 0 ; z < mri_dst->depth ; z++)
			{
				V3_Z(v1) = z ;
				MatrixMultiply(m_vox2vox, v1, v2) ;
				xd = nint(V3_X(v2)) ;  yd = nint(V3_Y(v2)) ;  zd = nint(V3_Z(v2)) ; 
				if (xd < 0 || xd >= mri_mask->width ||
						yd < 0 || yd >= mri_mask->height ||
						zd < 0 || zd >= mri_mask->depth )
					mask_val = mask+1 ; // allow it through
				else
					mask_val = nint(MRIgetVoxVal(mri_mask, xd, yd, zd,0)) ;
				for(f = 0; f < mri_dst->nframes; f++)
				{
					val = nint(MRIgetVoxVal(mri_src, x, y, z, f)) ;
          if(mask_val == mask) 
						val = out_val ;
					else 
						val = MRIgetVoxVal(mri_src, x, y, z, f) ;
					MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
				}
			}
		}
	}

	VectorFree(&v1) ; VectorFree(&v2) ; MatrixFree(&m_vox2vox) ;
	return(mri_dst) ;
}

/*--------------------------------------------------------------------------
  MRImask() - applies mask to an mri data set. If mri_mask == mask at a voxel,
  then that voxel is set to out_val. Eg, if your mask is binary (0 and 1),
  and you want to set voxels outside the mask to (ie, maskmri=0) 0, then
  MRImask(src,mask,dst,0,0). Handles multiple frames.
  --------------------------------------------------------------------------*/
MRI *
MRImask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int mask, float out_val)
{
  int     width, height, depth, nframes, x, y, z, f, mask_val;
  float   val ;

	if (mri_src->width != mri_mask->width ||
			mri_src->height != mri_mask->height ||
			mri_src->depth != mri_mask->depth)
		return(MRImaskDifferentGeometry(mri_src, mri_mask, mri_dst,mask,out_val));

  MRIcheckVolDims(mri_src, mri_mask);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  nframes = mri_src->nframes ;

  if(!mri_dst) mri_dst = MRIclone(mri_src, NULL) ;

  if(mri_src->type != mri_dst->type)
    ErrorReturn
      (NULL,
       (ERROR_UNSUPPORTED, "MRImask: src and dst must be same type")) ;

  for (z = 0 ; z < depth ; z++){
    for (y = 0 ; y < height ; y++){
      for (x = 0 ; x < width ; x++){
				if (x == Gx && y == Gy && z == Gz) DiagBreak() ;
				mask_val = MRIgetVoxVal(mri_mask, x, y, z, 0) ;
				for(f = 0; f < nframes; f++){
          if(mask_val == mask) val = out_val ;
					else val = MRIgetVoxVal(mri_src, x, y, z, f) ;
					MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*------------------------------------------------------------------
  MRImaskInvert() - changes the 1s to 0s and vice versa.
  ------------------------------------------------------------------*/
MRI *MRImaskInvert(MRI *mask, MRI *outmask)
{
  int c,r,s,f;
  double v;

  MRIcheckVolDims(mask, outmask);

  if(outmask == NULL) outmask = MRIcopy(mask,NULL);

  for(c=0; c < mask->width; c++){
    for(r=0; r < mask->height; r++){
      for(s=0; s < mask->depth; s++){
        for(f=0; f < mask->nframes; f++){
          v = MRIgetVoxVal(mask,c,r,s,f);
          if(v > 0.5) MRIsetVoxVal(outmask,c,r,s,f,0);
          else        MRIsetVoxVal(outmask,c,r,s,f,1);
        }
      }
    }
  }
  return(outmask);
}

/*------------------------------------------------------------------
  MRInMask() - count the number of voxels in the mask. "In the mask"
  means > 0.5.
  ------------------------------------------------------------------*/
int MRInMask(MRI *mask)
{
  int c,r,s,f,nmask;
  double v;

  nmask = 0;
  for(c=0; c < mask->width; c++){
    for(r=0; r < mask->height; r++){
      for(s=0; s < mask->depth; s++){
        for(f=0; f < mask->nframes; f++){
          v = MRIgetVoxVal(mask,c,r,s,f);
          if(v > 0.5) nmask++;
        }
      }
    }
  }
  return(nmask);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIthresholdMask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                 float mask_threshold,float out_val)
{
  int     width, height, depth, x, y, z;
  BUFTYPE *pdst, *psrc, *pmask, val, mask_val ;
  float   fval, fmask ;

  MRIcheckVolDims(mri_src, mri_mask);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_src->type == MRI_UCHAR && mri_mask->type == MRI_UCHAR)
    {
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              pdst = &MRIvox(mri_dst, 0, y, z) ;
              psrc = &MRIvox(mri_src, 0, y, z) ;
              pmask = &MRIvox(mri_mask, 0, y, z) ;
              for (x = 0 ; x < width ; x++)
                {
                  val = *psrc++ ;
                  mask_val = *pmask++ ;
                  if (mask_val >= mask_threshold)
                    *pdst++ = val ;
                  else
                    *pdst++ = out_val ;
                }
            }
        }
    }
  else
    {
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++)
                {
                  fval = MRIgetVoxVal(mri_src, x, y, z, 0) ;
                  fmask = MRIgetVoxVal(mri_mask, x, y, z, 0) ;
                  if (fmask < mask_threshold)
                    fval = out_val ;
                  MRIsetVoxVal(mri_dst, x, y, z, 0,fval) ;
                }
            }
        }
    }
  return(mri_dst) ;
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
MRI *
MRImaskThreshold(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, float threshold,
                 int out_label)
{
  BUFTYPE   *pmask, *pdst, *psrc, out_val, mask, in_val ;
  int       width, height, depth, x, y, z, nchanged, noff, non ;
  float     nvox ;

  MRIcheckVolDims(mri_src, mri_mask);

  if (mri_mask->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "MRI3Dthreshold: mask must be MRI_FLOAT")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
  */


  noff = non = nchanged = 0 ;
  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pmask = &MRIvox(mri_mask, 0, y, z) ;
          psrc = &MRIvox(mri_src, 0, y, z) ;
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              if (DEBUG_POINT(x,y,z))
                DiagBreak() ;  /* marked as 255, should be 127 */
              out_val = 0 ;
              mask = *pmask++ ;   /* value from inverse morphed volume */
              in_val = *psrc++ ;
              if (mask < 100-threshold)        /* probably off */
                out_val = 0 ;
              else  if (mask > threshold)      /* probably on */
                out_val = out_label ;
              else                    /* not sure, use original fill val */
                out_val = in_val ;
              if (out_val != in_val)
                {
                  if (out_val)
                    non++ ;
                  else
                    noff++ ;
                  nchanged++ ;
                }
              *pdst++ = out_val ;
            }
        }
    }

  nvox = (float)(width * height * depth) ;
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)nchanged/nvox) ;
  fprintf(stderr, "%8d of %8d voxels off     - %2.1f%%.\n",
          noff, (int)nvox,   100.0f*(float)noff/nvox) ;
  fprintf(stderr, "%8d of %8d voxels on      - %2.1f%%.\n",
          nchanged, (int)nvox, 100.0f*(float)non/nvox) ;
  return(mri_dst) ;
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
int
MRIgrowLabel(MRI *mri, MRI *mri_filled, int in_label, int out_label)
{
  int      x, y, z, width, height, depth, xi, yi, zi, xk, yk, zk, nfilled,
    total_filled, xmin, xmax, ymin, ymax, zmin, zmax ;

  MRIcheckVolDims(mri, mri_filled);

  width = mri->width ; height = mri->height ; depth = mri->depth ;
  total_filled = 0 ;

  xmin = width ; ymin = height ; zmin = depth ;
  xmax = ymax = zmax = 0 ;
  if (in_label) for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              if (MRIvox(mri, x, y, z) == in_label)
                {
                  if (x > xmax)
                    xmax = x ;
                  if (y > ymax)
                    ymax = y ;
                  if (z > zmax)
                    zmax = z ;
                  if (x < xmin)
                    xmin = x ;
                  if (y < ymin)
                    ymin = y ;
                  if (z < zmin)
                    zmin = z ;
                }
            }
        }
    }
  else   /* filling bg - do outside first (hack, but much faster)*/
    {
      /* find bounding box for real data */
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++)
                {
                  if (MRIvox(mri, x, y, z))
                    {
                      if (x > xmax)
                        xmax = x ;
                      if (y > ymax)
                        ymax = y ;
                      if (z > zmax)
                        zmax = z ;
                      if (x < xmin)
                        xmin = x ;
                      if (y < ymin)
                        ymin = y ;
                      if (z < zmin)
                        zmin = z ;
                    }
                }
            }
        }

      /* fill everything outside it */
      for (z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++)
                {
                  if (z <= zmin || z >= zmax ||
                      y <= ymin || y >= ymax ||
                      x <= xmin || x >= xmax)
                    {
                      total_filled++ ;
                      MRIvox(mri_filled, x, y, z) = out_label ;
                    }
                }
            }
        }
    }

#if 0
  xmin = ymin = zmin = 0 ;
  xmax = width-1 ; ymax = height-1 ; zmax = depth-1 ;
#endif

  do
    {
      nfilled = 0 ;
      for (z = zmin ; z <= zmax ; z++)
        {
          for (y = ymin ; y <= ymax ; y++)
            {
              for (x = xmin ; x <= xmax ; x++)
                {
                  if (MRIvox(mri_filled, x, y, z) == out_label)
                    {
                      for (zk = -1 ; zk <= 1 ; zk++)
                        {
                          zi = z+zk ;
                          if (zi < 0 || zi >= depth)
                            continue ;
                          for (yk = -1 ; yk <= 1 ; yk++)
                            {
                              yi = y+yk ;
                              if (yi < 0 || yi >= height)
                                continue ;
                              for (xk = -1 ; xk <= 1 ; xk++)
                                {
                                  xi = x+xk ;
                                  if (xi < 0 || xi >= width)
                                    continue ;
                                  if ((MRIvox(mri, xi, yi, zi) ==
                                       in_label) &&
                                      (MRIvox(mri_filled, xi, yi, zi) !=
                                       out_label))
                                    {
                                      nfilled++ ;
                                      MRIvox(mri_filled, xi, yi, zi) =
                                        out_label ;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
      total_filled += nfilled ;
      fprintf(stderr,
              "%d voxels filled, total = %d.\n", nfilled, total_filled);
    } while (nfilled > 0) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
int
MRIturnOnFG(MRI *mri, MRI *mri_fg, MRI *mri_bg)
{
  int    x, y, z, width, height, depth ;

  MRIcheckVolDims(mri, mri_fg);

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              if (x == 157 && y == 140 && z == 56)
                DiagBreak() ;
              if (MRIvox(mri_fg, x, y, z) > 0)
                MRIvox(mri, x, y, z) = MRIvox(mri_fg, x, y, z) ;
            }
        }
    }

  /* now remove islands */
  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              if (x == 157 && y == 140 && z == 56)
                DiagBreak() ;
              if ((MRIvox(mri_fg, x, y, z) == 0) && (MRIvox(mri,x,y,z) > 0))
                {
                  MRIvox(mri, x, y, z) = 0 ;
                  MRIvox(mri_bg, x, y, z) = 1 ;
                }
            }
        }
    }
  return(NO_ERROR) ;
}
static BUFTYPE findLabel(MRI *mri, int x0, int y0, int z0) ;
static BUFTYPE
findLabel(MRI *mri, int x, int y, int z)
{
  int   xi, yi, zi, xk, yk, zk, label, width, height, depth, counts[256],
    max_label, max_count, whalf = 1 ;

  if (x == 148 && y == 104 && z == 136)
    DiagBreak() ;

  do
    {
      memset(counts, 0, sizeof(counts)) ;

      width = mri->width ; height = mri->height ; depth = mri->depth ;
      for (zk = -whalf ; zk <= whalf ; zk++)
        {
          zi = z + zk ;
          if (zi < 0 || zi >= depth)
            continue ;
          for (yk = -whalf ; yk <= whalf ; yk++)
            {
              yi = y + yk ;
              if (yi < 0 || yi >= height)
                continue ;
              for (xk = -whalf ; xk <= whalf ; xk++)
                {
                  xi = x + xk ;
                  if (xi < 0 || xi >= width)
                    continue ;
                  label = MRIvox(mri, xi, yi, zi) ;
                  if (label <= WM_MIN_VAL)
                    continue ;
                  counts[label]++ ;
                }
            }
        }

      max_count = max_label = 0 ;
      for (label = 0 ; label <= 255 ; label++)
        if (counts[label] > max_count)
          {
            max_count = counts[label] ;
            max_label = label ;
          }
      whalf++ ;
    } while (max_count == 0 && whalf < 5) ;

  return(max_label) ;
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  turn off all voxels which are set in the bg image
  ----------------------------------------------------------------------*/
int
MRIturnOffBG(MRI *mri, MRI *mri_bg)
{
  int    x, y, z, width, height, depth ;

  MRIcheckVolDims(mri, mri_bg);

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              if (x == 157 && y == 140 && z == 56)
                DiagBreak() ;
              if (MRIvox(mri_bg, x, y, z) > 0)
                MRIvox(mri, x, y, z) = 0 ;
            }
        }
    }

  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              if (x == 157 && y == 140 && z == 56)
                DiagBreak() ;
              if ((MRIvox(mri_bg, x, y, z) == 0) &&
                  (MRIvox(mri, x, y, z) == 0))
                MRIvox(mri, x, y, z) = findLabel(mri, x, y, z) ;
            }
        }
    }
  return(NO_ERROR) ;
}
#if 0
static BUFTYPE
findLabel(MRI *mri, int x, int y, int z)
{
  int   xi, yi, zi, xk, yk, zk, left_count, right_count, label,
    width, height, depth ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;
  left_count = right_count = 0 ;
  for (zk = -1 ; zk <= 1 ; zk++)
    {
      zi = z + zk ;
      if (zi < 0 || zi >= depth)
        continue ;
      for (yk = -1 ; yk <= 1 ; yk++)
        {
          yi = y + yk ;
          if (yi < 0 || yi >= height)
            continue ;
          for (xk = -1 ; xk <= 1 ; xk++)
            {
              xi = x + xk ;
              if (xi < 0 || xi >= width)
                continue ;
              label = MRIvox(mri, xi, yi, zi) ;
              if (label == MRI_LEFT_HEMISPHERE)
                left_count++ ;
              else if (label == MRI_RIGHT_HEMISPHERE)
                right_count++ ;
            }
        }
    }
  return ((left_count > right_count) ?
          MRI_LEFT_HEMISPHERE :
          MRI_RIGHT_HEMISPHERE) ;
}
#endif
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Set an MRI intensity values to 0
  ------------------------------------------------------*/
int
MRIsetValues(MRI *mri, float val)
{
  int   width, depth, height, x, y, z, frame, nframes ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  nframes = mri->nframes ;

  for (frame = 0 ; frame < nframes ; frame++)
    {
      for (x = 0 ; x < width  ; x++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (z = 0 ; z < depth ; z++)
                {
                  MRIsetVoxVal(mri, x, y, z, frame, val) ;
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
  ------------------------------------------------------*/
MRI *
MRIsoapBubbleLabel(MRI *mri_src, MRI *mri_label, MRI *mri_dst, int label,
                   int niter)
{
  int     width, height, depth, x, y, z, x0, y0, z0, xi, yi, zi,
    xmin, xmax, ymin, ymax, zmin, zmax, i, n;
  BUFTYPE *plabel ;
  MRI     *mri_tmp = NULL ;
  float   mean ;

  MRIcheckVolDims(mri_src, mri_label);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  xmax = 0 ; xmin = width-1 ;
  ymax = 0 ; ymin = height-1 ;
  zmax = 0 ; zmin = depth-1 ;
  for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          plabel = &MRIvox(mri_label, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              if (*plabel++ == label)
                {
                  if (x-1 < xmin)
                    xmin = x-1 ;
                  if (x+1 > xmax)
                    xmax = x+1 ;
                  if (y-1 < ymin)
                    ymin = y-1 ;
                  if (y+1 > ymax)
                    ymax = y+1 ;
                  if (z-1 < zmin)
                    zmin = z-1 ;
                  if (z+1 > zmax)
                    zmax = z+1 ;
                }
            }
        }
    }
  xmin = MAX(0, xmin) ; xmax = MIN(width-1, xmax) ;
  ymin = MAX(0, ymin) ; ymax = MIN(height-1, ymax) ;
  zmin = MAX(0, zmin) ; zmax = MIN(depth-1, zmax) ;

  for (i = 0 ; i < niter ; i++)
    {
      mri_tmp = MRIcopy(mri_dst, mri_tmp) ; /* will allocate first time */
      for (z = zmin ; z <= zmax ; z++)
        {
          for (y = ymin ; y <= ymax ; y++)
            {
              for (x = xmin ; x <= xmax ; x++)
                {
                  if (x == 97 && y == 131 && z == 181)
                    DiagBreak() ;
                  if (x == 97 && y == 131 && z == 181 && i > 171)
                    DiagBreak() ;
                  if (MRIvox(mri_label, x, y, z) == label)
                    {
                      mean = MRIvox(mri_tmp, x, y, z) ; n = 1 ;
                      for (z0 = -1 ; z0 <= 1 ; z0++)
                        {
                          zi = mri_tmp->zi[z+z0] ;
                          for (y0 = -1 ; y0 <= 1 ; y0++)
                            {
                              yi = mri_tmp->yi[y+y0] ;
                              for (x0 = -1 ; x0 <= 1 ; x0++)
                                {
                                  xi = mri_tmp->xi[x+x0] ;
                                  mean += MRIvox(mri_tmp, xi,yi,zi) ;
                                  n++ ;
                                }
                            }
                        }
                      MRIvox(mri_dst, x,y,z) = mean / n ;
                    }
                }
            }
        }
    }
  MRIfree(&mri_tmp) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Copy the label from one volume to another
  ------------------------------------------------------*/
int
MRIcopyLabel(MRI *mri_src, MRI *mri_dst, int label)
{
  int     width, height, depth, x, y, z, nvox ;
  BUFTYPE *psrc, *pdst ;

  MRIcheckVolDims(mri_src, mri_dst);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (mri_src->type == MRI_UCHAR && mri_dst->type == MRI_UCHAR)
    {
      for (nvox = z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              psrc = &MRIvox(mri_src, 0, y, z) ;
              pdst = &MRIvox(mri_dst, 0, y, z) ;
              for (x = 0 ; x < width ; x++, pdst++)
                {
                  if (*psrc++ == label)
                    {
                      *pdst = label ;
                      nvox++ ;
                    }
                }
            }
        }
    }
  else
    {
      for (nvox = z = 0 ; z < depth ; z++)
        {
          for (y = 0 ; y < height ; y++)
            {
              for (x = 0 ; x < width ; x++, pdst++)
                {
                  if (MRIgetVoxVal(mri_src, x, y, z, 0) == label)
                    {
                      MRIsetVoxVal(mri_dst, x, y, z, 0, label) ;
                      nvox++ ;
                    }
                }
            }
        }
    }
  return(nvox) ;
}
int
MRIlabelOverlap(MRI *mri1, MRI *mri2, int label)
{
  int     width, height, depth, x, y, z, nvox ;
  BUFTYPE *pbuf1, *pbuf2 ;

  MRIcheckVolDims(mri1, mri2);

  width = mri1->width ; height = mri1->height ; depth = mri1->depth ;

  for (nvox = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pbuf1 = &MRIvox(mri1, 0, y, z) ;
          pbuf2 = &MRIvox(mri2, 0, y, z) ;
          for (x = 0 ; x < width ; x++, pbuf1++, pbuf2++)
            {
              if (*pbuf1 == label && *pbuf2 == label)
                nvox++ ;
            }
        }
    }
  return(nvox) ;
}
int
MRItotalVoxelsOn(MRI *mri, int thresh)
{
  int     width, height, depth, x, y, z, nvox ;
  BUFTYPE *pbuf ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (nvox = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pbuf = &MRIvox(mri, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              if (*pbuf++ > thresh)
                nvox++ ;
            }
        }
    }
  return(nvox) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Count and return the # of voxels with a given label.
  ------------------------------------------------------*/
int
MRIvoxelsInLabel(MRI *mri, int label)
{
  int     width, height, depth, x, y, z, nvox, val ;
  BUFTYPE *pbuf ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (nvox = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          pbuf = &MRIvox(mri, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
            {
              val = nint(MRIgetVoxVal(mri, x, y, z, 0)) ;
              if (val == label)
                nvox++ ;
            }
        }
    }
  return(nvox) ;
}

MRI *
MRInot(MRI *mri_src, MRI *mri_dst)
{
  int   x, y, z, out ;
  float val ;

  MRIcheckVolDims(mri_src, mri_dst);

  if (mri_dst == NULL)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ;  y < mri_src->height ; y++)
        {
          for (z = 0 ; z < mri_src->depth ; z++)
            {
              if (x == Gx && y == Gy && z == Gz)
                DiagBreak() ;
              val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
              out = !nint(val) ;
              MRIsetVoxVal(mri_dst, x, y, z, 0, out) ;
            }
        }
    }
  return(mri_dst) ;
}
int
MRIcopyLabeledVoxels(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, int label)
{
  int     width, height, depth, x, y, z, nvox ;
  float   val ;

  MRIcheckVolDims(mri_src, mri_labeled);

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (nvox = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
            {
              if (x == Gx && y == Gy && z == Gz)
                DiagBreak() ;
              if (MRIgetVoxVal(mri_labeled, x, y, z, 0) == label)
                {
                  val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
                  MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
                  nvox++ ;
                }
            }
        }
    }
  return(nvox) ;
}
MRI *
MRIsetLabelValues(MRI *mri_src, MRI *mri_label, MRI *mri_dst,
                  int label, float val)
{
  int     x, y, z, l ;

  MRIcheckVolDims(mri_src, mri_label);

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (x = 0 ; x < mri_dst->width; x++)
    {
      for (y = 0 ; y < mri_dst->height; y++)
        {
          for (z = 0 ; z < mri_dst->depth; z++)
            {
              l = MRIgetVoxVal(mri_label, x, y, z, 0) ;
              if (l == label)
                MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
            }
        }
    }
  return(mri_dst) ;
}

/*---------------------------------------------------------------
  MRIsphereMask() - creates a binary spherical mask centered at
  voxel (c0,r0,s0) with radius voxradius (measured in voxels).
  Voxels in the sphere will have a value of val, those outside
  will be 0. All frames will have the same mask.
  ---------------------------------------------------------------*/
MRI *MRIsphereMask(int ncols, int nrows, int nslices, int nframes,
                   int c0, int r0, int s0, double voxradius, double val,
                   MRI *mri)
{
  int r,c,s,f;
  double d2, r2, v;

  r2 = voxradius*voxradius;

  if(mri==NULL) mri = MRIallocSequence(ncols,nrows,nslices,MRI_FLOAT,nframes);

  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
        d2 = (c-c0)*(c-c0) + (r-r0)*(r-r0) + (s-s0)*(s-s0);
        if(d2 > r2) v = 0;
        else        v = val;
        for(f=0; f < mri->nframes; f++) MRIsetVoxVal(mri,c,r,s,f,v);
      }
    }
  }
  return(mri);
}
