/*
 *       FILE NAME:   mrihisto.c
 *
 *       DESCRIPTION: utilities for MRI  data structure histograms
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
#include "nr.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define DEBUG_POINT(x,y,z)  (((x) == 15)&&((y)==6)&&((z)==15))

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/
static HISTOGRAM *mriHistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo,
                                     MRI_REGION *region);
/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIapplyHistogram(MRI *mri_src, MRI *mri_dst, HISTOGRAM *histo)
{
  int       width, height, depth, x, y, z ;
  BUFTYPE   *psrc, *pdst, sval, dval ;

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
        sval = *psrc++ ;
        if (sval >= histo->nbins)
          sval = histo->nbins - 1 ;
        dval = histo->counts[sval] ;
        *pdst++ = dval ;
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
HISTOGRAM *
MRIhistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region)
{
  int               overlap ;
  float             fmin, fmax, bin_size ;
  BUFTYPE           bmin, bmax ;
  static MRI        *mri_prev = NULL ;
  static HISTOGRAM  h_prev ;
  static MRI_REGION reg_prev ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED,"MRIhistogramRegion: must by type UCHAR"));

  fmin = MRIvalRange(mri, &fmin, &fmax) ;
  bmin = (BUFTYPE)fmin ; bmax = (BUFTYPE)fmax ;
  if (!nbins)
    nbins = bmax - bmin + 1 ;

  if (!histo)
    histo = HISTOalloc(nbins) ;
  else
    histo->nbins = nbins ;

  HISTOclear(histo, histo) ;
  bin_size = (fmax - fmin + 1) / (float)nbins ;

  if (!mri_prev)   /* first invocation, initialize state machine */
  {
    HISTOclear(&h_prev, &h_prev) ;
    REGIONclear(&reg_prev) ;
  }

/*
   note that the overlap only works with subsequent windows advancing only 
   in the x direction.
 */
  /* check to see if regions overlap */
  overlap = ((mri == mri_prev) &&
             (region->x > reg_prev.x) && 
             (region->y == reg_prev.y) &&
             (region->z == reg_prev.z)) ;

  if (overlap)   /* take advantage of overlapping regions */
  {
    MRI_REGION  reg_left, reg_right ;
    HISTOGRAM   histo_left, histo_right ;

    HISTOclear(&histo_left, &histo_left) ;
    HISTOclear(&histo_right, &histo_right) ;
    histo_left.nbins = histo_right.nbins = 256 ;
    REGIONcopy(&reg_prev, &reg_left) ;
    reg_left.dx = region->x - reg_left.x ;
    mriHistogramRegion(mri, 0, &histo_left, &reg_left) ;
    REGIONcopy(&reg_prev, &reg_right) ;
    reg_right.x = reg_prev.x + reg_prev.dx ;
    reg_right.dx = region->x + region->dx - reg_right.x ;
    mriHistogramRegion(mri, 0, &histo_right, &reg_right) ;

    HISTOsubtract(&h_prev, &histo_left, histo) ;
    HISTOadd(histo, &histo_right, histo) ;
  }
  else
    mriHistogramRegion(mri, nbins, histo, region) ;
  
  mri_prev = mri ;
  HISTOcopy(histo, &h_prev) ;
  REGIONcopy(region, &reg_prev) ;
  return(histo) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static HISTOGRAM *
mriHistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region)
{
  int        width, height, depth, x, y, z, bin_no, x0, y0, z0 ;
  float      fmin, fmax, bin_size ;
  BUFTYPE    val, *psrc, bmin, bmax ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED,"MRIhistogramRegion: must by type UCHAR"));

  fmin = MRIvalRange(mri, &fmin, &fmax) ;
  bmin = (BUFTYPE)fmin ; bmax = (BUFTYPE)fmax ;
  if (!nbins)
    nbins = bmax - bmin + 1 ;

  if (!histo)
    histo = HISTOalloc(nbins) ;
  else
    histo->nbins = nbins ;

  HISTOclear(histo, histo) ;

  bin_size = (fmax - fmin + 1) / (float)nbins ;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

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

  for (bin_no = 0 ; bin_no < nbins ; bin_no++)
    histo->bins[bin_no] = (bin_no+1)*bin_size ;
  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, x0, y, z) ;
      for (x = x0 ; x < width ; x++)
      {
        val = *psrc++ ;
        bin_no = (int)((float)(val - bmin) / (float)bin_size) ;
        histo->counts[bin_no]++ ;
      }
    }
  }
  return(histo) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIhistoEqualizeRegion(MRI *mri_src, MRI *mri_dst, int low,MRI_REGION *region)
{
  HISTOGRAM  histo_eq ;

  MRIgetEqualizeHistoRegion(mri_src, &histo_eq, low, region, 1) ;
  mri_dst = MRIapplyHistogramToRegion(mri_src, mri_dst, &histo_eq, region) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIapplyHistogramToRegion(MRI *mri_src, MRI *mri_dst, 
                          HISTOGRAM *histo, MRI_REGION *region)
{
  int       width, height, depth, x, y, z, x0, y0, z0 ;
  BUFTYPE   *psrc, *pdst, sval, dval ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

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

  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, x0, y, z) ;
      pdst = &MRIvox(mri_dst, x0, y, z) ;
      for (x = x0 ; x < width ; x++)
      {
        sval = *psrc++ ;
        if (sval >= histo->nbins)
          sval = histo->nbins - 1 ;
        dval = histo->counts[sval] ;
        *pdst++ = dval ;
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
HISTOGRAM *
MRIgetEqualizeHistoRegion(MRI *mri, HISTOGRAM *histo_eq, int low, 
                          MRI_REGION *region, int norm)
{
  int       i, total, total_pix, *pc, *pdst, nbins ;
  HISTOGRAM histo ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, 
                 "MRIgetEqualizeHistoRegion: unsupported type %d",mri->type));
  MRIhistogramRegion(mri, 0, &histo, region) ;
  nbins = histo.nbins ;
  if (!histo_eq)
    histo_eq = HISTOalloc(nbins) ;
  else
  {
    histo_eq->nbins = nbins ;
    HISTOclear(histo_eq, histo_eq) ;
  }

  for (pc = &histo.counts[0], total_pix = 0, i = low ; i < nbins ; i++)
    total_pix += *pc++ ;

  if (total_pix) 
  {
    pc = &histo.counts[0] ;
    pdst = &histo_eq->counts[0] ;

    for (total = 0, i = low ; i < nbins ; i++)
    {
      total += *pc++ ;
      if (norm)
        *pdst++ = nint(255.0f * (float)total / (float)total_pix) ;
      else
        *pdst++ = total ;
    }
  }

  return(histo_eq) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
HISTOGRAM *
MRIgetEqualizeHisto(MRI *mri, HISTOGRAM *histo_eq, int low, int norm)
{
  int       i, total, total_pix ;
  HISTOGRAM *histo ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIgetEqualizeHisto: unsupported type %d",
                 mri->type)) ;
  histo = MRIhistogram(mri, 0) ;
  if (!histo_eq)
    histo_eq = HISTOalloc(histo->nbins) ;

  for (total_pix = 0, i = low ; i < histo->nbins ; i++)
    total_pix += histo->counts[i] ;

  for (total = 0, i = low ; i < histo->nbins ; i++)
  {
    total += histo->counts[i] ;
    if (norm)
      histo_eq->counts[i] = nint(255.0f * (float)total / (float)total_pix) ;
    else
      histo_eq->counts[i] = total ;
  }

  HISTOfree(&histo) ;
  return(histo_eq) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIhistoEqualize(MRI *mri_src, MRI *mri_dst, int low)
{
  HISTOGRAM  histo_eq ;

  MRIgetEqualizeHisto(mri_src, &histo_eq, low, 1) ;
  mri_dst = MRIapplyHistogram(mri_src, mri_dst, &histo_eq) ;
/*  MRIcrunch(mri_dst, mri_dst) ;*/
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIcrunch(MRI *mri_src, MRI *mri_dst)
{
  HISTOGRAM  *histo ;
  int        b, deleted ;

  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIcrunch: unsupported type %d",
                 mri_src->type)) ;
  histo = MRIhistogram(mri_src, 0) ;

  deleted = 0 ;
  for (b = 0 ; b < histo->nbins ; b++)
  {
    if (!histo->counts[b])
      deleted++ ;
    histo->counts[b] = b-deleted ;
  }

  histo->nbins -= deleted ;
  mri_dst = MRIapplyHistogram(mri_src, mri_dst, histo) ;
  HISTOfree(&histo) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
HISTOGRAM *
MRIhistogram(MRI *mri, int nbins)
{
  int        width, height, depth, x, y, z, bin_no ;
  HISTOGRAM  *histo ;
  float      fmin, fmax, bin_size ;
  BUFTYPE    val, *psrc, bmin, bmax ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,"MRIhistogram: must by type UCHAR"));

  fmin = MRIvalRange(mri, &fmin, &fmax) ;
  bmin = (BUFTYPE)fmin ; bmax = (BUFTYPE)fmax ;
  if (!nbins)
    nbins = bmax - bmin + 1 ;

  histo = HISTOalloc(nbins) ;

  bin_size = (fmax - fmin + 1) / (float)nbins ;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (bin_no = 0 ; bin_no < nbins ; bin_no++)
    histo->bins[bin_no] = (bin_no+1)*bin_size ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        bin_no = (int)((float)(val - bmin) / (float)bin_size) ;
        histo->counts[bin_no]++ ;
      }
    }
  }
  return(histo) ;
}
