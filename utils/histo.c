/**
 * @file  histo.c
 * @brief 1d histogram utilities.
 *
 * Utilities for computing and analyzing 1 dimensional histograms.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/08/06 16:53:16 $
 *    $Revision: 1.54 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
#include "histo.h"
#include "diag.h"
#include "macros.h"
#include "fio.h"

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfree(HISTOGRAM **phisto)
{
  HISTOGRAM *histo ;

  histo = *phisto ;
  *phisto = NULL ;
  if (histo)
  {
    if (histo->bins)
    {
      free(histo->bins) ;
      histo->bins = NULL;
    }
    else
      DiagBreak() ;
    if (histo->counts)
    {
      free(histo->counts) ;
      histo->counts = NULL;
    }
    else
      DiagBreak() ;
    free(histo) ;
  }
  else
    DiagBreak() ;
  
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOdump(HISTOGRAM *histo, FILE *fp)
{
  int  bin_no ;

  if (!histo)
    fprintf(stderr, "NULL histogram") ;
  else
  {
    fprintf(fp, "nbins = %d\n", histo->nbins) ;
    for (bin_no = 0 ; bin_no < histo->nbins ; bin_no++)
      if (histo->counts[bin_no])
        fprintf(fp, "bin[%d] = %2.1f = %2.2f\n",
                bin_no, histo->bins[bin_no], histo->counts[bin_no]) ;
  }
  return(NO_ERROR) ;
}


HISTOGRAM *
HISTOrealloc(HISTOGRAM *histo, int nbins)
{
  if (histo->bins)
    free(histo->bins) ;
  if (histo->counts)
    free(histo->counts) ;
  histo->bins = (float *)calloc(nbins, sizeof(float)) ;
  histo->counts = (float *)calloc(nbins, sizeof(float)) ;
  if (!histo->counts || !histo->bins)
    ErrorExit(ERROR_NOMEMORY, 
              "HISTOrealloc(%d): could not allocate histogram",
              nbins) ;
  histo->nbins = nbins ;

  return(histo) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOalloc(int nbins)
{
  HISTOGRAM *histo ;

  histo = (HISTOGRAM *)calloc(1, sizeof(HISTOGRAM)) ;
  if (!histo)
    ErrorExit(ERROR_NO_MEMORY, "HISTOalloc(%d): allocation failed", nbins) ;

  histo->bins = (float *)calloc(nbins, sizeof(float)) ;
  histo->counts = (float *)calloc(nbins, sizeof(float)) ;
  // fprintf(stderr, "histo->bins and ->counts allocated %d bins\n", nbins);
  if (!histo->counts || !histo->bins)
    ErrorExit(ERROR_NOMEMORY, 
              "HISTOalloc(%d): could not allocate histogram",
              nbins) ;
  histo->nbins = nbins ;

  return(histo) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOcrunch(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  int  b, deleted ;

  if (!histo_dst)
    histo_dst = HISTOalloc(histo_src->nbins) ;

  deleted = 0 ;
  for (b = 0 ; b < histo_src->nbins ; b++)
  {
    if (!histo_src->counts[b])
      deleted++ ;
    else
    {
      histo_dst->counts[b-deleted] = histo_src->counts[b] ;
      histo_dst->bins[b-deleted] = histo_src->bins[b] ;
    }
  }

  histo_dst->nbins -= deleted ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOcopy(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  if (!histo_dst)
    histo_dst = HISTOalloc(histo_src->nbins) ;
  histo_dst->nbins = histo_src->nbins ;
  histo_dst->bin_size = histo_src->bin_size ;
  memcpy(histo_dst->counts, 
         histo_src->counts, 
         sizeof(*histo_src->counts)*histo_src->nbins) ;
  memcpy(histo_dst->bins, 
         histo_src->bins, 
         sizeof(*histo_src->bins)*histo_src->nbins) ;
  histo_dst->min = histo_src->min ;
  histo_dst->max = histo_src->max ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Treat histo_src as a function mapping input intensities
  to output intensities, and numerically invert it.
  ------------------------------------------------------*/
HISTOGRAM *
HISTOinvert(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, int max_dst)
{
  int       b, bdst, val, max_count ;

  histo_dst = HISTOclear(histo_src, histo_dst) ;

  if (!max_dst)
    max_dst = 255 ;

  for (max_count = b = 0 ; b < histo_src->nbins ; b++)
  {
    val = histo_src->counts[b];
    if (val > max_count)
      max_count = val ;
  }

  for (b = 0 ; b < histo_src->nbins ; b++)
  {
    val = histo_src->counts[b] ;
    bdst = nint((float)max_dst * (float)val / (float)max_count) ;
    if (bdst > max_dst)
      bdst = max_dst ;
    if (!histo_dst->bins[bdst])
    {
      histo_dst->counts[bdst] = b ;
      histo_dst->bins[bdst] = 1 ;
    }
  }

  histo_dst->nbins = max_dst ;
  /*
     fill in zeros in the inverse histogram - they correspond to
     flat regions in the forward map (i.e. multiple-valued).
  */
  for (val = b = 0 ; b < histo_dst->nbins ; b++)
  {
    if (histo_dst->counts[b] > 0)
      val = histo_dst->counts[b] ;
    else
      histo_dst->counts[b] = val ;
  }

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, int max_out)
{
  int    max_count, b ;
  float  scale ;

  if (!max_out)
    max_out = 255 ;

  if (!histo_dst)
    histo_dst = HISTOalloc(histo_src->nbins) ;

  for (max_count = b = 0 ; b < histo_src->nbins ; b++)
    if (histo_src->counts[b] > max_count)
      max_count = histo_src->counts[b] ;

  scale = (float)max_out / (float)max_count ;
  for (b = 0 ; b < histo_src->nbins ; b++)
    histo_dst->counts[b] = (scale * (float)histo_src->counts[b]) ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOclear(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  if (!histo_dst)
    histo_dst = HISTOalloc(histo_src->nbins) ;

  memset(histo_dst->counts, 0, histo_dst->nbins*sizeof(*histo_dst->counts)) ;
  memset(histo_dst->bins, 0, histo_dst->nbins*sizeof(*histo_dst->bins)) ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOclearCounts(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  if (!histo_dst)
  {
    histo_dst = HISTOalloc(histo_src->nbins) ;
    memcpy(histo_dst->bins, histo_src->bins, sizeof(*histo_src->bins)) ;
  }

  memset(histo_dst->counts, 0, sizeof(histo_dst->counts)) ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  treating histo1 and histo2 as mappings from intensities to
  intensities (via the ->counts[] field), compose the two
  mappings into a composite transformation (first histo1 then
  histo2).
  ------------------------------------------------------*/
HISTOGRAM *
HISTOcompose(HISTOGRAM *histo1, HISTOGRAM *histo2, HISTOGRAM *histo_dst)
{
  int  b, val ;

  if (!histo_dst)
    histo_dst = HISTOalloc(histo1->nbins) ;

  for (b = 0 ; b < histo1->nbins ; b++)
  {
    val = histo1->counts[b] ;
    if (val >= histo2->nbins)
      val = histo2->nbins - 1 ;
    histo_dst->counts[b] = histo2->counts[val] ;
  }

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOfillZeros(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  int  b, val ;

  histo_dst = HISTOcopy(histo_src, histo_dst) ;
  /*
     fill in zeros in the inverse histogram - they correspond to
     flat regions in the forward map (i.e. multiple-valued).
  */
  for (val = b = 0 ; b < histo_dst->nbins ; b++)
  {
    if (histo_dst->counts[b] > 0)
      val = histo_dst->counts[b] ;
    else
      histo_dst->counts[b] = val ;
  }

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Make a compose two mappings - one forward and one
  inverse. This is useful in the context of histogram
  specification in which the forward is the equalization
  histogram of the source image, and the inverse is
  the equalization histogram of the template image.
  ------------------------------------------------------*/
#define MAX_STRETCH 2

HISTOGRAM *
HISTOcomposeInvert(HISTOGRAM *histo_fwd, HISTOGRAM *histo_inv,
                   HISTOGRAM *histo_dst)
{
  int   b, binv, val, max_fwd, max_inv, stretch ;
  float ffwd, next = 0.0f, prev ;

  if (!histo_dst)
    histo_dst = HISTOalloc(histo_fwd->nbins) ;
  else
  {
    HISTOclear(histo_dst, histo_dst) ;
    histo_dst->nbins = histo_fwd->nbins ;
  }

  for (max_fwd = b = 0 ; b < histo_fwd->nbins ; b++)
    if (histo_fwd->counts[b] > max_fwd)
      max_fwd = histo_fwd->counts[b] ;
  for (max_inv = b = 0 ; b < histo_inv->nbins ; b++)
    if (histo_inv->counts[b] > max_inv)
      max_inv = histo_inv->counts[b] ;

  if (!max_inv || !max_fwd)
    return(histo_dst) ;

  for (binv = b = 0 ; b < histo_fwd->nbins ; b++)
  {
    val = histo_fwd->counts[b] ;
    ffwd = (float)val / (float)max_fwd ;

    for (stretch = 0 ; stretch < MAX_STRETCH && binv < histo_dst->nbins ;
         stretch++, binv++)
    {
      next = (float)histo_inv->counts[binv] / (float)max_inv ;
      if (next >= ffwd)
        break ;
    }
    if (binv > 0)  /* see whether bin to left or to right is closer */
    {
      prev = (float)histo_inv->counts[binv-1] / (float)max_inv ;
      if (fabs(next-ffwd) > abs(prev-ffwd))  /* 'prev' bin is closer */
        binv-- ;
    }
    histo_dst->counts[b] = binv++ ;
  }

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOadd(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst)
{
  int  b ;
  float *pc1, *pc2, *pcdst ;

  if (!histo_dst)
    histo_dst = HISTOalloc(h1->nbins) ;
  else
    histo_dst->nbins = h1->nbins ;

  pc1 = &h1->counts[0] ;
  pc2 = &h2->counts[0] ;
  pcdst = &histo_dst->counts[0];
  for (b = 0 ; b < h1->nbins ; b++)
    *pcdst++ = *pc1++ + *pc2++ ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOmul(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst)
{
  int  b ;
  float *pc1, *pc2, *pcdst ;

  if (!histo_dst)
    histo_dst = HISTOalloc(h1->nbins) ;
  else
    histo_dst->nbins = h1->nbins ;

  pc1 = &h1->counts[0] ;
  pc2 = &h2->counts[0] ;
  pcdst = &histo_dst->counts[0];
  for (b = 0 ; b < h1->nbins ; b++)
    *pcdst++ = *pc1++ * *pc2++ ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOsubtract(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst)
{
  int   b ;
  float *pc1, *pc2, *pcdst ;

  if (!histo_dst)
    histo_dst = HISTOalloc(h1->nbins) ;
  else
    histo_dst->nbins = h1->nbins ;

  pc1 = &h1->counts[0] ;
  pc2 = &h2->counts[0] ;
  pcdst = &histo_dst->counts[0];
  for (b = 0 ; b < h1->nbins ; b++)
    *pcdst++ = *pc1++ - *pc2++ ;

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *
HISTOclearBins(HISTOGRAM *histo_src, 
               HISTOGRAM *histo_dst, 
               int min_val, 
               int max_val)
{
  int b ;

  if (!histo_src)
    return(NULL) ;

  if (!histo_dst || histo_dst != histo_src)
    histo_dst = HISTOcopy(histo_src, histo_dst) ;

  for (b = 0 ; b < histo_dst->nbins ; b++)
  {
    if (histo_dst->bins[b] >= min_val && histo_dst->bins[b] <= max_val)
    {
      histo_dst->counts[b] = 0 ;
    }
  }

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MAX_LEN 2000
HISTOGRAM *
HISTOsmooth(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,float sigma)
{
  float     norm, two_sigma, fx, k, kernel[MAX_LEN], total ;
  int       x, half, len, b, kx, b1, nbins ;

  nbins = histo_src->nbins ;
  if (!histo_dst)
  {
    histo_dst = HISTOalloc(nbins) ;
    histo_dst->bin_size = histo_src->bin_size ;
    histo_dst->min = histo_src->min ;
    histo_dst->max = histo_src->max ;
  }
  else
  {
    if (histo_dst->nbins < histo_src->nbins)
    {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTOrealloc(histo_dst, nbins);
    }
    histo_dst->nbins = nbins ;
    histo_dst->bin_size = histo_src->bin_size ;
  }
  /* build the kernel in k */
  len = (int)nint(8.0f * sigma)+1 ;
  if (ISEVEN(len))   /* ensure it's even */
    len++ ;
  if (MAX_LEN && (MAX_LEN < len))
    len = MAX_LEN ;
  half = len/2 ;

  norm = 0.0f ;
  two_sigma = 2.0f * sigma ;

  for (norm = 0.0f, x = 0 ; x < len ; x++)
  {
    fx = (float)(x-half) ;
    if (fabs(fx) <= two_sigma)
      k = (float)exp((double)(-fx*fx/(two_sigma*sigma))) ;
    else if (two_sigma < fabs(fx) && fabs(fx) <= 4.0f*sigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) *
          (float)pow(4.0f - fabs(fx)/(double)sigma, 4.0) ;
    else
      k = 0 ;

    kernel[x] = k ;
    norm += k ;
  }
  for (x = 0 ; x < len ; x++)
    kernel[x] /= norm ;

  for (b = 0 ; b < nbins ; b++)
  {
    for (total = 0.0f, x = 0 ; x < len ; x++)
    {
      kx = x - half ;
      b1 = b + kx ;
      if (b1 >= nbins || b1 < 0)
        continue ;
      total += kernel[x] * (float)histo_src->counts[b1] ;
    }
    histo_dst->counts[b] = total ;
    histo_dst->bins[b] = histo_src->bins[b] ;
  }

  return(histo_dst) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindLastPeakRelative(HISTOGRAM *h, int wsize, float min_pct)
{
  int  peak, b, bw, nbins, whalf, other_val ;
  float max_count, min_count, center_val ;

  peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins) ;
  if (peak < 0)
    return(-1) ;
  max_count = h->counts[peak] ;

  min_count = min_pct * max_count ;
  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = nbins-1 ; b >= 0 ; b--)
  {
    center_val = h->counts[b] ;
    if (center_val <= min_count)
      continue ;
    peak = 1 ;
    for (bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins)
        continue ;
      other_val = h->counts[bw] ;
      if (other_val > center_val)
        peak = 0 ;
    }
    if (peak)
      return(b) ;
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindFirstPeakRelative(HISTOGRAM *h, int wsize, float min_pct)
{
  int  peak, b, bw, nbins, whalf, other_val ;
  float max_count, min_count, center_val ;

  peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins) ;
  if (peak < 0)
    return(-1) ;
  max_count = h->counts[peak] ;

  min_count = min_pct * max_count ;
  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = 0 ; b < nbins ; b++)
  {
    center_val = h->counts[b] ;
    if (center_val <= min_count)
      continue ;
    peak = 1 ;
    for (bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins)
        continue ;
      other_val = h->counts[bw] ;
      if (other_val > center_val)
        peak = 0 ;
    }
    if (peak)
      return(b) ;
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindLastPeak(HISTOGRAM *h, int wsize, float min_pct)
{
  int  peak, b, bw, nbins, whalf ;
  float min_count, max_count, center_val, other_val ;

  for (max_count = b = 0 ; b < h->nbins ; b++)
  {
    center_val = h->counts[b];
    if (center_val > max_count)
      max_count = center_val ;
  }

  if (!max_count)
    return(-1) ;

  min_count = (min_pct * (float)max_count) ;
  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = nbins-1 ; b >= 0 ; b--)
  {
    center_val = h->counts[b] ;
    if (center_val <= min_count)
      continue ;
    peak = 1 ;
    for (bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins)
        continue ;
      other_val = h->counts[bw] ;
      if (other_val > center_val)
        peak = 0 ;
    }
    if (peak)
      return(b) ;
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindFirstPeak(HISTOGRAM *h, int wsize, float min_pct)
{
  int  peak, b, bw, nbins, whalf ;
  float center_val, max_count, other_val, min_count ;

  for (max_count = b = 0 ; b < h->nbins ; b++)
  {
    center_val = h->counts[b];
    if (center_val > max_count)
      max_count = center_val ;
  }

  if (!max_count)
    return(-1) ;

  min_count = min_pct * max_count ;
  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = 0 ; b < nbins ; b++)
  {
    center_val = h->counts[b] ;
    if (center_val <= min_count)
      continue ;
    peak = 1 ;
    for (bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins)
        continue ;
      other_val = h->counts[bw] ;
      if (other_val > center_val)
        peak = 0 ;
#if 0
      else if (center_val > other_val)
        peak++ ;
#endif
    }
    if (peak)
      return(b) ;
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindValley(HISTOGRAM *h, int wsize, int I0, int I1)
{
  int  valley, b, bw, nbins, whalf, center_val, max_count, other_val, b0, b1 ;

  // find max, and find bins corresponding to I0 and I1
  b0 = b1 = -1 ;
  for (max_count = b = 0 ; b < h->nbins ; b++)
  {
    if (h->bins[b] >= I0 && b0 == -1)
      b0 = b > 0 ? b-1 : b ;
    if (h->bins[b] >= I1 && b1 == -1)
      b1 = b < h->nbins-1 ? b+1 : b ;
    center_val = h->counts[b];
    if (center_val > max_count)
      max_count = center_val ;
  }

  if (!max_count)
    return(-1) ;

  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  if (b0 < 0)
    b0 = 0 ;
  if ((b1 < 0) || (b1 >= nbins))
    b1 = nbins - 1 ;

  /*
    check to see if the value at b is smaller than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = b0 ; b <= b1 ; b++)
  {
    center_val = h->counts[b] ;
    valley = 1 ;
    for (bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins)
        continue ;
      other_val = h->counts[bw] ;
      if (other_val < center_val)
        valley = 0 ;  /* something is smaller than current - not minimum */
    }
    if (valley)
      return(b) ;
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
/* only a peak if it is at least MIN_STD intensity units away
   from the mean in the wsize neighborhood.
*/
#define MIN_STD   1.9

int
HISTOfindLastPeakInRegion(HISTOGRAM *h, int wsize, float min_pct, int I0,
                          int I1)
{
  int    peak, b, bw, nbins, whalf, b0, b1 ;
  float  mean_count, min_count, max_count, other_val, center_val, total ;


  // find max, and find bins corresponding to I0 and I1
  b0 = b1 = -1 ;
  for (max_count = b = 0 ; b < h->nbins ; b++)
  {
    if (h->bins[b] >= I0 && b0 == -1)
      b0 = b > 0 ? b-1 : b ;
    if (h->bins[b] >= I1 && b1 == -1)
      b1 = b < h->nbins-1 ? b+1 : b ;
    center_val = h->counts[b];
    if (center_val > max_count)
      max_count = center_val ;
  }
  if (b1 == -1)
    b1 = h->nbins-1 ;

  if (!max_count)
    return(-1) ;

  min_count = min_pct * max_count ;
  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  /*
    check to see if the value at b is bigger than anything else within
    a whalf x whalf window on either side.
  */
  for (b = b1 ; b >= b0 ; b--)
  {
    center_val = h->counts[b] ;
    peak = 1 ;
    for (total = 0, bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins)
        continue ;
      other_val = h->counts[bw] ;
      total += other_val ;
      if (other_val > center_val)
      {
        peak = 0 ;
        break ;
      }
    }
    /* if average height in peak is greater than min_count accept it */
    if (peak)
    {
      if ((float)total/(float)wsize >= min_count)
      {
        mean_count = (float)total / (float)wsize ;
        if (((float)center_val-mean_count) >= MIN_STD)
          return(b) ;
      }
    }
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindFirstPeakInRegion(HISTOGRAM *h, 
                           int wsize, float min_pct,
                           int b0, int b1)
{
  int   peak, b, bw, nbins, whalf ;
  float center_val, max_count, other_val, min_count, total ;

  for (max_count = b = 0 ; b < h->nbins ; b++)
  {
    center_val = h->counts[b];
    if (center_val > max_count)
      max_count = center_val ;
  }

  if (!max_count)
    return(-1) ;

  min_count = min_pct * max_count ;
  whalf = (wsize-1)/2 ;
  nbins = h->nbins ;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = b0 ; b <= b1 ; b++)
  {
    center_val = h->counts[b] ;
    if (center_val <= min_count)
      continue ;
    peak = 1 ;
    for (total = 0, bw = b-whalf ; bw <= b+whalf ; bw++)
    {
      if (bw < 0 || bw >= nbins )
        continue ;
      other_val = h->counts[bw] ;
      if (other_val > center_val)
        peak = 0 ;
    }
    /* if average height in peak is greater than min_count accept it */
    if (peak && ((float)total/(float)wsize >= min_count))
      return(b) ;
  }

  return(-1) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOfindHighestPeakInRegion(HISTOGRAM *h, int b0, int b1)
{
  int  b, nbins, max_count_bin ;
  float val, max_count ;

  nbins = h->nbins ;

  if (b0 < 0)
    b0 = 0 ;
  if (b1 >= h->nbins)
    b1 = h->nbins-1 ;
  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  max_count = 0 ;
  max_count_bin = -1 ;
  for (b = b0 ; b <= b1 ; b++)
  {
    val = h->counts[b] ;
    if (val > max_count)
    {
      max_count_bin = b ;
      max_count= val ;
    }
  }

  return(max_count_bin) ;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
HISTOplot(HISTOGRAM *histo, char *fname)
{
  FILE *fp ;
  int  bin_no, bmin, bmax ;

  fp = fopen(fname, "w") ;
  if (fp == NULL)
    return(ERROR_NOFILE);

  for (bmin = 0 ; bmin < histo->nbins ; bmin++)
    if (histo->counts[bmin] > 0)
      break ;
  for (bmax = histo->nbins-1 ; bmax > bmin ; bmax--)
    if (histo->counts[bmax] > 0)
      break ;

  for (bin_no = bmin ; bin_no <= bmax ; bin_no++)
    fprintf(fp, "%f  %f\n", histo->bins[bin_no], histo->counts[bin_no]) ;
  fclose(fp) ;

  return(NO_ERROR) ;
}


/* at least one value on each side which is below the central value */
/*
   the standard deviation of the peak around it's mean must be at least
   two or above.
*/
int
HISTOcountPeaksInRegion(HISTOGRAM *h, int wsize, float min_pct,
                        int *peaks, int max_peaks, int b0, int b1)
{
  int npeaks = 0, peak_index ;

  do
  {
    peak_index = HISTOfindLastPeakInRegion(h, wsize, min_pct, b0, b1);
    peaks[npeaks++] = peak_index ;
    b1 = peak_index - wsize/2 ;
  }
  while (peak_index >= 0 && npeaks < max_peaks) ;

  return(npeaks) ;
}


int
HISTOaddSample(HISTOGRAM *histo, float val, float bmin, float bmax)
{
  int    bin_no ;
  float  bin_size ;

  bin_size = (bmax - bmin) / ((float)histo->nbins-1) ;
  bin_no = nint((val - bmin) / bin_size) ;
  histo->counts[bin_no]++ ;

  return(bin_no) ;
}


int
HISTOfindPreviousValley(HISTOGRAM *h, int b0)
{
  int    b ;
  float  prev_val, val ;

  prev_val = h->counts[b0] ;
  if (b0 <= 0)
    return(b0) ;
  for (b = b0-1 ; b >= 0 ; b--)
  {
    val = h->counts[b] ;
    if (val > prev_val)
      return(b) ;
    prev_val = val ;
  }

  return(-1) ;
}


int
HISTOfindNextValley(HISTOGRAM *h, int b0)
{
  int  b ;
  float prev_val, val ;

  prev_val = h->counts[b0] ;
  if (b0 >= h->nbins)
    return(b0) ;
  for (b = b0+1 ; b < h->nbins ; b++)
  {
    val = h->counts[b] ;
    if (val > prev_val)
      return(b) ;
    prev_val = val ;
  }

  return(-1) ;
}


int
HISTOfindNextPeak(HISTOGRAM *h, int b0, int whalf)
{
  int   b, bk, peak = 0 ;
  float val ;

  if (b0 > (h->nbins-2))
    return(b0) ;
  for (b = b0+1 ; b <= h->nbins-whalf ; b++)
  {
    val = h->counts[b] ;
    for (peak = 1, bk = b-whalf ; bk <= b+whalf ; bk++)
      if (h->counts[bk] > val)
      {
        peak = 0 ;
        break ;
      }
    if (peak)
      break ;
  }

  return(peak ? b : -1) ;
}

int
HISTOfindPreviousPeak(HISTOGRAM *h, int b0, int whalf)
{
  int   b, bk, peak = 0 ;
  float val ;

  if (b0 > (h->nbins-2))
    return(b0) ;
  for (b = b0-1 ; b >= whalf ; b--)
  {
    val = h->counts[b] ;
    for (peak = 1, bk = b-whalf ; bk <= b+whalf ; bk++)
      if (h->counts[bk] > val)
      {
        peak = 0 ;
        break ;
      }
    if (peak)
      break ;
  }

  return(peak ? b : -1) ;
}

int
HISTOfindStartOfPeak(HISTOGRAM *h, int b0, float pct_peak)
{
  int  b, b1 ;
  float val, thresh ;

  b1 = HISTOfindPreviousValley(h, b0) ;

  thresh = h->counts[b0]*pct_peak ;
  if (b0 <= 0)
    return(b0) ;
  for (b = b0-1 ; b >= b1 ; b--)
  {
    val = h->counts[b] ;
    if (val < thresh)
      return(b) ;
  }

  return(b) ;
}


int
HISTOfindEndOfPeak(HISTOGRAM *h, int b0, float pct_peak)
{
  int  b, b1 ;
  float val, thresh ;

  b1 = HISTOfindNextValley(h, b0) ;

  thresh = h->counts[b0]*pct_peak ;
  if (b0 >= h->nbins)
    return(b0) ;
  for (b = b0+1 ; b <= b1 ; b++)
  {
    val = h->counts[b] ;
    if (val < thresh)
      return(b) ;
  }

  return(b) ;
}

int
HISTOfindCurrentPeak(HISTOGRAM *histo, int b0, int wsize, float min_pct)
{
  int  b, whalf, bw, peak, nbins = histo->nbins ;
  float next_count, prev_count, other_val, center_val ;
  float max_count, min_count ;

  peak = HISTOfindHighestPeakInRegion(histo, 0, histo->nbins) ;
  if (peak < 0)
    return(-1) ;
  max_count = histo->counts[peak] ;
  min_count = min_pct * max_count ;

  whalf = (wsize-1)/2 ;
  for (next_count = prev_count = 0, bw = b0-whalf ; bw <= b0+whalf ; bw++)
  {
    if (bw < 0)
      continue ;
    if (bw >= nbins)
      continue ;
    if (bw < b0)
      prev_count += histo->counts[bw] ;
    else if (bw > b0)
      next_count += histo->counts[bw] ;
  }

  if (next_count > prev_count)  /* search forwards */
  {
    for (b = b0-whalf ; b < histo->nbins ; b++)
    {
      if (b < 0)
        continue ;
      center_val = histo->counts[b] ;
      if (center_val < min_count)
        continue ;
      peak = 1 ;
      for (bw = b ; bw <= b+whalf ; bw++)
      {
        if (bw < 0 || bw >= nbins)
          continue ;
        other_val = histo->counts[bw] ;
        if (other_val > center_val)
        {
          peak = 0 ;
          break ;
        }
      }
      if (peak)
      {
        int bv ;

        bv = HISTOfindNextValley(histo, 0) ;
#if 0
        if (bv >= 0 && center_val*min_pct < histo->counts[bv])
          continue ;
#endif
        return(b) ;
      }
    }
  }
  else   /* search backwards */
  {
    for (b = b0 ; b >= 0 ; b--)
    {
      if (b >= histo->nbins)
        continue ;
      center_val = histo->counts[b] ;
      if (center_val < min_count)
        continue ;
      peak = 1 ;
      for (bw = b-whalf ; bw <= b+whalf ; bw++)
      {
        if (bw < 0 || bw >= nbins)
          continue ;
        other_val = histo->counts[bw] ;
        if (other_val > center_val)
        {
          peak = 0 ;
          break ;
        }
      }
      if (peak)
      {
        int bv ;

        bv = HISTOfindPreviousValley(histo, 0) ;
#if 0
        if (bv >= 0 && center_val*min_pct < histo->counts[bv])
          continue ;
#endif
        return(b) ;
      }
    }
  }

  return(-1) ;
}


int
HISTOfillHoles(HISTO *h)
{
  int b ;

  for (b = 1 ; b < h->nbins-1 ; b++)
  {
    if (h->counts[b] == 0)
      h->counts[b] = (h->counts[b-1] + h->counts[b+1]) / 2 ;
  }
  if (h->counts[0] == 0)
    h->counts[0] = h->counts[1] ;
  if (h->counts[h->nbins-1] == 0)
    h->counts[h->nbins-1] = h->counts[h->nbins-2] ;

  return(NO_ERROR) ;
}


int
HISTOtotalInRegion(HISTO *h, int b0, int b1)
{
  int b, total ;

  for (total = 0, b = b0 ; b <= b1 ; b++)
  {
    total += h->counts[b] ;
  }

  return(total) ;
}


int
HISTOclearZeroBin(HISTOGRAM *h)
{
  int b ;

  if (h->bins[0] > h->bin_size)   // zero bin not in range
    return(NO_ERROR) ;
  for (b = 0 ; b < h->nbins ; b++)
  {
    if (h->bins[b] >= 0)
      break ;
  }
  h->counts[b] = 0 ;

  return(NO_ERROR) ;
}


int
HISTOfindBin(HISTOGRAM *h, float val)
{
  int b ;

  for (b = h->nbins-1 ; b > 0 ; b--)
    if (h->bins[b-1] < val)
      return(b) ;

  return(0) ;
}


HISTO *
HISTOclearBG(HISTOGRAM *hsrc, HISTOGRAM *hdst, int *pbg_end)
{
  int   b0, nv ;
  float min_count ;

  hdst = HISTOcopy(hsrc, hdst) ;
  min_count = hsrc->counts[0]*0.1 ;
  for (b0 = 1 ; b0 < hsrc->nbins ; b0++)
    if (hdst->counts[b0] < min_count)
      break ;

  /* ignore any valleys that are still part of bg noise, by waiting
   until amplitude has fallen below threshold */
  nv = HISTOfindNextValley(hsrc, b0 == 0 ? 0 : b0-1) ;
  printf("clearing  bg bins 0->%d (%d)\n", nv, b0) ;
  HISTOclearBins(hsrc, hdst, 0, nv) ;
  *pbg_end = nv ;
  return(hdst) ;
}

#if 0
static double
histoComputeLinearFitError(HISTOGRAM *h1, HISTOGRAM *h2, double a, double b)
{
  int    b1, b2, h2_done[256] ;
  double error, sse, c1, c2;

  if (h2->nbins > 256)
    ErrorExit(ERROR_UNSUPPORTED, 
              "histoComputeLinearFitError: only 256 bins allowed") ;
  memset(h2_done, 0, sizeof(h2_done)) ;
  for (sse = 0.0, b1 = 0 ; b1 < h1->nbins ; b1++)
  {
    b2 = nint(b1*a+b) ;
    if ((b2 < 0) || (b2 > h2->nbins-1))
      c2 = 0 ;
    else
    {
      c2 = h2->counts[b2] ;
      h2_done[b2] = 1 ;
    }
    c1 = h1->counts[b1] ;
    error = (c2 - c1) ;
    sse += error*error ;
  }

  // inverse map
  for (b2 = 0 ; b2 < h2->nbins ; b2++)
  {
    if (h2_done[b2])
      continue ;
    b1 = nint((b2-b)/a) ;
    if ((b1 < 0) || (b1 > h1->nbins-1))
      c1 = 0 ;
    else
      c1 = h1->counts[b1] ;
    c2 = h2->counts[b2] ;
    error = (c2 - c1) ;
    sse += error*error ;
  }

  return(sse) ;
}
#endif


static double
histoComputeLinearFitError(HISTOGRAM *h1, HISTOGRAM *h2, double a, double b)
{
  //this one uses correlation instead of SSE, and should be more robust!
  int    b1, b2 ;
  double error, sse, c1, c2;

  for (sse = 0.0, b1 = 0 ; b1 < h1->nbins ; b1++)
  {
    b2 = nint(b1*a+b) ;
    if ((b2 < 0) || (b2 > h2->nbins-1))
      c2 = 0 ;
    else
      c2 = h2->counts[b2] ;
    c1 = h1->counts[b1] ;
    error = c1*c2;
    sse -= error;
  }

  // printf("sse = %g\n", sse);
  // inverse map
  for (b2 = 0 ; b2 < h2->nbins ; b2++)
  {
    b1 = nint((b2-b)/a) ;
    if ((b1 < 0) || (b1 > h1->nbins-1))
      c1 = 0 ;
    else
      c1 = h1->counts[b1] ;
    c2 = h2->counts[b2] ;
    error = c1*c2;
    sse -= error;
  }
  // printf("a= %g, b= %g, sse = %g\n", a, b, sse);

  return(sse) ;
}


#define NSTEPS 100
int
HISTOfindLinearFit(HISTOGRAM *h1, HISTOGRAM *h2, double amin, double amax,
                   double bmin, double bmax, float *pa, float *pb)
{
  double   a, b, sse, min_sse, min_a, min_b, astep, bstep ;

  min_sse = histoComputeLinearFitError(h1, h2, 1.0, 0.0) ;
  min_a = 1.0 ;
  min_b = 0.0 ;
  astep = MIN(amin, (amax-amin)/NSTEPS) ;
  astep = MAX(astep, 0.01);
  bstep = (bmax-bmin)/NSTEPS ;
  bstep = MAX(bstep, 0.01);
  for (a = amin ; a <= amax ; a += astep)
    for (b = bmin ; b <= bmax ; b += bstep)
    {
      sse = histoComputeLinearFitError(h1, h2, a, b) ;
      if (sse < min_sse)
      {
        min_sse = sse ;
        min_a = a ;
        min_b = b ;
      }
    }

  *pa = min_a ;
  *pb = min_b ;

  return(NO_ERROR) ;
}


HISTO *
HISTOlinearScale(HISTOGRAM *hsrc, HISTOGRAM *hdst, float scale, float offset)
{
  int    b ;

  hdst = HISTOcopy(hsrc, hdst) ;

  for (b = 0 ; b < hdst->nbins ; b++)
    hdst->bins[b] = hdst->bins[b]*scale+offset ;

  return(hdst) ;
}


float
HISTOthreshSum(HISTOGRAM *h_mask, HISTOGRAM *h_src, float m_thresh)
{
  int   b ;
  float total ;

  for (total = 0.0, b = 0 ; b < h_mask->nbins ; b++)
  {
    if (h_mask->counts[b] > m_thresh)
      total += h_src->counts[b] ;
  }

  return(total) ;
}


HISTOGRAM *
HISTOmakePDF(HISTO *h_src, HISTO *h_dst)
{
  int    b ;
  float  total ;

  if (h_dst == NULL)
    h_dst = HISTOcopy(h_src, NULL) ;

  for (total = 0.0, b = 0 ; b < h_dst->nbins ; b++)
    total += h_dst->counts[b] ;

  if (total > 0)
    for (b = 0 ; b < h_dst->nbins ; b++)
      h_dst->counts[b]/=total ;

  return(h_dst) ;
}


/*-------------------------------------------------------------------
  HISTObins() - allocates histogram and assigns bin centers uniformly
  between min and max. The first bin is centered at min. The last
  bin is centered at max. This makes the bin size (max-min)/(nbins-1).
  -------------------------------------------------------------------*/
HISTO *HISTObins(int nbins, double min, double max)
{
  HISTO *h;
  int n;

  h = HISTOalloc(nbins);
  h->bin_size = (max-min)/(nbins-1);
  for (n=0; n < nbins; n++)  h->bins[n] = min + h->bin_size*n;

  return(h);
}


/*-------------------------------------------------------------------
  HISTOcount() - builds histogram based on samples. Must have already
  allocated hist and set bin centers.
  -------------------------------------------------------------------*/
int HISTOcount(HISTO *h, double *samples, int nsamples)
{
  int n,bin;

  for (n=0; n<nsamples; n++)
  {
    bin = HISTOvalToBin(h, samples[n]);
    h->counts[bin] ++;
  }

  return(0);
}


/*----------------------------------------------------------
  HISTOvalToBin() - returns the histogram bin number for
  the given value.
  ----------------------------------------------------------*/
int HISTOvalToBin(HISTO *h, double val)
{
  int bin, nthbin;
  double d,dmin;

  bin = 0;
  dmin = fabs(h->bins[0]-val);
  for (nthbin=0; nthbin < h->nbins; nthbin++)
  {
    d = fabs(h->bins[nthbin]-val);
    if (dmin > d)
    {
      dmin = d;
      bin = nthbin;
    }
  }

  return(bin);
}


int
HISTOvalToBinDirect(HISTOGRAM *histo, float val)
{
  int bin_no ;
  bin_no = nint((float)(val - histo->min) / (float)histo->bin_size) ;

  return(bin_no) ;
}


float
HISTOvalToCount(HISTOGRAM *histo, float val)
{
  int bin_no ;
  if (histo == NULL)
    return(0.0) ;
  bin_no = nint((float)(val - histo->min) / (float)histo->bin_size) ;
  if (bin_no < 0 || bin_no >= histo->nbins)
    return(0.0) ;

  return(histo->counts[bin_no]) ;
}
HISTOGRAM *
HISTOinit(HISTOGRAM *h, int nbins, double mn, double mx)
{
  int b ;

  if (h == NULL)
    h = HISTOalloc(nbins) ;
  else
    HISTOclear(h, h);

  h->min = mn ; h->max = mx ;
  h->bin_size = (mx - mn) / (nbins-1) ;
  for (b = 0 ; b < nbins ; b++)
    h->bins[b] = mn + h->bin_size*(float)b ;
  return(h) ;
}

int
HISTOwriteInto(HISTOGRAM *h, FILE *fp)
{
  int b ;

  fwriteInt(h->nbins, fp) ;
  fwriteFloat(h->bin_size, fp) ;
  fwriteFloat(h->min, fp) ;
  fwriteFloat(h->max, fp) ;
  for (b = 0 ; b < h->nbins ; b++)
    fwriteFloat(h->bins[b], fp) ;

  for (b = 0 ; b < h->nbins ; b++)
    fwriteFloat(h->counts[b], fp) ;
  return(NO_ERROR) ;
}

HISTOGRAM *
HISTOreadFrom(FILE *fp)
{
  int       b, nbins ;
  HISTOGRAM *h ;

  nbins = freadInt(fp) ;
  h = HISTOalloc(nbins) ;
  h->bin_size = freadFloat(fp) ;
  h->min = freadFloat(fp) ;
  h->max = freadFloat(fp) ;
  for (b = 0 ; b < h->nbins ; b++)
    h->bins[b] = freadFloat(fp) ;

  for (b = 0 ; b < h->nbins ; b++)
    h->counts[b] = freadFloat(fp) ;
  return(h) ;
}
