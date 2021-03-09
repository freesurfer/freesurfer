/**
 * @brief 1d histogram utilities.
 *
 * Utilities for computing and analyzing 1 dimensional histograms.
 */
/*
 * Original Author: Bruce Fischl
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#include <cmath>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "fio.h"
#include "histo.h"
#include "macros.h"
#include "proto.h"
#include "utils.h"

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfree(HISTOGRAM **phisto)
{
  HISTOGRAM *histo;

  histo = *phisto;
  *phisto = NULL;
  if (histo) {
    if (histo->bins) {
      free(histo->bins);
      histo->bins = NULL;
    }
    else
      DiagBreak();
    if (histo->counts) {
      free(histo->counts);
      histo->counts = NULL;
    }
    else
      DiagBreak();
    free(histo);
  }
  else
    DiagBreak();

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOdump(HISTOGRAM *histo, FILE *fp)
{
  int bin_no;

  if (!histo)
    fprintf(stderr, "NULL histogram");
  else {
    fprintf(fp, "nbins = %d\n", histo->nbins);
    for (bin_no = 0; bin_no < histo->nbins; bin_no++)
      if (histo->counts[bin_no])
        fprintf(fp, "bin[%d] = %2.1f = %2.2f\n", bin_no, histo->bins[bin_no], histo->counts[bin_no]);
  }
  return (NO_ERROR);
}

HISTOGRAM *HISTOrealloc(HISTOGRAM *histo, int nbins)
{
  if (histo == NULL) return (HISTOalloc(nbins));

  if (histo->bins) free(histo->bins);
  if (histo->counts) free(histo->counts);
  histo->bins = (float *)calloc(nbins, sizeof(float));
  histo->counts = (float *)calloc(nbins, sizeof(float));
  if (!histo->counts || !histo->bins)
    ErrorExit(ERROR_NOMEMORY, "HISTOrealloc(%d): could not allocate histogram", nbins);
  histo->max_bins = histo->nbins = nbins;

  return (histo);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOalloc(int nbins)
{
  HISTOGRAM *histo;

  histo = (HISTOGRAM *)calloc(1, sizeof(HISTOGRAM));
  if (!histo) ErrorExit(ERROR_NO_MEMORY, "HISTOalloc(%d): allocation failed", nbins);

  histo->bins = (float *)calloc(nbins, sizeof(float));
  histo->counts = (float *)calloc(nbins, sizeof(float));
  // fprintf(stderr, "histo->bins and ->counts allocated %d bins\n", nbins);
  if (!histo->counts || !histo->bins) ErrorExit(ERROR_NOMEMORY, "HISTOalloc(%d): could not allocate histogram", nbins);
  histo->max_bins = histo->nbins = nbins;

  return (histo);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOcrunch(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  int b, deleted;

  if (!histo_dst) histo_dst = HISTOalloc(histo_src->nbins);

  deleted = 0;
  for (b = 0; b < histo_src->nbins; b++) {
    if (!histo_src->counts[b])
      deleted++;
    else {
      histo_dst->counts[b - deleted] = histo_src->counts[b];
      histo_dst->bins[b - deleted] = histo_src->bins[b];
    }
  }

  histo_dst->nbins -= deleted;

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOcopy(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  if (!histo_dst) histo_dst = HISTOalloc(histo_src->nbins);
  histo_dst->nbins = histo_src->nbins;
  histo_dst->bin_size = histo_src->bin_size;
  // use memmove, not memcpy, as src and dst could overlap
  memmove(histo_dst->counts, histo_src->counts, sizeof(*histo_src->counts) * histo_src->nbins);
  memmove(histo_dst->bins, histo_src->bins, sizeof(*histo_src->bins) * histo_src->nbins);
  histo_dst->min = histo_src->min;
  histo_dst->max = histo_src->max;

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Treat histo_src as a function mapping input intensities
  to output intensities, and numerically invert it.
  ------------------------------------------------------*/
HISTOGRAM *HISTOinvert(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, int max_dst)
{
  int b, bdst, val, max_count;

  histo_dst = HISTOclear(histo_src, histo_dst);

  if (!max_dst) max_dst = 255;

  for (max_count = b = 0; b < histo_src->nbins; b++) {
    val = histo_src->counts[b];
    if (val > max_count) max_count = val;
  }

  for (b = 0; b < histo_src->nbins; b++) {
    val = histo_src->counts[b];
    bdst = nint((float)max_dst * (float)val / (float)max_count);
    if (bdst > max_dst) bdst = max_dst;
    if (!histo_dst->bins[bdst]) {
      histo_dst->counts[bdst] = b;
      histo_dst->bins[bdst] = 1;
    }
  }

  histo_dst->nbins = max_dst;
  /*
     fill in zeros in the inverse histogram - they correspond to
     flat regions in the forward map (i.e. multiple-valued).
  */
  for (val = b = 0; b < histo_dst->nbins; b++) {
    if (histo_dst->counts[b] > 0)
      val = histo_dst->counts[b];
    else
      histo_dst->counts[b] = val;
  }

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, int max_out)
{
  int max_count, b;
  float scale;

  if (!max_out) max_out = 255;

  if (!histo_dst) histo_dst = HISTOalloc(histo_src->nbins);

  for (max_count = b = 0; b < histo_src->nbins; b++)
    if (histo_src->counts[b] > max_count) max_count = histo_src->counts[b];

  scale = (float)max_out / (float)max_count;
  for (b = 0; b < histo_src->nbins; b++) histo_dst->counts[b] = (scale * (float)histo_src->counts[b]);

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOclear(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  if (!histo_dst) histo_dst = HISTOalloc(histo_src->nbins);

  memset(histo_dst->counts, 0, histo_dst->nbins * sizeof(*histo_dst->counts));
  memset(histo_dst->bins, 0, histo_dst->nbins * sizeof(*histo_dst->bins));

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOclearCounts(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  if (!histo_dst) {
    histo_dst = HISTOalloc(histo_src->nbins);
    memmove(histo_dst->bins, histo_src->bins, sizeof(*histo_src->bins));
  }

  memset(histo_dst->counts, 0, sizeof(histo_dst->counts[0]));

  return (histo_dst);
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
HISTOGRAM *HISTOcompose(HISTOGRAM *histo1, HISTOGRAM *histo2, HISTOGRAM *histo_dst)
{
  int b, val;

  if (!histo_dst) histo_dst = HISTOalloc(histo1->nbins);

  for (b = 0; b < histo1->nbins; b++) {
    val = histo1->counts[b];
    if (val >= histo2->nbins) val = histo2->nbins - 1;
    histo_dst->counts[b] = histo2->counts[val];
  }

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOfillZeros(HISTOGRAM *histo_src, HISTOGRAM *histo_dst)
{
  int b, val;

  histo_dst = HISTOcopy(histo_src, histo_dst);
  /*
     fill in zeros in the inverse histogram - they correspond to
     flat regions in the forward map (i.e. multiple-valued).
  */
  for (val = b = 0; b < histo_dst->nbins; b++) {
    if (histo_dst->counts[b] > 0)
      val = histo_dst->counts[b];
    else
      histo_dst->counts[b] = val;
  }

  return (histo_dst);
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

HISTOGRAM *HISTOcomposeInvert(HISTOGRAM *histo_fwd, HISTOGRAM *histo_inv, HISTOGRAM *histo_dst)
{
  int b, binv, val, max_fwd, max_inv, stretch;
  float ffwd, next = 0.0f, prev;

  if (!histo_dst)
    histo_dst = HISTOalloc(histo_fwd->nbins);
  else {
    HISTOclear(histo_dst, histo_dst);
    histo_dst->nbins = histo_fwd->nbins;
  }

  for (max_fwd = b = 0; b < histo_fwd->nbins; b++)
    if (histo_fwd->counts[b] > max_fwd) max_fwd = histo_fwd->counts[b];
  for (max_inv = b = 0; b < histo_inv->nbins; b++)
    if (histo_inv->counts[b] > max_inv) max_inv = histo_inv->counts[b];

  if (!max_inv || !max_fwd) return (histo_dst);

  for (binv = b = 0; b < histo_fwd->nbins; b++) {
    val = histo_fwd->counts[b];
    ffwd = (float)val / (float)max_fwd;

    for (stretch = 0; stretch < MAX_STRETCH && binv < histo_dst->nbins; stretch++, binv++) {
      next = (float)histo_inv->counts[binv] / (float)max_inv;
      if (next >= ffwd) break;
    }
    if (binv > 0) /* see whether bin to left or to right is closer */
    {
      prev = (float)histo_inv->counts[binv - 1] / (float)max_inv;
      if (fabs(next - ffwd) > abs(prev - ffwd)) /* 'prev' bin is closer */
        binv--;
    }
    histo_dst->counts[b] = binv++;
  }

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOadd(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst)
{
  int b;
  float *pc1, *pc2, *pcdst;

  if (!histo_dst) {
    histo_dst = HISTOalloc(h1->nbins);
    for (b = 0; b < h1->nbins; b++) histo_dst->bins[b] = h1->bins[b];
    histo_dst->min = h1->min;
    histo_dst->max = h1->max;
    histo_dst->bin_size = h1->bin_size;
    if (h2 == NULL) h2 = histo_dst;
  }
  else
    histo_dst->nbins = h1->nbins;

  pc1 = &h1->counts[0];
  pc2 = &h2->counts[0];
  pcdst = &histo_dst->counts[0];
  for (b = 0; b < h1->nbins; b++) *pcdst++ = *pc1++ + *pc2++;

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOmul(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst)
{
  int b;
  float *pc1, *pc2, *pcdst;

  if (!histo_dst)
    histo_dst = HISTOalloc(h1->nbins);
  else
    histo_dst->nbins = h1->nbins;

  pc1 = &h1->counts[0];
  pc2 = &h2->counts[0];
  pcdst = &histo_dst->counts[0];
  for (b = 0; b < h1->nbins; b++) *pcdst++ = *pc1++ * *pc2++;

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOsubtract(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst)
{
  int b;
  float *pc1, *pc2, *pcdst;

  if (!histo_dst)
    histo_dst = HISTOalloc(h1->nbins);
  else
    histo_dst->nbins = h1->nbins;

  pc1 = &h1->counts[0];
  pc2 = &h2->counts[0];
  pcdst = &histo_dst->counts[0];
  for (b = 0; b < h1->nbins; b++) *pcdst++ = *pc1++ - *pc2++;

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *HISTOclearBins(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, int min_val, int max_val)
{
  int b;

  if (!histo_src) return (NULL);

  if (!histo_dst || histo_dst != histo_src) histo_dst = HISTOcopy(histo_src, histo_dst);

  for (b = 0; b < histo_dst->nbins; b++) {
    if (histo_dst->bins[b] >= min_val && histo_dst->bins[b] <= max_val) {
      histo_dst->counts[b] = 0;
    }
  }

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MAX_LEN 2000
HISTOGRAM *HISTOsmooth(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, float sigma)
{
  float norm, two_sigma, fx, k, kernel[MAX_LEN], total;
  int x, half, len, b, kx, b1, nbins;

  nbins = histo_src->nbins;
  if (!histo_dst) {
    histo_dst = HISTOcopy(histo_src, NULL);
    HISTOclearCounts(histo_dst, histo_dst);
  }
  else {
    if (histo_dst->nbins < histo_src->nbins) {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTOrealloc(histo_dst, nbins);
    }
    histo_dst->nbins = nbins;
    histo_dst->bin_size = histo_src->bin_size;
  }
  /* build the kernel in k */
  len = (int)nint(8.0f * sigma) + 1;
  if (ISEVEN(len)) /* ensure it's odd */
    len++;
  if (MAX_LEN && (MAX_LEN < len)) len = MAX_LEN;
  half = len / 2;

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

    kernel[x] = k;
    norm += k;
  }
  for (x = 0; x < len; x++) kernel[x] /= norm;

  for (b = 0; b < nbins; b++) {
    for (norm = 0.0, total = 0.0f, x = 0; x < len; x++) {
      kx = x - half;
      b1 = b + kx;
      if (b1 >= nbins || b1 < 0) continue;

      norm += kernel[x];
      total += kernel[x] * (float)histo_src->counts[b1];
    }
    histo_dst->counts[b] = total / norm;
    histo_dst->bins[b] = histo_src->bins[b];
  }

  return (histo_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindLastPeakRelative(HISTOGRAM *h, int wsize, float min_pct)
{
  int peak, b, bw, nbins, whalf, other_val;
  float max_count, min_count, center_val;

  peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins);
  if (peak < 0) return (-1);
  max_count = h->counts[peak];

  min_count = min_pct * max_count;
  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = nbins - 1; b >= 0; b--) {
    center_val = h->counts[b];
    if (center_val <= min_count) continue;
    peak = 1;
    for (bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      if (other_val > center_val) peak = 0;
    }
    if (peak) return (b);
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindFirstPeakRelative(HISTOGRAM *h, int wsize, float min_pct)
{
  int peak, b, bw, nbins, whalf, other_val;
  float max_count, min_count, center_val;

  peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins);
  if (peak < 0) return (-1);
  max_count = h->counts[peak];

  min_count = min_pct * max_count;
  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = 0; b < nbins; b++) {
    center_val = h->counts[b];
    if (center_val <= min_count) continue;
    peak = 1;
    for (bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      if (other_val > center_val) peak = 0;
    }
    if (peak) return (b);
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindLastPeak(HISTOGRAM *h, int wsize, float min_pct)
{
  int peak, b, bw, nbins, whalf;
  float min_count, max_count, center_val, other_val;

  for (max_count = b = 0; b < h->nbins; b++) {
    center_val = h->counts[b];
    if (center_val > max_count) max_count = center_val;
  }

  if (!max_count) return (-1);

  min_count = (min_pct * (float)max_count);
  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = nbins - 1; b >= 0; b--) {
    center_val = h->counts[b];
    if (center_val <= min_count) continue;
    peak = 1;
    for (bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      if (other_val > center_val) peak = 0;
    }
    if (peak) return (b);
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindFirstPeak(HISTOGRAM *h, int wsize, float min_pct)
{
  int peak, b, bw, nbins, whalf;
  float center_val, max_count, other_val, min_count;

  for (max_count = b = 0; b < h->nbins; b++) {
    center_val = h->counts[b];
    if (center_val > max_count) max_count = center_val;
  }

  if (!max_count) return (-1);

  min_count = min_pct * max_count;
  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = 0; b < nbins; b++) {
    center_val = h->counts[b];
    if (center_val <= min_count) continue;
    peak = 1;
    for (bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      if (other_val > center_val) peak = 0;
#if 0
      else if (center_val > other_val)
        peak++ ;
#endif
    }
    if (peak) return (b);
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindValley(HISTOGRAM *h, int wsize, int I0, int I1)
{
  int valley, b, bw, nbins, whalf, center_val, max_count, other_val, b0, b1;

  // find max, and find bins corresponding to I0 and I1
  b0 = b1 = -1;
  for (max_count = b = 0; b < h->nbins; b++) {
    if (h->bins[b] >= I0 && b0 == -1) b0 = b > 0 ? b - 1 : b;
    if (h->bins[b] >= I1 && b1 == -1) b1 = b < h->nbins - 1 ? b + 1 : b;
    center_val = h->counts[b];
    if (center_val > max_count) max_count = center_val;
  }

  if (!max_count) return (-1);

  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  if (b0 < 0) b0 = 0;
  if ((b1 < 0) || (b1 >= nbins)) b1 = nbins - 1;

  /*
    check to see if the value at b is smaller than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = b0; b <= b1; b++) {
    center_val = h->counts[b];
    valley = 1;
    for (bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      if (other_val < center_val) valley = 0; /* something is smaller than current - not minimum */
    }
    if (valley) return (b);
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
/* only a peak if it is at least MIN_STD intensity units away
   from the mean in the wsize neighborhood.
*/
#define MIN_STD 1.9

int HISTOfindLastPeakInRegion(HISTOGRAM *h, int wsize, float min_pct, int I0, int I1)
{
  int peak, b, bw, nbins, whalf, b0, b1;
  float mean_count, min_count, max_count, other_val, center_val, total;

  // find max, and find bins corresponding to I0 and I1
  b0 = b1 = -1;
  for (max_count = b = 0; b < h->nbins; b++) {
    if (h->bins[b] >= I0 && b0 == -1) b0 = b > 0 ? b - 1 : b;
    if (h->bins[b] >= I1 && b1 == -1) b1 = b < h->nbins - 1 ? b + 1 : b;
    center_val = h->counts[b];
    if (center_val > max_count) max_count = center_val;
  }
  if (b1 == -1) b1 = h->nbins - 1;

  if (!max_count) return (-1);

  min_count = min_pct * max_count;
  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  /*
    check to see if the value at b is bigger than anything else within
    a whalf x whalf window on either side.
  */
  for (b = b1; b >= b0; b--) {
    center_val = h->counts[b];
    peak = 1;
    for (total = 0, bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      total += other_val;
      if (other_val > center_val) {
        peak = 0;
        break;
      }
    }
    /* if average height in peak is greater than min_count accept it */
    if (peak) {
      if ((float)total / (float)wsize >= min_count) {
        mean_count = (float)total / (float)wsize;
        if (((float)center_val - mean_count) >= MIN_STD) return (b);
      }
    }
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindFirstPeakInRegion(HISTOGRAM *h, int wsize, float min_pct, int b0, int b1)
{
  int peak, b, bw, nbins, whalf;
  float center_val, max_count, other_val, min_count, total;

  for (max_count = b = 0; b < h->nbins; b++) {
    center_val = h->counts[b];
    if (center_val > max_count) max_count = center_val;
  }

  if (!max_count) return (-1);

  min_count = min_pct * max_count;
  whalf = (wsize - 1) / 2;
  nbins = h->nbins;

  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  for (b = b0; b <= b1; b++) {
    center_val = h->counts[b];
    if (center_val <= min_count) continue;
    peak = 1;
    for (total = 0, bw = b - whalf; bw <= b + whalf; bw++) {
      if (bw < 0 || bw >= nbins) continue;
      other_val = h->counts[bw];
      if (other_val > center_val) peak = 0;
    }
    /* if average height in peak is greater than min_count accept it */
    if (peak && ((float)total / (float)wsize >= min_count)) return (b);
  }

  return (-1);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOfindHighestPeakInRegion(HISTOGRAM *h, int b0, int b1)
{
  // int  nbins;
  int b, max_count_bin;
  float val, max_count;

  // nbins = h->nbins;

  if (b0 < 0) b0 = 0;
  if (b1 >= h->nbins) b1 = h->nbins - 1;
  /*
    check to see if the value at b is bigger than anything else within
    a whalfxwhalf window on either side.
  */
  max_count = 0;
  max_count_bin = -1;
  for (b = b0; b <= b1; b++) {
    val = h->counts[b];
    if (val > max_count) {
      max_count_bin = b;
      max_count = val;
    }
  }

  return (max_count_bin);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int HISTOplot(HISTOGRAM *histo, const char *fname)
{
  FILE *fp;
  int bin_no, bmin, bmax;

  fp = fopen(fname, "w");
  if (fp == NULL) return (ERROR_NOFILE);

  for (bmin = 0; bmin < histo->nbins; bmin++)
    if (histo->counts[bmin] > 0) break;
  for (bmax = histo->nbins - 1; bmax > bmin; bmax--)
    if (histo->counts[bmax] > 0) break;

  for (bin_no = bmin; bin_no <= bmax; bin_no++)
    fprintf(fp, "%f  %10.10f\n", histo->bins[bin_no], histo->counts[bin_no]);
  fclose(fp);

  return (NO_ERROR);
}

/* at least one value on each side which is below the central value */
/*
   the standard deviation of the peak around it's mean must be at least
   two or above.
*/
int HISTOcountPeaksInRegion(HISTOGRAM *h, int wsize, float min_pct, int *peaks, int max_peaks, int b0, int b1)
{
  int npeaks = 0, peak_index;

  do {
    peak_index = HISTOfindLastPeakInRegion(h, wsize, min_pct, b0, b1);
    peaks[npeaks++] = peak_index;
    b1 = peak_index - wsize / 2;
  } while (peak_index >= 0 && npeaks < max_peaks);

  return (npeaks);
}

int HISTOaddSample(HISTOGRAM *histo, float val, float bmin, float bmax)
{
  int bin_no;
  float bin_size;

  if (bmin >= bmax)  // not specified by caller
  {
    bmin = histo->min;
    bmax = histo->max;
    bin_size = histo->bin_size;
  }
  else
    bin_size = (bmax - bmin) / ((float)histo->nbins - 1);

  bin_no = nint((val - bmin) / bin_size);
  if (bin_no < histo->nbins && bin_no >= 0)
  {
    double dp, dn, dbin ;
    int    bp, bn ;
    dbin = ((val - bmin) / bin_size) ;
    dp = dbin-floor(dbin) ;
    dn = ceil(dbin)-dbin ;
    bp = (int)floor(dbin) ;
    bn = (int)ceil(dbin) ;
    if (bp < 0)
      bp = 0 ;
    if (bn >= histo->nbins)
      bn = histo->nbins-1 ;
    histo->counts[bn] += dp ;
    histo->counts[bp]  += dn ;
  }
  else {
    if (bin_no < bmin)
      histo->counts[0]++;
    else if (bin_no >= histo->nbins)
      histo->counts[histo->nbins - 1]++;
  }

  return (bin_no);
}

int HISTOaddFractionalSample(HISTOGRAM *histo, float val, float bmin, float bmax, float frac)
{
  int bin_no;
  float bin_size;

  if (bmin >= bmax)  // not specified by caller
  {
    bmin = histo->min;
    bmax = histo->max;
    bin_size = histo->bin_size;
  }
  else
    bin_size = (bmax - bmin) / ((float)histo->nbins - 1);

  bin_no = nint((val - bmin) / bin_size);
  if (bin_no < histo->nbins && bin_no >= 0)
    histo->counts[bin_no] += frac;
  else
    DiagBreak();

  return (bin_no);
}

int HISTOfindPreviousValley(HISTOGRAM *h, int b0)
{
  int b;
  float prev_val, val;

  prev_val = h->counts[b0];
  if (b0 <= 0) return (b0);
  for (b = b0 - 1; b >= 0; b--) {
    val = h->counts[b];
    if (val > prev_val) return (b);
    prev_val = val;
  }

  return (-1);
}

int HISTOfindNextValley(HISTOGRAM *h, int b0)
{
  int b;
  float prev_val, val;

  prev_val = h->counts[b0];
  if (b0 >= h->nbins) return (b0);
  for (b = b0 + 1; b < h->nbins; b++) {
    val = h->counts[b];
    if (val > prev_val) return (b);
    prev_val = val;
  }

  return (-1);
}

int HISTOfindNextPeak(HISTOGRAM *h, int b0, int whalf)
{
  int b, bk, peak = 0;
  float val;

  if (b0 > (h->nbins - 2)) return (b0);
  for (b = b0 + 1; b <= h->nbins - whalf; b++) {
    val = h->counts[b];
    for (peak = 1, bk = b - whalf; bk <= b + whalf; bk++)
      if (h->counts[bk] > val) {
        peak = 0;
        break;
      }
    if (peak) break;
  }

  return (peak ? b : -1);
}

int HISTOfindPreviousPeak(HISTOGRAM *h, int b0, int whalf)
{
  int b, bk, peak = 0;
  float val;

  if (b0 > (h->nbins - 2)) return (b0);
  for (b = b0 - 1; b >= whalf; b--) {
    val = h->counts[b];
    for (peak = 1, bk = b - whalf; bk <= b + whalf; bk++)
      if (h->counts[bk] > val) {
        peak = 0;
        break;
      }
    if (peak) break;
  }

  return (peak ? b : -1);
}

int HISTOfindStartOfPeak(HISTOGRAM *h, int b0, float pct_peak)
{
  int b, b1;
  float val, thresh;

  b1 = HISTOfindPreviousValley(h, b0);

  thresh = h->counts[b0] * pct_peak;
  if (b0 <= 0) return (b0);
  for (b = b0 - 1; b >= b1; b--) {
    val = h->counts[b];
    if (val < thresh) return (b);
  }

  return (b);
}

int HISTOfindEndOfPeak(HISTOGRAM *h, int b0, float pct_peak)
{
  int b, b1;
  float val, thresh;

  b1 = HISTOfindNextValley(h, b0);

  thresh = h->counts[b0] * pct_peak;
  if (b0 >= h->nbins) return (b0);
  for (b = b0 + 1; b <= b1; b++) {
    val = h->counts[b];
    if (val < thresh) return (b);
  }

  if (h->nbins >= b) b = h->nbins - 1;
  return (b);
}

int HISTOfindCurrentPeak(HISTOGRAM *histo, int b0, int wsize, float min_pct)
{
  int b, whalf, bw, peak, nbins = histo->nbins;
  float next_count, prev_count, other_val, center_val;
  float max_count, min_count;

  peak = HISTOfindHighestPeakInRegion(histo, 0, histo->nbins);
  if (peak < 0) return (-1);
  max_count = histo->counts[peak];
  min_count = min_pct * max_count;

  whalf = (wsize - 1) / 2;
  for (next_count = prev_count = 0, bw = b0 - whalf; bw <= b0 + whalf; bw++) {
    if (bw < 0) continue;
    if (bw >= nbins) continue;
    if (bw < b0)
      prev_count += histo->counts[bw];
    else if (bw > b0)
      next_count += histo->counts[bw];
  }

  if (next_count > prev_count) /* search forwards */
  {
    for (b = b0 - whalf; b < histo->nbins; b++) {
      if (b < 0) continue;
      center_val = histo->counts[b];
      if (center_val < min_count) continue;
      peak = 1;
      for (bw = b; bw <= b + whalf; bw++) {
        if (bw < 0 || bw >= nbins) continue;
        other_val = histo->counts[bw];
        if (other_val > center_val) {
          peak = 0;
          break;
        }
      }
      if (peak) {
#if 1
        HISTOfindNextValley(histo, 0);
#else
        int bv = HISTOfindNextValley(histo, 0);
        if (bv >= 0 && center_val * min_pct < histo->counts[bv]) continue;
#endif
        return (b);
      }
    }
  }
  else /* search backwards */
  {
    for (b = b0; b >= 0; b--) {
      if (b >= histo->nbins) continue;
      center_val = histo->counts[b];
      if (center_val < min_count) continue;
      peak = 1;
      for (bw = b - whalf; bw <= b + whalf; bw++) {
        if (bw < 0 || bw >= nbins) continue;
        other_val = histo->counts[bw];
        if (other_val > center_val) {
          peak = 0;
          break;
        }
      }
      if (peak) {
#if 1
        HISTOfindPreviousValley(histo, 0);
#else
        int bv = HISTOfindPreviousValley(histo, 0);
        if (bv >= 0 && center_val * min_pct < histo->counts[bv]) continue;
#endif
        return (b);
      }
    }
  }

  return (-1);
}

int HISTOfillHoles(HISTO *h)
{
  int b, min_b, max_b, *filled, niter;
  double max_change, new_val, change;

  filled = (int *)calloc(h->nbins, sizeof(int));
  if (filled == NULL) ErrorExit(ERROR_NOMEMORY, "HISTOfillHoles: couldn't allocate %d-len buffer\n", h->nbins);

  min_b = h->nbins;
  max_b = 0;
  for (b = 1; b < h->nbins - 1; b++) {
    if (h->counts[b] > 0) {
      filled[b] = 1;
      if (b > max_b) max_b = b;
      if (b < min_b) min_b = b;
    }
  }
  niter = 0;
  do {
    max_change = 0;
    for (b = min_b; b <= max_b; b++) {
      if (filled[b] == 0)  // was a hole
      {
        new_val = (h->counts[b - 1] + h->counts[b + 1]) / 2;
        if (h->counts[b] == 0)
          change = 1;
        else
          change = fabs(new_val - h->counts[b]) / fabs(h->counts[b]);
        h->counts[b] = new_val;
        if (change > max_change) max_change = change;
      }
    }
  } while (max_change > 0.01 && niter++ < 100);

  if (h->counts[0] == 0) h->counts[0] = h->counts[1];
  if (h->counts[h->nbins - 1] == 0) h->counts[h->nbins - 1] = h->counts[h->nbins - 2];

  free(filled);
  return (NO_ERROR);
}

float HISTOtotalInRegion(HISTO *h, int b0, int b1)
{
  int b;
  float total = 0.0;

  for (b = b0; b <= b1; b++) {
    total += h->counts[b];
  }

  return (total);
}

float HISTOtotal(HISTO *h) { return HISTOtotalInRegion(h, 0, h->nbins - 1); }

int HISTOclearZeroBin(HISTOGRAM *h)
{
  int b;

  if (h->bins[0] > h->bin_size)  // zero bin not in range
    return (NO_ERROR);
  for (b = 0; b < h->nbins; b++) {
    if (h->bins[b] >= 0) break;
  }
  h->counts[b] = 0;

  return (NO_ERROR);
}

int HISTOfindBinWithCount(HISTOGRAM *h, float val)
{
  int b, min_b;
  double min_dist, dist;

  min_b = 0;
  min_dist = fabs(h->counts[0] - val);
  for (b = 1; b < h->nbins; b++) {
    dist = fabs(h->counts[b] - val);
    if (dist < min_dist) {
      min_dist = dist;
      min_b = b;
    }
  }

  return (min_b);
}

int HISTOfindBin(HISTOGRAM *h, float val)
{
  int b;

  if (h->bin_size == 1) return ((int)val - h->bins[0]);

  for (b = h->nbins - 1; b > 0; b--)
    if (h->bins[b - 1] < val) return (b);

  return (0);
}
double HISTOgetCount(HISTOGRAM *h, float bin_val)
{
  int b;

  b = HISTOfindBin(h, bin_val);
  if (b < 0 || b >= h->nbins) return (0);
  return (h->counts[b]);
}

HISTO *HISTOclearBG(HISTOGRAM *hsrc, HISTOGRAM *hdst, int *pbg_end)
{
  int b0, nv;
  float min_count;

  b0 = HISTOfindNextPeak(hsrc, 0, 4);
  hdst = HISTOcopy(hsrc, hdst);
  min_count = hsrc->counts[b0] * 0.1;
  for (; b0 < hsrc->nbins; b0++)
    if (hdst->counts[b0] < min_count) break;

  /* ignore any valleys that are still part of bg noise, by waiting
   until amplitude has fallen below threshold */
  nv = HISTOfindNextValley(hsrc, b0 == 0 ? 0 : b0 - 1);
  printf("clearing  bg bins 0->%d (%d)\n", nv, b0);
  HISTOclearBins(hsrc, hdst, 0, nv);
  *pbg_end = nv;
  return (hdst);
}

double HISTOcorrelate(HISTOGRAM *h1, HISTOGRAM *h2)
{
  int b1, b2, h2_done[256];
  double correlation, c1, c2, norm1, norm2;

  if (h2->nbins > 256) ErrorExit(ERROR_UNSUPPORTED, "histoComputeLinearFitCorrelation: only 256 bins allowed");
  memset(h2_done, 0, sizeof(h2_done));
  norm1 = norm2 = 0.0;
  for (correlation = 0, b1 = 0; b1 < h1->nbins; b1++) {
    b2 = HISTOfindBin(h2, h1->bins[b1]);
    if ((b2 < 0) || (b2 > h2->nbins - 1))
      c2 = 0;
    else {
      c2 = h2->counts[b2];
      h2_done[b2] = 1;
    }
    c1 = h1->counts[b1];
    correlation += c1 * c2;
    norm1 += c1 * c1;
    norm2 += c2 * c2;
  }

  if (FZERO(norm1)) norm1 = 1.0;
  if (FZERO(norm2)) norm2 = 1.0;
  correlation /= (sqrt(norm1) * sqrt(norm2));
  return (correlation);
}

static double histoComputeLinearFitCorrelation(HISTOGRAM *h1, HISTOGRAM *h2, double a, double b)
{
  // this one uses correlation instead of SSE, and should be more robust!
  int b1, b2;
  double correlation, c1, c2;

  for (correlation = 0.0, b1 = 0; b1 < h1->nbins; b1++) {
    b2 = nint(b1 * a + b);
    if ((b2 < 0) || (b2 > h2->nbins - 1))
      c2 = 0;
    else
      c2 = h2->counts[b2];
    c1 = h1->counts[b1];
    correlation += c1 * c2;
  }

  // printf("correlation = %g\n", correlation);
  // inverse map
  for (b2 = 0; b2 < h2->nbins; b2++) {
    b1 = nint((b2 - b) / a);
    if ((b1 < 0) || (b1 > h1->nbins - 1))
      c1 = 0;
    else
      c1 = h1->counts[b1];
    c2 = h2->counts[b2];
    correlation += c1 * c2;
  }
  // printf("a= %g, b= %g, correlation = %g\n", a, b, correlation);

  return (correlation);
}

#define NSTEPS 1000
double HISTOfindLinearFit(
    HISTOGRAM *h1, HISTOGRAM *h2, double amin, double amax, double bmin, double bmax, float *pa, float *pb)
{
  double a, b, correlation, best_correlation, best_a, best_b, astep, bstep;

  best_correlation = histoComputeLinearFitCorrelation(h1, h2, 1.0, 0.0);
  best_a = 1.0;
  best_b = 0.0;
  astep = MIN(amin, (amax - amin) / NSTEPS);
  astep = MAX(astep, 0.01);
  bstep = (bmax - bmin) / NSTEPS;
  bstep = MAX(bstep, 0.01);
  for (a = amin; a <= amax; a += astep)
    for (b = bmin; b <= bmax; b += bstep) {
      correlation = histoComputeLinearFitCorrelation(h1, h2, a, b);
      if (correlation > best_correlation) {
        best_correlation = correlation;
        best_a = a;
        best_b = b;
      }
    }

  *pa = best_a;
  *pb = best_b;

  return (best_correlation);
}

HISTO *HISTOlinearScale(HISTOGRAM *hsrc, HISTOGRAM *hdst, float scale, float offset)
{
  int b;

  hdst = HISTOcopy(hsrc, hdst);

  for (b = 0; b < hdst->nbins; b++) hdst->bins[b] = hdst->bins[b] * scale + offset;

  hdst->min = hsrc->min * scale + offset;
  hdst->max = hsrc->max * scale + offset;
  hdst->bin_size = hsrc->bin_size * scale;
  return (hdst);
}

float HISTOthreshSum(HISTOGRAM *h_mask, HISTOGRAM *h_src, float m_thresh)
{
  int b;
  float total;

  for (total = 0.0, b = 0; b < h_mask->nbins; b++) {
    if (h_mask->counts[b] > m_thresh) total += h_src->counts[b];
  }

  return (total);
}

HISTOGRAM *HISTOmakePDF(HISTO *h_src, HISTO *h_dst)
{
  int b;
  float total;

  if (h_dst == NULL) h_dst = HISTOcopy(h_src, NULL);

  for (total = 0.0, b = 0; b < h_dst->nbins; b++) total += h_dst->counts[b];

  if (total > 0)
    for (b = 0; b < h_dst->nbins; b++) {
      h_dst->counts[b] /= total;
      if (DZERO(h_dst->counts[b])) h_dst->counts[b] = 1.0 / (10 * total);
    }

  return (h_dst);
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
  h->bin_size = (max - min) / (nbins - 1);
  for (n = 0; n < nbins; n++) h->bins[n] = min + h->bin_size * n;
  h->min = min;
  h->max = max;

  return (h);
}

/*-------------------------------------------------------------------
  HISTOcount() - builds histogram based on samples. Must have already
  allocated hist and set bin centers.
  -------------------------------------------------------------------*/
int HISTOcount(HISTO *h, double *samples, int nsamples)
{
  int n, bin;

  for (n = 0; n < nsamples; n++) {
    bin = HISTOvalToBin(h, samples[n]);
    h->counts[bin]++;
  }

  return (0);
}

/*----------------------------------------------------------
  HISTOvalToBin() - returns the histogram bin number for
  the given value.
  ----------------------------------------------------------*/
int HISTOvalToBin(HISTO *h, double val)
{
  int bin, nthbin;
  double d, dmin;

  bin = 0;
  dmin = fabs(h->bins[0] - val);
  for (nthbin = 0; nthbin < h->nbins; nthbin++) {
    d = fabs(h->bins[nthbin] - val);
    if (dmin > d) {
      dmin = d;
      bin = nthbin;
    }
  }

  return (bin);
}

int HISTOvalToBinDirect(HISTOGRAM *histo, float val)
{
  int bin_no;
  bin_no = nint((float)(val - histo->min) / (float)histo->bin_size);

  return (bin_no);
}

float HISTOvalToCount(HISTOGRAM *histo, float val)
{
  int bin_no;
  if (histo == NULL) return (0.0);
  bin_no = nint((float)(val - histo->min) / (float)histo->bin_size);
  if (bin_no < 0 || bin_no >= histo->nbins) return (0.0);

  return (histo->counts[bin_no]);
}
HISTOGRAM *HISTOinit(HISTOGRAM *h, int nbins, double mn, double mx)
{
  int b;

  if (h == NULL)
    h = HISTOalloc(nbins);
  else
    HISTOclear(h, h);

  if (nbins > h->max_bins)
    ErrorExit(ERROR_BADPARM, "HISTOinit: nbins %d bigger than max bins %d\n", nbins, h->max_bins);
  if (mx <= mn) {
    mx = (double)(h->max_bins - 1);
    mn = 0;
  }
  h->min = mn;
  h->max = mx;
  h->bin_size = (mx - mn) / (nbins - 1);
  if (h->bin_size <= 0) h->bin_size = 1;
  for (b = 0; b < nbins; b++) h->bins[b] = mn + h->bin_size * (float)b;
  return (h);
}

int HISTOwriteInto(HISTOGRAM *h, FILE *fp)
{
  int b;

  fwriteInt(h->nbins, fp);
  fwriteFloat(h->bin_size, fp);
  fwriteFloat(h->min, fp);
  fwriteFloat(h->max, fp);
  for (b = 0; b < h->nbins; b++) fwriteFloat(h->bins[b], fp);

  for (b = 0; b < h->nbins; b++) fwriteFloat(h->counts[b], fp);
  return (NO_ERROR);
}

HISTOGRAM *HISTOreadFrom(FILE *fp)
{
  int b, nbins;
  HISTOGRAM *h;

  nbins = freadInt(fp);
  h = HISTOalloc(nbins);
  h->bin_size = freadFloat(fp);
  h->min = freadFloat(fp);
  h->max = freadFloat(fp);
  for (b = 0; b < h->nbins; b++) h->bins[b] = freadFloat(fp);

  for (b = 0; b < h->nbins; b++) h->counts[b] = freadFloat(fp);
  return (h);
}
double HISTOfindMedian(HISTOGRAM *h)
{
  double median, total, total2;
  int b;

  for (total = 0.0, b = 0; b < h->nbins; b++) total += h->counts[b];
  median = 0.0;
  for (total2 = 0.0, b = 0; b < h->nbins; b++) {
    if (total2 < total / 2 && total2 + h->counts[b] > total / 2) {
      double d1, d2;
      d1 = total / 2 - total2;
      d2 = total2 + h->counts[b] - total / 2;
      if (d1 > d2 || b == 0)
        median = h->bins[b];
      else
        median = h->bins[b - 1];
      break;
    }
    total2 += h->counts[b];
  }
  return (median);
}

HISTOGRAM *HISTOmakeReverseCDF(HISTOGRAM *hsrc, HISTOGRAM *hdst)
{
  int b, b2;
  double total = 0;
  HISTOGRAM *habs;

  habs = HISTOabs(hsrc, NULL);

  if (hdst == NULL) hdst = HISTOcopy(hsrc, NULL);
  for (b = hdst->nbins - 1; b >= 0; b--) {
    if (hdst->bins[b] >= 0) {
      for (total = 0.0, b2 = habs->nbins - 1; b2 >= 0; b2--)
        if (fabs(habs->bins[b2]) >= hdst->bins[b]) total += habs->counts[b2];
      hdst->counts[b] = total;
    }
    else
      hdst->counts[b] = 0;
  }

  for (b = 0; b < hdst->nbins; b++) hdst->counts[b] /= total;

  HISTOfree(&habs);
  return (hdst);
}
HISTOGRAM *HISTOmakeCDF(HISTOGRAM *hsrc, HISTOGRAM *hdst)
{
  int b, b2;
  double total = 0;
  HISTOGRAM *habs;

  habs = HISTOabs(hsrc, NULL);

  if (hdst == NULL) hdst = HISTOcopy(hsrc, NULL);
  for (b = 0; b < hdst->nbins; b++) {
    if (hdst->bins[b] >= 0) {
      for (total = 0.0, b2 = 0; b2 < habs->nbins; b2++)
        if (fabs(habs->bins[b2]) <= hdst->bins[b]) total += habs->counts[b2];
      hdst->counts[b] = total;
    }
    else
      hdst->counts[b] = 0;
  }

  for (b = 0; b < hdst->nbins; b++) hdst->counts[b] /= total;

  HISTOfree(&habs);
  return (hdst);
}

int HISTOrobustGaussianFit(HISTOGRAM *h, double max_percentile, double *poffset, double *psigma)
{
  int peak, bin, min_bin, max_bin, zbin, n;
  HISTOGRAM *habs, *hcdf, *hz, *hs;
  double thresh, sigma, delta_sigma, max_sigma, predicted_val, val, sqrt_2pi, sse, best_sigma, best_sse, error, scale,
      mean;

  if (h->bin_size <= 0) h->bin_size = 1;
  sqrt_2pi = sqrt(2 * M_PI);
  hs = HISTOsmooth(h, NULL, 2.0 / h->bin_size);
  peak = HISTOfindHighestPeakInRegion(hs, 0, hs->nbins);
  mean = *poffset = hs->bins[peak];
  hz = HISTOlinearScale(h, NULL, 1.0, -mean);
  habs = HISTOabs(hz, NULL);
  hcdf = HISTOmakeCDF(hz, NULL);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(habs, "ha.plt");
    HISTOplot(hcdf, "hc.plt");
    HISTOplot(h, "h.plt");
    HISTOplot(hs, "hs.plt");
    HISTOplot(hz, "hz.plt");
  }

  bin = HISTOfindNextValley(hs, peak + 10 / hs->bin_size);  // real valley must not be right next to peak
  if (bin < 0) {
    bin = HISTOfindBinWithCount(hcdf, max_percentile);
    thresh = hcdf->bins[bin];  // threshold for computing variance
  }
  else  // found a valid next valley. Now look backwards to find bulk of distribution
  {
    while (h->counts[bin] < h->counts[peak] * max_percentile &&
           bin > peak + 1)  // don't fit way out on tails of distribution
      bin--;
    thresh = h->bins[bin] - h->bins[peak];
  }
  min_bin = HISTOfindBin(hz, -thresh);
  max_bin = HISTOfindBin(hz, thresh);
  if (hz->nbins == 0) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "HISTOrobustGaussianFit: 0 length cdf"));
  while (min_bin >= max_bin)  // shouldn't happen
  {
    min_bin = MAX(0, min_bin - 1);
    max_bin = MIN(hz->nbins - 1, max_bin + 1);
  }
#define MIN_BINS 12  // need enough data to estimate parameters of gaussian
  if (max_bin - min_bin < MIN_BINS) {
    max_bin += (MIN_BINS - (max_bin - min_bin)) / 2;
    min_bin -= (MIN_BINS - (max_bin - min_bin)) / 2;
    min_bin = MAX(0, min_bin);
    max_bin = MIN(hz->nbins - 1, max_bin);
  }
  zbin = HISTOfindBin(hz, 0);
  delta_sigma = h->bin_size / 10;
  max_sigma = habs->nbins * habs->bin_size / 2;
  best_sse = -1;
  best_sigma = 0;
  for (sigma = delta_sigma; sigma <= max_sigma; sigma += delta_sigma) {
    scale = hz->counts[zbin] * sigma * sqrt_2pi;
    predicted_val = scale / (sigma * sqrt_2pi) * exp(0 / (2 * sigma * sigma));

    scale = 0;
    for (n = 0, bin = min_bin; bin <= max_bin; bin++) {
      val = hz->bins[bin];
      predicted_val = 1 / (sigma * sqrt_2pi) * exp(-val * val / (2 * sigma * sigma));
      if (predicted_val > 1e-3) {
        scale += hz->counts[bin] / predicted_val;
        n++;
      }
    }
    scale /= n;
    sse = 0;
    for (bin = min_bin; bin <= max_bin; bin++) {
      val = hz->bins[bin];
      predicted_val = scale / (sigma * sqrt_2pi) * exp(-val * val / (2 * sigma * sigma));
      error = predicted_val - hz->counts[bin];

      sse += error * error;
    }

    if (sse < best_sse || best_sse < 0) {
      best_sse = sse;
      best_sigma = sigma;
    }
  }
  *psigma = best_sigma;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOGRAM *hp;

    hp = HISTOcopy(hz, NULL);

    sigma = best_sigma;
    scale = hz->counts[zbin] * best_sigma * sqrt_2pi;
    scale = 0;
    for (n = 0, bin = min_bin; bin <= max_bin; bin++) {
      val = hz->bins[bin] - mean;
      predicted_val = 1 / (sigma * sqrt_2pi) * exp(-val * val / (2 * sigma * sigma));
      if (predicted_val > 1e-3) {
        scale += hz->counts[bin] / predicted_val;
        n++;
      }
    }
    scale /= n;
    for (bin = 0; bin < hp->nbins; bin++) {
      val = hp->bins[bin] - mean;
      predicted_val = scale / (best_sigma * sqrt_2pi) * exp(-val * val / (2 * best_sigma * best_sigma));
      hp->counts[bin] = predicted_val;
    }
    HISTOplot(hp, "hp.plt");
    HISTOfree(&hp);
  }

  HISTOfree(&hs);
  HISTOfree(&habs);
  HISTOfree(&hcdf);
  HISTOfree(&hz);
  return (NO_ERROR);
}

HISTOGRAM *HISTOgaussianCDF(HISTOGRAM *h, double mean, double sigma, int nbins)
{
  HISTOGRAM *hpdf;

  hpdf = HISTOgaussianPDF(h, mean, sigma, nbins);
  h = HISTOmakeCDF(hpdf, h);
  HISTOfree(&hpdf);
  return (h);
}

HISTOGRAM *HISTOgaussianPDF(HISTOGRAM *h, double mean, double sigma, int nbins)
{
  int bin;
  double val, sqrt_2pi;

  sqrt_2pi = sqrt(2 * M_PI);
  if (h == NULL) h = HISTOalloc(nbins);

#define NSIGMAS 10
  h->bin_size = (2 * NSIGMAS) * sigma / nbins;
  h->min = mean - NSIGMAS * sigma;

  for (bin = 0; bin < h->nbins; bin++) {
    h->bins[bin] = h->min + h->bin_size * bin;
    val = h->bins[bin] - mean;
    val = 1 / (sigma * sqrt_2pi) * exp(-val * val / (2 * sigma * sigma));
    h->counts[bin] = val;
  }

  HISTOmakePDF(h, h);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) HISTOplot(h, "h.plt");
  return (h);
}
HISTOGRAM *HISTOabs(HISTOGRAM *h, HISTOGRAM *habs)
{
  int b, pbins, pbin;
  double mx, mn;

  mx = 0;
  mn = 1e10;
  for (b = pbins = 0; b < h->nbins; b++) {
    if (fabs(h->bins[b]) > mx) mx = fabs(h->bins[b]);
    if (fabs(h->bins[b]) < mn) mn = fabs(h->bins[b]);
  }

  if (mn >= mx) mn = mx * .9;
  habs = HISTOrealloc(habs, nint(ceil((mx - mn) / h->bin_size) + 1));
  habs->bin_size = h->bin_size;
  habs->min = mn;
  habs->max = mx;
  for (b = 0; b < habs->nbins; b++) habs->bins[b] = b * habs->bin_size + mn;

  for (b = 0; b < h->nbins; b++) {
    pbin = HISTOfindBin(habs, fabs(h->bins[b]));
    habs->counts[pbin] += h->counts[b];
  }
  return (habs);
}
HISTOGRAM *HISTOscalarMul(HISTOGRAM *hsrc, float mul, HISTOGRAM *hdst)
{
  int b;

  if (!hdst)
    hdst = HISTOalloc(hsrc->nbins);
  else
    hdst->nbins = hsrc->nbins;

  for (b = 0; b < hsrc->nbins; b++) hdst->counts[b] = hsrc->counts[b] * mul;

  return (hdst);
}

HISTOGRAM *HISTOscalarAdd(HISTOGRAM *hsrc, float add, HISTOGRAM *hdst)
{
  int b;

  if (!hdst)
    hdst = HISTOalloc(hsrc->nbins);
  else
    hdst->nbins = hsrc->nbins;

  for (b = 0; b < hsrc->nbins; b++) hdst->counts[b] = hsrc->counts[b] + add;

  return (hdst);
}

double HISTOgetEntropy(HISTOGRAM *h)
{
  double total = (double)HISTOtotal(h);
  double entropy = 0.0, p;
  int b;

  for (b = 0; b < h->nbins; b++) {
    if (h->counts[b] > 0) {
      p = (double)h->counts[b] / total;
      entropy -= p * log(p);
    }
  }
  return (entropy);
}

HISTOGRAM *HISTOsoapBubbleZeros(HISTOGRAM *hsrc, HISTOGRAM *hdst, int niters)
{
  int n, b, *control;
  double *tmp;

  tmp = (double *)calloc(hsrc->nbins, sizeof(double));
  control = (int *)calloc(hsrc->nbins, sizeof(int));

  hdst = HISTOcopy(hsrc, hdst);

  for (b = 0; b < hsrc->nbins; b++)
    if (hsrc->counts[b] > 0) control[b] = 1;  // in case hsrc == hdst

  for (n = 0; n < niters; n++) {
    for (b = 0; b < hsrc->nbins; b++) {
      if (control[b])
        tmp[b] = hsrc->counts[b];
      else {
        if (b == 0)
          tmp[b] = (hdst->counts[b] + hdst->counts[b + 1]) / 2;
        else if (b == hdst->nbins - 1)
          tmp[b] = (hdst->counts[b] + hdst->counts[b - 1]) / 2;
        else
          tmp[b] = (hdst->counts[b + 1] + hdst->counts[b] + hdst->counts[b - 1]) / 3;
      }
    }
    for (b = 0; b < hsrc->nbins; b++) hdst->counts[b] = tmp[b];
  }

  free(tmp);
  free(control);
  return (hdst);
}
int HISTOfindMaxDerivative(HISTOGRAM *h, double min_count, double max_count, int whalf, int grad_dir)
{
  int index, i0, peak_index, num = 0;
  double prev_val, next_val, max_d, d;

  max_d = 0;
  peak_index = -1;
  for (index = 0; index < h->nbins; index++) {
    for (prev_val = 0.0, num = 0, i0 = MAX(0, index - whalf); i0 < index; i0++, num++) prev_val += h->counts[i0];
    if (num == 0) continue;
    prev_val /= num;
    if (prev_val < min_count || prev_val > max_count) continue;
    for (next_val = 0.0, num = 0, i0 = MIN(h->nbins - 1, index + 1); i0 < MIN(index + whalf, h->nbins - 1); i0++, num++)
      next_val += h->counts[i0];
    if (num == 0) continue;
    next_val /= num;
    if (next_val < min_count || next_val > max_count) continue;
    d = (next_val - prev_val);
    if (d > max_d) {
      max_d = d;
      peak_index = index;
    }
  }

  return (peak_index);
}

/**************************   2D Histogram stuff starts *************** */

HISTOGRAM2D *HISTO2Dalloc(int nbins1, int nbins2)
{
  HISTOGRAM2D *histo;
  int i;

  histo = (HISTOGRAM2D *)calloc(1, sizeof(HISTOGRAM2D));
  if (!histo) ErrorExit(ERROR_NO_MEMORY, "HISTO2Dalloc(%d, %d): allocation failed", nbins1, nbins2);

  histo->bins1 = (float *)calloc(nbins1, sizeof(float));
  histo->bins2 = (float *)calloc(nbins2, sizeof(float));
  histo->counts = (float **)calloc(nbins1, sizeof(float *));
  if (!histo->counts || !histo->bins1 || !histo->bins2)
    ErrorExit(ERROR_NOMEMORY, "HISTO2Dalloc(%d, %d): could not allocate histogram", nbins1, nbins2);
  histo->nbins1 = nbins1;
  histo->nbins2 = nbins2;
  for (i = 0; i < nbins1; i++) {
    histo->counts[i] = (float *)calloc(nbins2, sizeof(float));
    if (!histo->counts[i])
      ErrorExit(ERROR_NOMEMORY, "HISTO2Dalloc(%d, %d): could not allocate histogram slice %d", nbins1, nbins2, i);
  }

  return (histo);
}
int HISTO2Dfree(HISTOGRAM2D **phisto)
{
  HISTOGRAM2D *histo;
  int i;

  histo = *phisto;
  *phisto = NULL;
  if (histo) {
    if (histo->bins1) free(histo->bins1);
    if (histo->bins2) free(histo->bins2);
    if (histo->counts) {
      for (i = 0; i < histo->nbins1; i++) free(histo->counts[i]);
      free(histo->counts);
    }
    else
      DiagBreak();
    free(histo);
  }
  else
    DiagBreak();

  return (NO_ERROR);
}

int HISTO2Ddump(HISTOGRAM2D *histo, FILE *fp)
{
  int bin1, bin2;

  if (!histo)
    fprintf(stderr, "NULL histogram");
  else {
    fprintf(fp, "nbins = %d, %d\n", histo->nbins1, histo->nbins2);
    for (bin1 = 0; bin1 < histo->nbins1; bin1++) {
      for (bin2 = 0; bin2 < histo->nbins2; bin2++)
        if (histo->counts[bin1][bin2])
          fprintf(fp,
                  "bin[%d, %d] = %2.1f, %2.1f = %2.2f\n",
                  bin1,
                  bin2,
                  histo->bins1[bin1],
                  histo->bins2[bin2],
                  histo->counts[bin1][bin2]);
    }
  }
  return (NO_ERROR);
}
int HISTO2Dwrite(HISTOGRAM2D *h, char *fname)
{
  FILE *fp;
  int ret;

  fp = fopen(fname, "wb");
  ret = HISTO2DwriteInto(h, fp);
  fclose(fp);
  return (ret);
}
int HISTO2DwriteInto(HISTOGRAM2D *h, FILE *fp)
{
  int b1, b2;

  fwriteInt(h->nbins1, fp);
  fwriteInt(h->nbins2, fp);
  fwriteFloat(h->bin_size1, fp);
  fwriteFloat(h->bin_size2, fp);
  fwriteFloat(h->min1, fp);
  fwriteFloat(h->min2, fp);
  fwriteFloat(h->max1, fp);
  fwriteFloat(h->max2, fp);
  for (b1 = 0; b1 < h->nbins1; b1++) fwriteFloat(h->bins1[b1], fp);
  for (b2 = 0; b2 < h->nbins2; b2++) fwriteFloat(h->bins2[b2], fp);

  for (b1 = 0; b1 < h->nbins1; b1++)
    for (b2 = 0; b2 < h->nbins2; b2++) fwriteFloat(h->counts[b1][b2], fp);
  return (NO_ERROR);
}
HISTOGRAM2D *HISTO2Dread(char *fname)
{
  HISTOGRAM2D *histo;
  FILE *fp;

  fp = fopen(fname, "rb");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "HISTO2Dread(%s); could not open file", fname));

  histo = HISTO2DreadFrom(fp);
  fclose(fp);
  return (histo);
}

HISTOGRAM2D *HISTO2DreadFrom(FILE *fp)
{
  int b1, nbins1, b2, nbins2;
  HISTOGRAM2D *h;

  nbins1 = freadInt(fp);
  nbins2 = freadInt(fp);
  h = HISTO2Dalloc(nbins1, nbins2);
  h->bin_size1 = freadFloat(fp);
  h->bin_size2 = freadFloat(fp);
  h->min1 = freadFloat(fp);
  h->min2 = freadFloat(fp);
  h->max1 = freadFloat(fp);
  h->max2 = freadFloat(fp);
  for (b1 = 0; b1 < h->nbins1; b1++) h->bins1[b1] = freadFloat(fp);
  for (b2 = 0; b2 < h->nbins2; b2++) h->bins2[b2] = freadFloat(fp);

  for (b1 = 0; b1 < h->nbins1; b1++)
    for (b2 = 0; b2 < h->nbins2; b2++) h->counts[b1][b2] = freadFloat(fp);
  return (h);
}
HISTOGRAM2D *HISTO2Dclear(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst)
{
  int b;

  if (!histo_dst) histo_dst = HISTO2Dalloc(histo_src->nbins1, histo_src->nbins2);

  memset(histo_dst->bins1, 0, histo_dst->nbins1 * sizeof(*histo_dst->bins1));
  memset(histo_dst->bins2, 0, histo_dst->nbins2 * sizeof(*histo_dst->bins2));
  for (b = 0; b < histo_dst->nbins1; b++)
    memset(histo_dst->counts[b], 0, histo_dst->nbins2 * sizeof(*histo_dst->counts[b]));

  return (histo_dst);
}
HISTOGRAM2D *HISTO2Dinit(HISTOGRAM2D *h, int nbins1, int nbins2, double mn1, int mx1, double mn2, double mx2)
{
  int b;

  if (h == NULL)
    h = HISTO2Dalloc(nbins1, nbins2);
  else
    HISTO2Dclear(h, h);

  h->min1 = mn1;
  h->max1 = mx1;
  h->min2 = mn2;
  h->max2 = mx2;
  h->bin_size1 = (mx1 - mn1) / (nbins1 - 1);
  h->bin_size2 = (mx2 - mn2) / (nbins2 - 1);
  HISTO2Dclear(h, h);
  for (b = 0; b < h->nbins1; b++) h->bins1[b] = mn1 + h->bin_size1 * (float)b;
  for (b = 0; b < h->nbins2; b++) h->bins2[b] = mn2 + h->bin_size2 * (float)b;

  return (h);
}
HISTOGRAM2D *HISTO2Drealloc(HISTOGRAM2D *histo, int nbins1, int nbins2)
{
  int b;
  if (histo == NULL) return (HISTO2Dalloc(nbins1, nbins2));

  if (histo->bins1) free(histo->bins1);
  if (histo->bins2) free(histo->bins2);
  if (histo->counts) {
    for (b = 0; b < histo->nbins1; b++) free(histo->counts[b]);
    free(histo->counts);
  }
  histo->bins1 = (float *)calloc(nbins1, sizeof(float));
  histo->bins2 = (float *)calloc(nbins2, sizeof(float));
  histo->counts = (float **)calloc(nbins1, sizeof(float *));
  if (!histo->counts || !histo->bins1 || !histo->bins2)
    ErrorExit(ERROR_NOMEMORY, "HISTO2Drealloc(%d, %d): could not allocate histogram", nbins1, nbins2);
  histo->nbins1 = nbins1;
  histo->nbins2 = nbins2;
  for (b = 0; b < histo->nbins1; b++) histo->counts[b] = (float *)calloc(nbins2, sizeof(float));

  return (histo);
}
int HISTO2DaddFractionalSample(
    HISTOGRAM2D *histo, float val1, float val2, float bmin1, float bmax1, float bmin2, float bmax2, float frac)
{
  int bin_no1, bin_no2;
  float bin_size1, bin_size2;

  if (bmin1 >= bmax1)  // not specified by caller
  {
    bmin1 = histo->min1;
    bmax1 = histo->max1;
    bin_size1 = histo->bin_size1;
  }
  else
    bin_size1 = (bmax1 - bmin1) / ((float)histo->nbins1 - 1);

  bin_no1 = nint((val1 - bmin1) / bin_size1);
  if (bmin2 >= bmax2)  // not specified by caller
  {
    bmin2 = histo->min2;
    bmax2 = histo->max2;
    bin_size2 = histo->bin_size2;
  }
  else
    bin_size2 = (bmax2 - bmin2) / ((float)histo->nbins2 - 1);

  bin_no2 = nint((val2 - bmin2) / bin_size2);

  if ((bin_no1 < histo->nbins1 && bin_no1 >= 0) && (bin_no2 < histo->nbins2 && bin_no2 >= 0))
    histo->counts[bin_no1][bin_no2] += frac;
  else
    DiagBreak();

  return (NO_ERROR);
}
int HISTO2DaddSample(HISTOGRAM2D *histo, float val1, float val2, float bmin1, float bmax1, float bmin2, float bmax2)
{
  int bin_no1, bin_no2;
  float bin_size1, bin_size2;

  if (bmin1 >= bmax1)  // not specified by caller
  {
    bmin1 = histo->min1;
    bmax1 = histo->max1;
    bin_size1 = histo->bin_size1;
  }
  else
    bin_size1 = (bmax1 - bmin1) / ((float)histo->nbins1 - 1);
  bin_no1 = nint((val1 - bmin1) / bin_size1);

  if (bmin2 >= bmax2)  // not specified by caller
  {
    bmin2 = histo->min2;
    bmax2 = histo->max2;
    bin_size2 = histo->bin_size2;
  }
  else
    bin_size2 = (bmax2 - bmin2) / ((float)histo->nbins2 - 1);
  bin_no2 = nint((val2 - bmin2) / bin_size2);

  if ((bin_no1 < histo->nbins1 && bin_no1 >= 0) && (bin_no2 < histo->nbins2 && bin_no2 >= 0))
    histo->counts[bin_no1][bin_no2]++;
  else
    DiagBreak();

  return (NO_ERROR);
}
int HISTO2Dplot(HISTOGRAM2D *histo, char *fname)
{
  FILE *fp;
  int bin_no1, bmin1, bmax1, bin_no2, bmin2, bmax2;

  fp = fopen(fname, "w");
  if (fp == NULL) return (ERROR_NOFILE);

  bmin2 = 0;
  bmax2 = histo->nbins2 - 1;
  for (bmin1 = 0; bmin1 < histo->nbins1; bmin1++) {
    for (bmin2 = 0; bmin2 < histo->nbins2; bmin2++)
      if (histo->counts[bmin1][bmin2] > 0) break;
    if (bmin2 < histo->nbins2 && histo->counts[bmin1][bmin2] > 0) break;
  }
  for (bmax1 = histo->nbins1 - 1; bmax1 > bmin1; bmax1--) {
    for (bmax2 = histo->nbins2 - 1; bmax2 > bmin2; bmax2--)
      if (histo->counts[bmax1][bmax2] > 0) break;
    if (bmax2 > bmin2 && histo->counts[bmax1][bmax2] > 0) break;
  }

  for (bin_no1 = bmin1; bin_no1 <= bmax1; bin_no1++) fprintf(fp, "%f ", histo->bins1[bin_no1]);
  fprintf(fp, "\n");
  for (bin_no2 = bmin2; bin_no2 <= bmax2; bin_no2++) fprintf(fp, "%f ", histo->bins2[bin_no2]);
  fprintf(fp, "\n");
  for (bin_no1 = bmin1; bin_no1 <= bmax1; bin_no1++) {
    for (bin_no2 = bmin2; bin_no2 <= bmax2; bin_no2++) fprintf(fp, "%f ", histo->counts[bin_no1][bin_no2]);
    fprintf(fp, "\n");
  }

  fclose(fp);

  return (NO_ERROR);
}
HISTOGRAM2D *HISTO2Dcopy(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst)
{
  int b;

  if (!histo_dst) histo_dst = HISTO2Dalloc(histo_src->nbins1, histo_src->nbins2);
  histo_dst->nbins1 = histo_src->nbins1;
  histo_dst->nbins2 = histo_src->nbins2;
  histo_dst->bin_size1 = histo_src->bin_size1;
  histo_dst->bin_size2 = histo_src->bin_size2;
  // use memmove, not memcpy, as src and dst could overlap
  for (b = 0; b < histo_dst->nbins1; b++)
    memmove(histo_dst->counts[b], histo_src->counts[b], sizeof(*histo_src->counts[b]) * histo_src->nbins2);

  memmove(histo_dst->bins1, histo_src->bins1, sizeof(*histo_src->bins1) * histo_src->nbins1);
  memmove(histo_dst->bins2, histo_src->bins2, sizeof(*histo_src->bins2) * histo_src->nbins2);
  histo_dst->min1 = histo_src->min1;
  histo_dst->max1 = histo_src->max1;
  histo_dst->min2 = histo_src->min2;
  histo_dst->max2 = histo_src->max2;

  return (histo_dst);
}
HISTOGRAM2D *HISTO2DmakePDF(HISTO2D *h_src, HISTO2D *h_dst)
{
  int b1, b2;
  double total, unlikely;

  h_dst = HISTO2Dcopy(h_src, h_dst);

  for (total = 0.0, b1 = 0; b1 < h_dst->nbins1; b1++)
    for (b2 = 0; b2 < h_dst->nbins2; b2++) total += h_dst->counts[b1][b2];

  if (total > 0) {
    unlikely = 1.0 / (10 * total);
    for (b1 = 0; b1 < h_dst->nbins1; b1++)
      for (b2 = 0; b2 < h_dst->nbins2; b2++) {
        h_dst->counts[b1][b2] /= total;
        if (DZERO(h_dst->counts[b1][b2])) h_dst->counts[b1][b2] = unlikely;
      }
  }

  return (h_dst);
}
int HISTO2DfindBin1(HISTOGRAM2D *h, float val)
{
  int b;

  if (h->bin_size1 == 1) return ((int)val - h->min1);

  for (b = h->nbins1 - 1; b > 0; b--)
    if (h->bins1[b - 1] < val) return (b);

  return (0);
}
int HISTO2DfindBin2(HISTOGRAM2D *h, float val)
{
  int b;

  if (h->bin_size2 == 1) return ((int)val - h->min2);

  for (b = h->nbins2 - 1; b > 0; b--)
    if (h->bins2[b - 1] < val) return (b);

  return (0);
}
double HISTO2DgetCount(HISTOGRAM2D *h, float bin_val1, float bin_val2)
{
  int b1, b2;

  b1 = HISTO2DfindBin1(h, bin_val1);
  if (b1 < 0 || b1 >= h->nbins1) return (0);
  b2 = HISTO2DfindBin2(h, bin_val2);
  if (b2 < 0 || b2 >= h->nbins2) return (0);
  return (h->counts[b1][b2]);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MAX_LEN 2000
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MAX_LEN 2000
HISTOGRAM2D *HISTO2Dsmooth(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst, float sigma)
{
  float norm, two_sigma, fx, k, kernel[MAX_LEN], total;
  int x, half, len, b1, b2, kx, nbins1, nbins2, y, ky, b1k, b2k, i;

  nbins1 = histo_src->nbins1;
  nbins2 = histo_src->nbins2;
  if (!histo_dst) {
    histo_dst = HISTO2Dalloc(nbins1, nbins2);
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->min1 = histo_src->min1;
    histo_dst->max1 = histo_src->max1;
    histo_dst->bin_size2 = histo_src->bin_size2;
    histo_dst->min2 = histo_src->min2;
    histo_dst->max2 = histo_src->max2;
  }
  else {
    if ((histo_dst->nbins1 < histo_src->nbins1) || (histo_dst->nbins2 < histo_src->nbins2)) {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTO2Drealloc(histo_dst, nbins1, nbins2);
    }
    histo_dst->nbins1 = nbins1;
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->nbins2 = nbins2;
    histo_dst->bin_size2 = histo_src->bin_size2;
    HISTO2Dclear(histo_dst, histo_dst);
  }
  /* build the kernel in k */
  len = (int)nint(sqrt(2) * 8.0f * sigma) + 1;
  if (ISEVEN(len)) /* ensure it's odd */
    len++;
  if (MAX_LEN && (MAX_LEN < len)) len = MAX_LEN;
  half = len / 2;

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

    kernel[x] = k;
    norm += k;
  }
  for (x = 0; x < len; x++) kernel[x] /= norm;

  for (b1 = 0; b1 < nbins1; b1++) {
    histo_dst->bins1[b1] = histo_src->bins1[b1];
    for (b2 = 0; b2 < nbins2; b2++) {
      histo_dst->bins2[b2] = histo_src->bins2[b2];
      for (norm = 0.0, total = 0.0f, x = 0; x < len; x++)
        for (y = 0; y < len; y++) {
          kx = x - half;
          ky = y - half;
          b1k = b1 + kx;
          b2k = b2 + ky;
          if ((b1k >= nbins1 || b1k < 0) || (b2k >= nbins2 || b2k < 0)) continue;

          i = nint(sqrt(x * x + y * y));
          if (i > len - 1) i = len - 1;
          norm += kernel[i];
          total += kernel[i] * (float)histo_src->counts[b1k][b2k];
          if ((float)histo_src->counts[b1k][b2k] > 0.00027) DiagBreak();
          if (!std::isfinite(norm) || !std::isfinite(total) || !std::isfinite(kernel[i])) DiagBreak();
        }
      if (b1 == Gx && b2 == Gy) DiagBreak();
      if (DZERO(norm)) {
        DiagBreak();
        histo_dst->counts[b1][b2] = histo_src->counts[b1][b2];
        continue;
      }

      histo_dst->counts[b1][b2] = total / norm;
    }
  }

  return (histo_dst);
}
#if 1
HISTOGRAM2D *HISTO2DsmoothAnisotropic(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst, float sigma1, float sigma2)
{
  float norm, k, total;
  int x, half1, half2, len1, len2, b1, b2, kx, nbins1, nbins2, y, ky, b1k, b2k, dist;

  nbins1 = histo_src->nbins1;
  nbins2 = histo_src->nbins2;
  if (!histo_dst) {
    histo_dst = HISTO2Dalloc(nbins1, nbins2);
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->min1 = histo_src->min1;
    histo_dst->max1 = histo_src->max1;
    histo_dst->bin_size2 = histo_src->bin_size2;
    histo_dst->min2 = histo_src->min2;
    histo_dst->max2 = histo_src->max2;
  }
  else {
    if ((histo_dst->nbins1 < histo_src->nbins1) || (histo_dst->nbins2 < histo_src->nbins2)) {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTO2Drealloc(histo_dst, nbins1, nbins2);
    }
    histo_dst->nbins1 = nbins1;
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->nbins2 = nbins2;
    histo_dst->bin_size2 = histo_src->bin_size2;
    HISTO2Dclear(histo_dst, histo_dst);
  }

  /* build the kernel in k for bins1 */
  len1 = (int)nint(8.0f * sigma1) + 1;
  if (ISEVEN(len1)) /* ensure it's odd */
    len1++;
  half1 = len1 / 2;
  len2 = (int)nint(8.0f * sigma2) + 2;
  if (MAX_LEN && (MAX_LEN < len1)) len1 = MAX_LEN;
  if (MAX_LEN && (MAX_LEN < len2)) len2 = MAX_LEN;
  if (ISEVEN(len2)) /* ensure it's odd */
    len2++;
  half2 = len2 / 2;

  for (b1 = 0; b1 < nbins1; b1++) {
    histo_dst->bins1[b1] = histo_src->bins1[b1];
    for (b2 = 0; b2 < nbins2; b2++) {
      histo_dst->bins2[b2] = histo_src->bins2[b2];
      for (norm = 0.0, total = 0.0f, x = 0; x < len1; x++)
        for (y = 0; y < len2; y++) {
          kx = x - half1;
          ky = y - half2;
          b1k = b1 + kx;
          b2k = b2 + ky;
          if ((b1k >= nbins1 || b1k < 0) || (b2k >= nbins2 || b2k < 0)) continue;

          dist = (SQR(kx / sigma1) + SQR(ky / sigma2));
          k = exp(-0.5 * dist);
          norm += k;
          total += k * (float)histo_src->counts[b1k][b2k];
          if ((float)histo_src->counts[b1k][b2k] > 0.00027) DiagBreak();
          if (!std::isfinite(norm) || !std::isfinite(total) || !std::isfinite(k)) DiagBreak();
        }
      if (b1 == Gx && b2 == Gy) DiagBreak();
      if (DZERO(norm)) {
        DiagBreak();
        histo_dst->counts[b1][b2] = histo_src->counts[b1][b2];
        continue;
      }

      histo_dst->counts[b1][b2] = total / norm;
    }
  }

  return (histo_dst);
}
#else
#if 1
HISTOGRAM2D *HISTO2DsmoothAnisotropic(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst, float sigma1, float sigma2)
{
  float norm, two_sigma, fx, k, k1, k2, kernel1[MAX_LEN], kernel2[MAX_LEN], total;
  int x, half1, half2, len1, len2, b1, b2, kx, nbins1, nbins2, y, ky, b1k, b2k;

  nbins1 = histo_src->nbins1;
  nbins2 = histo_src->nbins2;
  if (!histo_dst) {
    histo_dst = HISTO2Dalloc(nbins1, nbins2);
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->min1 = histo_src->min1;
    histo_dst->max1 = histo_src->max1;
    histo_dst->bin_size2 = histo_src->bin_size2;
    histo_dst->min2 = histo_src->min2;
    histo_dst->max2 = histo_src->max2;
  }
  else {
    if ((histo_dst->nbins1 < histo_src->nbins1) || (histo_dst->nbins2 < histo_src->nbins2)) {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTO2Drealloc(histo_dst, nbins1, nbins2);
    }
    histo_dst->nbins1 = nbins1;
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->nbins2 = nbins2;
    histo_dst->bin_size2 = histo_src->bin_size2;
    HISTO2Dclear(histo_dst, histo_dst);
  }

  /* build the kernel in k for bins1 */
  len1 = (int)nint(sqrt(2) * 8.0f * sigma1) + 1;
  if (ISEVEN(len1)) /* ensure it's odd */
    len1++;
  if (MAX_LEN && (MAX_LEN < len1)) len1 = MAX_LEN;
  half1 = len1 / 2;

  norm = 0.0f;
  two_sigma = 2.0f * sigma1;

  for (norm = 0.0f, x = 0; x < len1; x++) {
    fx = (float)(x - half1);
    if (fabs(fx) <= two_sigma)
      k = (float)exp((double)(-fx * fx / (two_sigma * sigma1)));
    else if (two_sigma < fabs(fx) && fabs(fx) <= 4.0f * sigma1)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0f - fabs(fx) / (double)sigma1, 4.0);
    else
      k = 0;

    kernel1[x] = k;
    norm += k;
  }
  for (x = 0; x < len1; x++) kernel1[x] /= norm;

  /* build the kernel in k for bins2 */
  len2 = (int)nint(sqrt(2) * 8.0f * sigma2) + 1;
  if (ISEVEN(len2)) /* ensure it's odd */
    len2++;
  if (MAX_LEN && (MAX_LEN < len2)) len2 = MAX_LEN;
  half2 = len2 / 2;

  norm = 0.0f;
  two_sigma = 2.0f * sigma2;

  for (norm = 0.0f, x = 0; x < len2; x++) {
    fx = (float)(x - half2);
    if (fabs(fx) <= two_sigma)
      k = (float)exp((double)(-fx * fx / (two_sigma * sigma2)));
    else if (two_sigma < fabs(fx) && fabs(fx) <= 4.0f * sigma2)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0f - fabs(fx) / (double)sigma2, 4.0);
    else
      k = 0;

    kernel2[x] = k;
    norm += k;
  }
  for (x = 0; x < len1; x++) kernel2[x] /= norm;

  for (b1 = 0; b1 < nbins1; b1++) {
    histo_dst->bins1[b1] = histo_src->bins1[b1];
    for (b2 = 0; b2 < nbins2; b2++) {
      histo_dst->bins2[b2] = histo_src->bins2[b2];
      for (norm = 0.0, total = 0.0f, x = 0; x < len1; x++)
        for (y = 0; y < len2; y++) {
          kx = x - half1;
          ky = y - half2;
          b1k = b1 + kx;
          b2k = b2 + ky;
          if ((b1k >= nbins1 || b1k < 0) || (b2k >= nbins2 || b2k < 0)) continue;

          k1 = kernel1[x];
          k2 = kernel2[y];
          k = k1 * k2;
          norm += k;
          total += k * (float)histo_src->counts[b1k][b2k];
          if ((float)histo_src->counts[b1k][b2k] > 0.00027) DiagBreak();
          if (!std::isfinite(norm) || !std::isfinite(total) || !std::isfinite(k)) DiagBreak();
        }
      if (b1 == Gx && b2 == Gy) DiagBreak();
      if (DZERO(norm)) {
        DiagBreak();
        histo_dst->counts[b1][b2] = histo_src->counts[b1][b2];
        continue;
      }

      histo_dst->counts[b1][b2] = total / norm;
    }
  }

  return (histo_dst);
}
#else
HISTOGRAM2D *HISTO2DsmoothAnisotropic(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst, float sigma1, float sigma2)
{
  HISTO2D *h_tmp;
  h_tmp = HISTO2DsmoothBins1(histo_src, NULL, sigma1);
  histo_dst = HISTO2DsmoothBins2(h_tmp, histo_dst, sigma2);
  HISTO2Dfree(&h_tmp);
  return (histo_dst);
}
#endif
#endif
HISTOGRAM2D *HISTO2DsmoothBins1(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst, float sigma)
{
  float norm, two_sigma, fx, k, kernel[MAX_LEN], total;
  int x, half, len, b1, b2, kx, nbins1, nbins2, b1k;

  nbins1 = histo_src->nbins1;
  nbins2 = histo_src->nbins2;
  if (!histo_dst) {
    histo_dst = HISTO2Dalloc(nbins1, nbins2);
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->min1 = histo_src->min1;
    histo_dst->max1 = histo_src->max1;
    histo_dst->bin_size2 = histo_src->bin_size2;
    histo_dst->min2 = histo_src->min2;
    histo_dst->max2 = histo_src->max2;
  }
  else {
    if ((histo_dst->nbins1 < histo_src->nbins1) || (histo_dst->nbins2 < histo_src->nbins2)) {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTO2Drealloc(histo_dst, nbins1, nbins2);
    }
    histo_dst->nbins1 = nbins1;
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->nbins2 = nbins2;
    histo_dst->bin_size2 = histo_src->bin_size2;
    HISTO2Dclear(histo_dst, histo_dst);
  }
  /* build the kernel in k */
  len = (int)nint(8.0f * sigma) + 1;
  if (ISEVEN(len)) /* ensure it's odd */
    len++;
  if (MAX_LEN && (MAX_LEN < len)) len = MAX_LEN;
  half = len / 2;

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

    kernel[x] = k;
    norm += k;
  }
  for (x = 0; x < len; x++) kernel[x] /= norm;

  for (b1 = 0; b1 < nbins1; b1++) {
    histo_dst->bins1[b1] = histo_src->bins1[b1];
    for (b2 = 0; b2 < nbins2; b2++) {
      histo_dst->bins2[b2] = histo_src->bins2[b2];
      for (norm = 0.0, total = 0.0f, x = 0; x < len; x++) {
        kx = x - half;
        b1k = b1 + kx;
        if (b1k >= nbins1 || b1k < 0) continue;

        norm += kernel[x];
        total += kernel[x] * (float)histo_src->counts[b1k][b2];
        if (!std::isfinite(norm) || !std::isfinite(total) || !std::isfinite(kernel[x])) DiagBreak();
      }
      if (b1 == Gx && b2 == Gy) DiagBreak();
      if (DZERO(norm)) {
        DiagBreak();
        histo_dst->counts[b1][b2] = histo_src->counts[b1][b2];
        continue;
      }

      histo_dst->counts[b1][b2] = total / norm;
    }
  }

  return (histo_dst);
}

HISTOGRAM2D *HISTO2DsmoothBins2(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst, float sigma)
{
  float norm, two_sigma, fx, k, kernel[MAX_LEN], total;
  int x, kx, half, len, b1, b2, nbins1, nbins2, b2k;

  nbins1 = histo_src->nbins1;
  nbins2 = histo_src->nbins2;
  if (!histo_dst) {
    histo_dst = HISTO2Dalloc(nbins1, nbins2);
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->min1 = histo_src->min1;
    histo_dst->max1 = histo_src->max1;
    histo_dst->bin_size2 = histo_src->bin_size2;
    histo_dst->min2 = histo_src->min2;
    histo_dst->max2 = histo_src->max2;
  }
  else {
    if ((histo_dst->nbins1 < histo_src->nbins1) || (histo_dst->nbins2 < histo_src->nbins2)) {
      // fprintf(stderr, "realloc: histo_dst->nbins = %d, "
      //"histo_src->nbins = %d\n",
      //         histo_dst->nbins, histo_src->nbins);
      HISTO2Drealloc(histo_dst, nbins1, nbins2);
    }
    histo_dst->nbins1 = nbins1;
    histo_dst->bin_size1 = histo_src->bin_size1;
    histo_dst->nbins2 = nbins2;
    histo_dst->bin_size2 = histo_src->bin_size2;
    HISTO2Dclear(histo_dst, histo_dst);
  }
  /* build the kernel in k */
  len = (int)nint(8.0f * sigma) + 1;
  if (ISEVEN(len)) /* ensure it's odd */
    len++;
  if (MAX_LEN && (MAX_LEN < len)) len = MAX_LEN;
  half = len / 2;

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

    kernel[x] = k;
    norm += k;
  }
  for (x = 0; x < len; x++) kernel[x] /= norm;

  for (b1 = 0; b1 < nbins1; b1++) {
    histo_dst->bins1[b1] = histo_src->bins1[b1];
    for (b2 = 0; b2 < nbins2; b2++) {
      histo_dst->bins2[b2] = histo_src->bins2[b2];
      for (norm = 0.0, total = 0.0f, x = 0; x < len; x++) {
        kx = x - half;
        b2k = b2 + kx;
        if (b2k >= nbins2 || b2k < 0) continue;

        norm += kernel[x];
        total += kernel[x] * (float)histo_src->counts[b1][b2k];
        if (!std::isfinite(norm) || !std::isfinite(total) || !std::isfinite(kernel[x])) DiagBreak();
      }
      if (b1 == Gx && b2 == Gy) DiagBreak();
      if (DZERO(norm)) {
        DiagBreak();
        histo_dst->counts[b1][b2] = histo_src->counts[b1][b2];
        continue;
      }

      histo_dst->counts[b1][b2] = total / norm;
    }
  }

  return (histo_dst);
}

HISTOGRAM2D *HISTO2DsoapBubbleZeros(HISTOGRAM2D *hsrc, HISTOGRAM2D *hdst, int niters)
{
  int n, b1, b2, **control, num;
  double **tmp, val;

  tmp = (double **)calloc(hsrc->nbins1, sizeof(double *));
  control = (int **)calloc(hsrc->nbins1, sizeof(int *));
  for (b1 = 0; b1 < hsrc->nbins1; b1++) {
    tmp[b1] = (double *)calloc(hsrc->nbins2, sizeof(double));
    control[b1] = (int *)calloc(hsrc->nbins2, sizeof(int));
  }

  hdst = HISTO2Dcopy(hsrc, hdst);

  for (b1 = 0; b1 < hsrc->nbins1; b1++)
    for (b2 = 0; b2 < hsrc->nbins2; b2++)
      if (hsrc->counts[b1][b2] > 0) control[b1][b2] = 1;  // in case hsrc == hdst

  // first set values of all zero bins that are adjacent to non-zero to the avg of the nbrs
  for (b1 = 0; b1 < hsrc->nbins1; b1++)
    for (b2 = 0; b2 < hsrc->nbins2; b2++) {
      if (control[b1][b2]) continue;

      val = 0;
      num = 0;
      if (b1 > 0) {
        if (control[b1 - 1][b2]) {
          val += hdst->counts[b1 - 1][b2];
          num++;
        }
      }
      if (b2 > 0) {
        if (control[b1][b2 - 1]) {
          val += hdst->counts[b1][b2 - 1];
          num++;
        }
      }
      if (b2 < hdst->nbins2 - 1) {
        if (control[b1][b2 + 1]) {
          val += hdst->counts[b1][b2 + 1];
          num++;
        }
      }
      if (b1 < hdst->nbins1 - 1) {
        if (control[b1 + 1][b2]) {
          val += hdst->counts[b1 + 1][b2];
          num++;
        }
      }
      if (num > 0) hdst->counts[b1][b2] = val / num;
    }

  for (n = 0; n < niters; n++) {
    for (b1 = 0; b1 < hsrc->nbins1; b1++)
      for (b2 = 0; b2 < hsrc->nbins2; b2++) {
        if (control[b1][b2])
          tmp[b1][b2] = hsrc->counts[b1][b2];
        else {
          val = hsrc->counts[b1][b2];
          num = 1;
          if (b1 > 0) {
            val += hdst->counts[b1 - 1][b2];
            num++;
          }
          if (b2 > 0) {
            val += hdst->counts[b1][b2 - 1];
            num++;
          }
          if (b2 < hdst->nbins2 - 1) {
            val += hdst->counts[b1][b2 + 1];
            num++;
          }
          if (b1 < hdst->nbins1 - 1) {
            val += hdst->counts[b1 + 1][b2];
            num++;
          }
          tmp[b1][b2] = val / num;
        }
      }
    for (b1 = 0; b1 < hsrc->nbins1; b1++)
      for (b2 = 0; b2 < hsrc->nbins2; b2++) hdst->counts[b1][b2] = tmp[b1][b2];
  }

  for (b1 = 0; b1 < hsrc->nbins1; b1++) {
    free(tmp[b1]);
    free(control[b1]);
  }
  free(tmp);
  free(control);
  return (hdst);
}
float HISTOcomputeFWHM(HISTOGRAM *h, int peak)
{
  int width, max_width;
  float thresh, fwhm;

  thresh = h->counts[peak] / 2;
  max_width = MIN(peak, h->nbins - (peak + 1));

  for (width = 1; width < max_width; width++) {
    if (h->counts[peak - width] < thresh && h->counts[peak + width] < thresh) break;
  }

  fwhm = 2 * width * h->bin_size;
  return (fwhm);
}

double HISTOrmsDifference(HISTOGRAM *h1, HISTOGRAM *h2)
{
  int b;
  double rms, error;

  for (rms = 0.0, b = 0; b < h1->nbins; b++) {
    error = h1->counts[b] - h2->counts[b];
    rms += error * error;
  }
  return (sqrt(rms / h1->nbins));
}
double HISTOearthMoversDistance(HISTOGRAM *h1, HISTOGRAM *h2)
{
  double emd = 0.0;

  h1 = HISTOcopy(h1, NULL);
  HISTOmakePDF(h1, h1);
  h2 = HISTOcopy(h2, NULL);
  HISTOmakePDF(h2, h2);

  HISTOfree(&h1);
  HISTOfree(&h2);  // these were copied at the beginning, not the caller's version
  return (emd);
}

double HISTOksDistance(HISTOGRAM *h1, HISTOGRAM *h2)
{
  double ks_dist, dist, c2;
  int b;

  h1 = HISTOcopy(h1, NULL);
  HISTOmakeCDF(h1, h1);
  h2 = HISTOcopy(h2, NULL);
  HISTOmakeCDF(h2, h2);

  for (ks_dist = 0.0, b = 0; b < h1->nbins; b++) {
    c2 = HISTOgetCount(h2, h1->bins[b]);
    dist = fabs(h1->counts[b] - c2);
    if (dist > ks_dist) ks_dist = dist;
  }

  HISTOfree(&h1);
  HISTOfree(&h2);  // these were copied at the beginning, not the caller's version
  return (ks_dist);
}

HISTOGRAM *HISTOeraseRightmostPeak(HISTOGRAM *hsrc, HISTOGRAM *hdst, int whalf, float min_pct, int min_val, int max_val)
{
  int bmin, bmax, b;

  hdst = HISTOcopy(hsrc, hdst);

  bmin = HISTOfindBin(hdst, min_val);
  bmax = HISTOfindBin(hdst, max_val);
  if (bmax <= 0 || bmax >= hdst->nbins) bmax = hdst->nbins - 1;
  for (b = bmax; b >= bmin; b--) {
    if (HISTOisPeak(hdst, b, whalf)) break;
  }

  if (b >= bmin)  // find a peak
    HISTOerase(hdst, MAX(0, b - 1), bmax);
  return (hdst);
}
int HISTOerase(HISTOGRAM *h, int bmin, int bmax)
{
  int b;

  bmin = MAX(0, bmin);
  bmax = MIN(h->nbins - 1, bmax);
  for (b = bmin; b <= bmax; b++) h->counts[b] = 0;
  if (bmax == h->nbins - 1) h->nbins = bmin - 1;
  return (NO_ERROR);
}
int HISTOisPeak(HISTOGRAM *h, int bin, int whalf)
{
  int b, bmin, bmax;
  float val;

  bmin = MAX(0, bin - whalf);
  bmax = MIN(h->nbins - 1, bin + whalf);
  val = h->counts[bin];
  if (bmin > bmax) return (0);
  for (b = bmin; b <= bmax; b++) {
    if (b == bin) continue;
    if (h->counts[b] >= val) return (0);
  }
  return (1);  // h->counts[bin] > all other vals in 2*whalf+1
}

/*!
  \fn int HISTOwrite(HISTOGRAM *histo, char *fname) 
  \brief Writes histogram into the given filename. This uses HISTOwriteInto()
  which writes in binary format. 
 */
int HISTOwrite(HISTOGRAM *histo, char *fname) 
{
  FILE *fp;
  fp = fopen(fname,"w");
  if(fp==NULL){
    printf("ERROR: HISTOwrite(): opening %s\n",fname);
    return(1);
  }
  HISTOwriteInto(histo, fp);
  fclose(fp);
  return(0);
}

/*!
  \fn int HISTOwriteTxt(HISTOGRAM *histo, char *fname) 
  \brief Writes histogram into the given filename in ascii format.
 */
int HISTOwriteTxt(HISTOGRAM *histo, const char *fname) 
{
  int n;
  FILE *fp;
  fp = fopen(fname,"w");
  if(fp==NULL){
    printf("ERROR: HISTOwrite(): opening %s\n",fname);
    return(1);
  }
  for(n=0; n < histo->nbins; n++)
    fprintf(fp,"%4d %10.5f %10.5f\n",n,histo->bins[n],histo->counts[n]);
  fclose(fp);
  return(0);
}

/*!
  \fn int HISTOsumNorm(HISTOGRAM *h)
  \brief Normalize the histogram to have sum=1
 */
int HISTOsumNorm(HISTOGRAM *h)
{
  int n;
  double sum=0;
  for(n=0; n < h->nbins; n++) sum += h->counts[n];
  if(sum != 0){
    for(n=0; n < h->nbins; n++) h->counts[n] /= sum;
  }
  return(0);
}

/*!
  \fn HISTOGRAM *HISTOcumsum(HISTOGRAM *h, HISTOGRAM *hout)
  \brief Compute the cumulative sum of the histogram. If the histogram is 
  is a PDF, then this creates a CDF.
 */
HISTOGRAM *HISTOcumsum(HISTOGRAM *h, HISTOGRAM *hout)
{
  int n;
  double sum;
  if(hout==NULL) hout = HISTOcopy(h,hout);
  sum = 0;
  for(n=0; n < h->nbins; n++) {
    sum += h->counts[n];
    hout->counts[n] = sum;
  }
  return(hout);
}


