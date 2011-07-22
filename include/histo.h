/**
 * @file  histo.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/07/22 12:50:27 $
 *    $Revision: 1.43 $
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


#ifndef HISTO_H
#define HISTO_H

#include <stdio.h>

#ifndef UCHAR
#define UCHAR  unsigned char
#endif

#define MAX_BINS  10000
typedef struct
{
  int     nbins ;
  float   *bins ;   /* upper end of the range that maps to this bin */
  float   *counts ; /* # of voxels which map to this bin */
  float   bin_size ;
  float   min ;
  float   max ;     // min and max vals in the histo
}
HISTOGRAM, HISTO ;

int       HISTOfree(HISTOGRAM **phisto) ;
int       HISTOdump(HISTOGRAM *histo, FILE *fp) ;
int       HISTOwrite(HISTOGRAM *histo, char *fname) ;
HISTOGRAM *HISTOinit(HISTOGRAM *h, int nbins, double mn, double mx) ;
HISTOGRAM *HISTOalloc(int nbins) ;
HISTOGRAM *HISTOrealloc(HISTOGRAM *histo, int nbins) ;
HISTOGRAM *HISTOcrunch(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
double     HISTOcorrelate(HISTOGRAM *h1, HISTOGRAM *h2) ;
HISTOGRAM *HISTOcopy(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOinvert(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,int max_dst);
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,
                          int max_out) ;
HISTOGRAM *HISTOclear(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
int       HISTOclearZeroBin(HISTOGRAM *h) ;
HISTOGRAM *HISTOclearCounts(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOfillZeros(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOcompose(HISTOGRAM *histo1, HISTOGRAM *histo2,
                        HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,
                          int max_val) ;
HISTOGRAM *HISTOcomposeInvert(HISTOGRAM *histo_fwd, HISTOGRAM *histo_inv,
                              HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOmul(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOscalarMul(HISTOGRAM *hsrc, float mul, HISTOGRAM *hdst) ;
HISTOGRAM *HISTOscalarAdd(HISTOGRAM *hsrc, float add, HISTOGRAM *hdst) ;
HISTOGRAM *HISTOadd(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOsubtract(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOclearBins(HISTOGRAM *h1, HISTOGRAM *h2, int b0, int b1) ;
HISTOGRAM *HISTOsmooth(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,float sigma);
int       HISTOfindLastPeak(HISTOGRAM *h, int wsize, float min_pct) ;
int       HISTOfindLastPeakRelative(HISTOGRAM *h, int wsize, float min_pct) ;
int       HISTOfindFirstPeak(HISTOGRAM *h, int wsize, float min_pct) ;
int       HISTOfindFirstPeakRelative(HISTOGRAM *h, int wsize, float min_pct);
int       HISTOfindValley(HISTOGRAM *h, int wsize, int b0, int b1) ;
int       HISTOfindNextValley(HISTOGRAM *h, int b0) ;
int       HISTOfindNextPeak(HISTOGRAM *h, int b0, int whalf) ;
int       HISTOfindPreviousValley(HISTOGRAM *h, int b0) ;
int       HISTOfindEndOfPeak(HISTOGRAM *h, int b0, float pct_peak) ;
int       HISTOfindStartOfPeak(HISTOGRAM *h, int b0, float pct_peak) ;
int       HISTOfindLastPeakInRegion(HISTOGRAM *h, int wsize, float min_pct,
                                    int b0, int b1) ;
int       HISTOcountPeaksInRegion(HISTOGRAM *h, int wsize, float min_pct,
                                  int *peaks, int max_peaks, int b0, int b1) ;
int       HISTOfindFirstPeakInRegion(HISTOGRAM *h, int wsize, float min_pct,
                                     int b0, int b1) ;
int       HISTOfindHighestPeakInRegion(HISTOGRAM *h, int b0, int b1);
int       HISTOplot(HISTOGRAM *histo, char *fname) ;
int       HISTOaddSample(HISTOGRAM *histo, float val, float bmin, float bmax) ;
int       HISTOfindCurrentPeak(HISTOGRAM *histo,
                               int b0,
                               int wsize,
                               float min_pct) ;
int       HISTOfillHoles(HISTO *h) ;
float     HISTOtotal(HISTO *h) ;
float     HISTOtotalInRegion(HISTO *h, int b0, int b1) ;
HISTOGRAM *HISTOmakeCDF(HISTOGRAM *hsrc, HISTOGRAM *hdst) ;
int       HISTOfindBin(HISTOGRAM *h, float val) ;
int       HISTOfindBinWithCount(HISTOGRAM *h, float val) ;
HISTO     *HISTOclearBG(HISTOGRAM *hsrc, HISTOGRAM *hdst, int *pbg_end) ;
int       HISTOfindPreviousPeak(HISTOGRAM *h, int b0, int whalf) ;
HISTO     *HISTOlinearScale(HISTOGRAM *hsrc, HISTOGRAM *hdst, float scale,
                            float offset) ;

double    HISTOfindLinearFit(HISTOGRAM *h1,
                             HISTOGRAM *h2,
                             double amin, double amax,
                             double bmin, double bmax,
                             float *pa, float *pb) ;
float HISTOthreshSum(HISTOGRAM *h_mask, HISTOGRAM *h_src, float m_thresh) ;
HISTOGRAM *HISTOmakePDF(HISTO *h_src, HISTO *h_dst) ;

int HISTOcount(HISTO *h, double *samples, int nsamples);
int HISTOvalToBin(HISTO *h, double val);
HISTO *HISTObins(int nbins, double min, double max);
int HISTOvalToBinDirect(HISTOGRAM *histo, float val) ;
float HISTOvalToCount(HISTOGRAM *histo, float val) ;
int  HISTOwriteInto(HISTOGRAM *h, FILE *fp) ;
HISTOGRAM* HISTOreadFrom(FILE *fp) ;
double HISTOfindMedian(HISTOGRAM *h) ;
int    HISTOrobustGaussianFit(HISTOGRAM *h, double max_percentile, 
                              double *poffset, double *psigma) ;
HISTOGRAM *HISTOabs(HISTOGRAM *h, HISTOGRAM *habs) ;
HISTOGRAM *HISTOgaussianPDF(HISTOGRAM *h, double mean, double sigma, int nbins) ;
HISTOGRAM *HISTOgaussianCDF(HISTOGRAM *h, double mean, double sigma, int nbins) ;
double    HISTOgetCount(HISTOGRAM *h, float bin_val);
double    HISTOgetEntropy(HISTOGRAM *h);
HISTOGRAM *HISTOsoapBubbleZeros(HISTOGRAM *hsrc, HISTOGRAM *hdst, int niters) ;
int       HISTOfindMaxDerivative(HISTOGRAM *h, double min_count, double max_count, int whalf, 
                                 int grad_dir) ;
#endif
