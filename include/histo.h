/*
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


#ifndef HISTO_H
#define HISTO_H

#include <stdio.h>

#define MAX_BINS  10000
typedef struct
{
  int     max_bins ;  // total number allocated - should never change
  int     nbins ;     // starts the same as max_bins, but can be less if the full range is not used
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
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,int max_out) ;
HISTOGRAM *HISTOclear(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
int       HISTOclearZeroBin(HISTOGRAM *h) ;
HISTOGRAM *HISTOclearCounts(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOfillZeros(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOcompose(HISTOGRAM *histo1, HISTOGRAM *histo2,
                        HISTOGRAM *histo_dst) ;
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
int       HISTOplot(HISTOGRAM *histo, const char *fname) ;
int       HISTOaddFractionalSample(HISTOGRAM *histo, float val, float bmin, float bmax, float frac);
int       HISTOaddSample(HISTOGRAM *histo, float val, float bmin, float bmax) ;
int       HISTOfindCurrentPeak(HISTOGRAM *histo,
                               int b0,
                               int wsize,
                               float min_pct) ;
int       HISTOfillHoles(HISTO *h) ;
float     HISTOtotal(HISTO *h) ;
float     HISTOtotalInRegion(HISTO *h, int b0, int b1) ;
HISTOGRAM *HISTOmakeReverseCDF(HISTOGRAM *hsrc, HISTOGRAM *hdst) ;
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
double    HISTOrmsDifference(HISTOGRAM *h1, HISTOGRAM *h2) ;
double    HISTOearthMoversDistance(HISTOGRAM *h1, HISTOGRAM *h2) ;
double    HISTOksDistance(HISTOGRAM *h1, HISTOGRAM *h2)  ;
HISTOGRAM    *HISTOeraseRightmostPeak(HISTOGRAM *hsrc, HISTOGRAM *hdst, int whalf, float min_pct, int min_val, int max_val) ;
int       HISTOerase(HISTOGRAM *h, int bmin, int bmax)  ;
int       HISTOisPeak(HISTOGRAM *h, int bin, int whalf) ;



typedef struct
{
  int     nbins1 ;
  int     nbins2 ;
  float   *bins1 ;   /* upper end of the range that maps to this bin */
  float   *bins2 ;   /* upper end of the range that maps to this bin */
  float   **counts ; /* # of voxels which map to this bin */
  float   bin_size1 ;
  float   bin_size2 ;
  float   min1 ;
  float   min2 ;
  float   max1 ;     // min and max vals in the histo
  float   max2 ;     // min and max vals in the histo
}
HISTOGRAM2D, HISTO2D ;

HISTOGRAM2D *HISTO2Dalloc(int nbins1, int nbins2) ;

int         HISTO2Dfree(HISTOGRAM2D **phisto) ;
int         HISTO2Ddump(HISTOGRAM2D *histo, FILE *fp) ;
int         HISTO2Dwrite(HISTOGRAM2D *h, char *fname) ;
int         HISTO2DwriteInto(HISTOGRAM2D *histo, FILE *fp) ;
HISTOGRAM2D *HISTO2DreadFrom(FILE *fp) ;
HISTOGRAM2D *HISTO2Dinit(HISTOGRAM2D *h, int nbins1, int nbins2, double mn1, int mx1, double mn2, double mx2) ;
HISTOGRAM2D *HISTO2Drealloc(HISTOGRAM2D *histo, int nbins1, int nbins2) ;
int         HISTO2DaddFractionalSample(HISTOGRAM2D *histo, float val1, float val2, float bmin1, float bmax1, float bmin2, float bmax2, float frac) ;

int         HISTO2DaddSample(HISTOGRAM2D *histo, float val1, float val2, float bmin1, float bmax1, float bmin2, float bmax2) ;
int         HISTO2Dplot(HISTOGRAM2D *histo, char *fname) ;
HISTOGRAM2D *HISTO2DmakePDF(HISTO2D *h_src, HISTO2D *h_dst) ;
HISTOGRAM2D *HISTO2Dclear(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst) ;
HISTOGRAM2D *HISTO2Dcopy(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst) ;
double       HISTO2DgetCount(HISTOGRAM2D *h, float bin_val1, float bin_val2);
int          HISTO2DfindBin1(HISTOGRAM2D *h, float val) ;
int          HISTO2DfindBin2(HISTOGRAM2D *h, float val) ;
HISTOGRAM2D *HISTO2Dsmooth(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst,float sigma) ;
HISTOGRAM2D *HISTO2DsmoothAnisotropic(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst,float sigma1, float sigma2) ;
HISTOGRAM2D *HISTO2DsmoothBins1(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst,float sigma) ;
HISTOGRAM2D *HISTO2DsmoothBins2(HISTOGRAM2D *histo_src, HISTOGRAM2D *histo_dst,float sigma) ;
HISTOGRAM2D *HISTO2Dread(char *fname) ;
HISTOGRAM2D *HISTO2DsoapBubbleZeros(HISTOGRAM2D *hsrc, HISTOGRAM2D *hdst, int niters) ;
float       HISTOcomputeFWHM(HISTOGRAM *h, int peak) ;
int HISTOwriteTxt(HISTOGRAM *histo, const char *fname) ;
int HISTOsumNorm(HISTOGRAM *histo);
HISTOGRAM *HISTOcumsum(HISTOGRAM *h, HISTOGRAM *hout);

#endif
