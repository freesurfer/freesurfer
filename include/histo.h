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
} HISTOGRAM, HISTO ;

int       HISTOfree(HISTOGRAM **phisto) ;
int       HISTOdump(HISTOGRAM *histo, FILE *fp) ;
int       HISTOwrite(HISTOGRAM *histo, char *fname) ;
HISTOGRAM *HISTOalloc(int nbins) ;
HISTOGRAM *HISTOrealloc(HISTOGRAM *histo, int nbins) ;
HISTOGRAM *HISTOcrunch(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
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
int       HISTOfindCurrentPeak(HISTOGRAM *histo, int b0, int wsize, float min_pct) ;
int       HISTOfillHoles(HISTO *h) ;
int       HISTOtotal(HISTO *h) ;
int       HISTOtotalInRegion(HISTO *h, int b0, int b1) ;
int       HISTOfindBin(HISTOGRAM *h, float val) ;
HISTO     *HISTOclearBG(HISTOGRAM *hsrc, HISTOGRAM *hdst, int *pbg_end) ;
int       HISTOfindPreviousPeak(HISTOGRAM *h, int b0, int whalf) ;

#endif
