#ifndef HISTO_H
#define HISTO_H

#include <stdio.h>

#ifndef UCHAR
#define UCHAR  unsigned char
#endif

#define MAX_BINS  1024
typedef struct
{
  int     nbins ;
  float   bins[MAX_BINS] ;   /* upper end of the range that maps to this bin */
  float   counts[MAX_BINS] ; /* # of voxels which map to this bin */
} HISTOGRAM, HISTO ;

int       HISTOfree(HISTOGRAM **phisto) ;
int       HISTOdump(HISTOGRAM *histo, FILE *fp) ;
int       HISTOwrite(HISTOGRAM *histo, char *fname) ;
HISTOGRAM *HISTOalloc(int nbins) ;
HISTOGRAM *HISTOcrunch(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOcopy(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOinvert(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,int max_dst);
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, 
                          int max_out) ;
HISTOGRAM *HISTOclear(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
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
int       HISTOfindValley(HISTOGRAM *h, int wsize, int b0, int b1) ;
int       HISTOfindLastPeakInRegion(HISTOGRAM *h, int wsize, float min_pct, 
                                    int b0, int b1) ;
int       HISTOcountPeaksInRegion(HISTOGRAM *h, int wsize, float min_pct, 
                                  int *peaks, int max_peaks, int b0, int b1) ;
int       HISTOfindFirstPeakInRegion(HISTOGRAM *h, int wsize, float min_pct, 
                                     int b0, int b1) ;
int       HISTOfindHighestPeakInRegion(HISTOGRAM *h, int b0, int b1);
int       HISTOplot(HISTOGRAM *histo, char *fname) ;
int       HISTOaddSample(HISTOGRAM *histo, float val, float bmin, float bmax) ;

#endif
