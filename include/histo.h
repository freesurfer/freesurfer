#ifndef HISTO_H
#define HISTO_H

#include <stdio.h>

#ifndef UCHAR
#define UCHAR  unsigned char
#endif

#define MAX_BINS  256
typedef struct
{
  int     nbins ;
  UCHAR   bins[MAX_BINS] ;   /* upper end of the range that maps to this bin */
  int     counts[MAX_BINS] ; /* # of voxels which map to this bin */
} HISTOGRAM, HISTO ;

int       HISTOfree(HISTOGRAM **phisto) ;
int       HISTOdump(HISTOGRAM *histo, FILE *fp) ;
HISTOGRAM *HISTOalloc(int nbins) ;
HISTOGRAM *HISTOcrunch(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOcopy(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOinvert(HISTOGRAM *histo_src, HISTOGRAM *histo_dst,int max_dst);
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, 
                          int max_out) ;
HISTOGRAM *HISTOclear(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOfillZeros(HISTOGRAM *histo_src, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOcompose(HISTOGRAM *histo1, HISTOGRAM *histo2, 
                        HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOnormalize(HISTOGRAM *histo_src, HISTOGRAM *histo_dst, 
                          int max_val) ;
HISTOGRAM *HISTOcomposeInvert(HISTOGRAM *histo_fwd, HISTOGRAM *histo_inv, 
                        HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOadd(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst) ;
HISTOGRAM *HISTOsubtract(HISTOGRAM *h1, HISTOGRAM *h2, HISTOGRAM *histo_dst) ;

#endif
