/**************************************************************
 * Name:    minmaxrc.h
 * Author:  Douglas N. Greve, 5/14/96
 * Purpose: hips-based routines to find the minimum and maximum
 *          pixel values in an image as well as their row and
 *          column.
 ***************************************************************/
#ifndef MINMAXRC_INC
#define MINMAXRC_INC
int h_minmaxrc(struct header *phdSrc, 
       Pixelval *Min, int MinPoint[2],
       Pixelval *Max, int MaxPoint[2]);
int h_minmaxrc_b(struct header *phdSrc, 
         byte *Min, int MinPoint[2],
         byte *Max, int MaxPoint[2]);
int h_minmaxrc_s(struct header *phdSrc, 
         short *Min, int MinPoint[2],
         short *Max, int MaxPoint[2]);
int h_minmaxrc_i(struct header *phdSrc, 
         int *Min, int MinPoint[2],
         int *Max, int MaxPoint[2]);
int h_minmaxrc_f(struct header *phdSrc, 
         float *Min, int MinPoint[2],
         float *Max, int MaxPoint[2]);
#endif /* from #ifndef MINMAXRC_INC */
