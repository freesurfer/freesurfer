/*******************************************************************
 * Name:    rescale.h
 * Author:  Douglas N. Greve, 5/14/96
 * Purpose: rescales the pixels of an image to within the given
 *          maximum and minimum.  As a bonus, it also returns
 *          the row and column of the min and max. 
 *******************************************************************/
#ifndef RESCALE_INC
#define RESCALE_INC
int h_rescale(struct header *phdSrc, 
        float NewMin, float NewMax,
        int MinPoint[2], int MaxPoint[2],
        struct header *phdDst);
int h_rescale_fb(struct header *phdSrc, 
     float Slope, float Offset,
     struct header *phdDst);
int h_rescale_ff(struct header *phdSrc, 
     float Slope, float Offset,
     struct header *phdDst);
#endif /** from #ifndef RESCALE_INC *****/
