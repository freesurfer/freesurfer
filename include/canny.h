/**
 * @file  canny.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.3 $
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


/*****************************************************
 * Name:    canny.h
 *****************************************************/
#ifndef CANNY_H

#define NOEDGE 0
#define POSSIBLE_EDGE 254
#define EDGE 255

/*----------------*/
/* CONVOLVE IMAGE */
/*--------------- */

/* edge detection parameters */
#define DEFAULT_SIGMA        1.0
#define DEFAULT_WINDOWSIZE   9
#define DEFAULT_LFRAC        0.5
#define DEFAULT_HFRAC        0.9

/* POSSIBLE VALUES OF THE FILTERTYPE PARMETER */

#define NOSYMMETRY      0
#define SYMMETRIC       1
#define ANITSYMMETRIC   2

/* POSSIBLE VALUES OF THE BOUNDERY PARAMETER */

#define ZERO      0    /* ZEROES THE REGION OUTSIDE OF THE IMAGE */
#define WRAP      1    /* THE FILTER IS CIRCULAR */
#define MAKESMALL 2    /* RETURNS A SMALLER IMAGE. */
#define EXTEND    3    /* EXTENDS THE BOUNDARY PIXELS */
#define MASKIT    4    /* MASKS THE BOUNDARY TO ZERO */

/* POSSIBLE VALUES OF THE STATUS PARAMETER */

#define SUCCESS        1    /* SUCCESSFULL EXECUTION */
#define NO_SUCH_FILTER 2  /* IF THE FILTER TYPE PARAMETER IS NOT ONE OF THE ALLOWED VALS. */
#define EVENWINDOWSIZE  3  /* IF THE FILTER IS EVENSIZED */
#define NOSUCHDIRECTION 4 /* Direction is not XDIR of YDIR */
#define NOSUCHBOUNDERY 5  /* Nonexistant boundery option specified */

/* POSSIBLE VALUES OF THE DIRECTION PARAMETER */

#define XDIR  1
#define YDIR  2

#include <hips.h> /* dng */

/* prototypes */
void canny(int *magmax, int *hthresh, int *lthresh, int *image, int *xsize, int *ysize, short *shortim,
           int *windowsize, double *sigma, int *bordermode, double *hfrac, double *lfrac, int *pflag,
           short *gx, short *gy, short *mag, int *hist, int *histsize, unsigned char *nms,
           unsigned char *edgemap, float *gm, float *gmp,short *temp) ;
void cleanup(unsigned char *map, int xsize, int ysize) ;
void find_edges(unsigned char *map, short *mag, int xsize, int ysize, int maxmag, float hpixel_fraction,
                float lpixel_fraction, int *hgram, int hsize, int *actual_hthresh, int *actual_lthresh);
void follow_edges(unsigned char *edgemapptr, short *edgemagptr) ;
void clear_borders(unsigned char *charimage, int xsize, int ysize) ;
void gauss_filter(short *inimage, int inx, int iny, int direction, int boundary, int masksize,
                  float sigma, short *grad, int *outx, int *outy,
                  float *gmask, float *gprimemask, short *tempimage) ;
void copyimage(int *charimage, int ncols, int nrows, short *shortimage) ;
void thin(unsigned char *edges, int height, int width) ;

int h_canny(struct header *Isrc, struct header *Idst, double sigma,
            int mask_size,double lfrac,double hfrac,int dothin); /* dng */

#endif


