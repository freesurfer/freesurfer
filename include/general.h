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


#ifndef GENERAL_H
#define GENERAL_H

void order(int *small, int *big) ;

#if 1
/* putting these prototypes in causes the canny edge detecter to fail, so
   leave them out for now */
short nearestshort(float x) ;
void  fporder(float *big, float *small) ;
float fpclip(float x, float bound1, float bound2) ;
int nearestint(float x) ;
#endif

#endif
