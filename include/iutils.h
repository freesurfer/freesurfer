/**
 * @file  iutils.h
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


#ifndef IUTILS_H
#define IUTILS_H

#if 1
/*
  these prototypes cause the canny edge detector to fail - I don't know why,
  and I'm not going to spend the time to find out.
*/
void crop_image(unsigned char *imageptr, int *cols,int *rows, int cropcornerax,
                int cropcorneray, int cropcornerbx, int cropcornerby) ;
#endif
void get_histogram_threshold(int hgram[], int histsize, int pixelmax,
                             int pixelmin, float fraction,
                             int zflag, float *ht, float *lt) ;
void histogram(short *theimage, int xsize, int ysize, int pixelmax,
               int pixelmin, int hgram[], int histsize) ;

#endif
