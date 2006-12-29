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
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
