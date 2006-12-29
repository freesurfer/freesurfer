/**
 * @file  minmaxrc.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
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
