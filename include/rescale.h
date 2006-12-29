/**
 * @file  rescale.h
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
