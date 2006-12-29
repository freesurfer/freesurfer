/**
 * @file  mosaic.h
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


#ifndef _MOSAIC_H
#define _MOSAIC_H

int VolSS2MosSS(int cvol, int rvol, int svol,
                int ncvol, int nrvol,
                int ncmos, int nrmos,
                int *cmos, int *rmos,
                int *OutOfBounds);

int MosSS2VolSS(int cmos,  int rmos,
                int ncmos, int nrmos,
                int ncvol, int nrvol, int nsvol,
                int *cvol, int *rvol, int *svol,
                int *OutOfBounds);

int CheckMosaic(void);

#endif
