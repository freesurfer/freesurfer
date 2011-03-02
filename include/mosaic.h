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
 *    $Date: 2011/03/02 00:04:10 $
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
