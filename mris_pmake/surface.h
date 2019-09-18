/**
 * @file  surface.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/27 21:18:07 $
 *    $Revision: 1.2 $
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

#ifndef __SURFACE_H__
#define __SURFACE_H__

#include "mri.h"
#include "mrisurf.h"

void  mark_geodesic(MRIS *surf, int vno_i, int vno_f, int mark);

#endif
