/**
 * @file  retinotopy.h
 * @brief Utilities for retinotopy analysis.
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Doug Greve (and Marty and Anders, for now)
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2008/01/15 23:47:39 $
 *    $Revision: 1.3 $
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



#ifndef RETINOTOPY_INC
#define RETINOTOPY_INC

#include "mri.h"
#include "mrisurf.h"

#ifdef X
#undef X
#endif

void RETcompute_fieldsign(MRIS *mris);
void RETcompute_angles(MRIS *mris);
int RETlogMap(MRIS *surf, double k, double a, double xc0, double yc0);
int RETinvLogMapFunc(double xc, double yc, double xc0, double yc0, 
		     double a, double k, double *r, double *theta);
float RETcircsubtract(float a,float b) ;
int RETreverseSign(MRI *mri);

#endif
