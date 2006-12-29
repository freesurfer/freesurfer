/**
 * @file  fmarching3dnband.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:44 $
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


/*****************************************************/
/*          File: fmarching3d.h                      */
/*   Purpose: Header file for fmarching3d functions  */
/*         Add NarrowBand Threshold                  */

#ifndef FMARCHING3D
#define FMARCHING3D

#include <math.h>
#include "heap_vol.h"
#include "mri.h"
#define ALIVE 1
#define NBAND 2
#define FAWAY 3


float ReCompute(float,float,float,float,float,float,unsigned char,unsigned char,
                unsigned char,unsigned char,unsigned char,unsigned char);
void fmarching3d(MRI *, MRI *,  float Thred);

#endif


