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
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.4 $
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


