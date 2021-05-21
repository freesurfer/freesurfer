/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


