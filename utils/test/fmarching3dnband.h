/*****************************************************/
/*          File: fmarching3d.h                      */
/*   Purpose: Header file for fmarching3d functions  */
/*         Add NarrowBand Threshold                  */

#ifndef FMARCHING3D
#define FMARCHING3D

#include <math.h>
#include "heap.h"
#include "mri.h"
#define ALIVE 1
#define NBAND 2
#define FAWAY 3 


float ReCompute(float,float,float,float,float,float,unsigned char,unsigned char,
         unsigned char,unsigned char,unsigned char,unsigned char);
void fmarching3d(MRI *, MRI *,  float Thred);

#endif


