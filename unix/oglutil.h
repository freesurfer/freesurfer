#ifndef OGLUTIL_H
#define OGLUTIL_H

#include "mrisurf.h"

int    OGLUinit(MRI_SURFACE *mris, long frame_xsize, long frame_ysize) ;
int    OGLUcompile(MRI_SURFACE *mris, int *marked_vertices, int flags, 
                   float cslope) ;

#define SCALE_FACTOR       0.55f

extern double oglu_fov ;

#define PATCH_FLAG   0x0001
#define TP_FLAG      0x0002   /* show tangent plane and principal directions */
                              /*  of marked vertices */

#endif
